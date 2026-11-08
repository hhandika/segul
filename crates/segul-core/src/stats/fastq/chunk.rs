//! Read FASTQ files in large blocks that end on a record boundary.
//!
//! Each block holds whole 4-line records, so workers can parse blocks
//! independently. The compression format is detected from the first bytes
//! of the file, not its extension.
use std::fs::File;
use std::io::{BufReader, Read};
use std::num::NonZero;
use std::path::Path;
use std::sync::mpsc::{Receiver, Sender, SyncSender};

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use noodles::bgzf;

/// Target size of a block of records for a file read by many threads.
pub(crate) const CHUNK_SIZE: usize = 8 * 1024 * 1024;
/// Target size of a block when many files are read at once, each by
/// fewer threads. Small enough that blocks in flight for every file stay
/// within a few hundred megabytes.
pub(crate) const MIN_CHUNK_SIZE: usize = 1024 * 1024;
// Threads per file that get the largest blocks.
const FULL_CHUNK_THREADS: usize = 8;

/// Block size for a file read by `threads` threads.
pub(crate) fn chunk_size_for(threads: usize) -> usize {
    (CHUNK_SIZE / FULL_CHUNK_THREADS * threads).clamp(MIN_CHUNK_SIZE, CHUNK_SIZE)
}

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
// Header length needed to find the BGZF extra field.
const BGZF_HEADER_LEN: usize = 14;
const FLAG_EXTRA: u8 = 0x04;
const LINE_FEED: u8 = b'\n';
const LINES_PER_RECORD: usize = 4;

/// Compression of a FASTQ file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Compression {
    Plain,
    Gzip,
    /// Blocked gzip, which can be decompressed in parallel.
    Bgzf,
}

impl Compression {
    /// Detect the compression from the magic bytes at the start of a file.
    pub(crate) fn detect(path: &Path) -> Result<Self> {
        let mut file =
            File::open(path).with_context(|| format!("Failed opening {}", path.display()))?;
        let mut header = [0u8; BGZF_HEADER_LEN];
        let mut len = 0;
        while len < header.len() {
            match file.read(&mut header[len..])? {
                0 => break,
                n => len += n,
            }
        }
        Ok(Self::from_header(&header[..len]))
    }

    fn from_header(header: &[u8]) -> Self {
        if header.len() < 2 || header[..2] != GZIP_MAGIC {
            return Self::Plain;
        }
        let is_bgzf = header.len() >= BGZF_HEADER_LEN
            && header[3] & FLAG_EXTRA != 0
            && &header[12..14] == b"BC";
        if is_bgzf { Self::Bgzf } else { Self::Gzip }
    }
}

/// Open a FASTQ file as a stream of decompressed bytes.
///
/// `decode_threads` is the number of threads used for BGZF files.
pub(crate) fn open_fastq(path: &Path, decode_threads: usize) -> Result<Box<dyn Read + Send>> {
    let compression = Compression::detect(path)?;
    let file = File::open(path).with_context(|| format!("Failed opening {}", path.display()))?;
    let reader: Box<dyn Read + Send> = match compression {
        Compression::Plain => Box::new(file),
        Compression::Gzip => Box::new(MultiGzDecoder::new(BufReader::new(file))),
        Compression::Bgzf => {
            let workers = NonZero::new(decode_threads).unwrap_or(NonZero::<usize>::MIN);
            Box::new(bgzf::io::MultithreadedReader::with_worker_count(
                workers, file,
            ))
        }
    };
    Ok(reader)
}

/// A block of whole FASTQ records.
pub(crate) struct Chunk {
    /// Index of the first read in the block, counted from the start of the file.
    pub(crate) first_read: u64,
    pub(crate) data: Vec<u8>,
}

/// Splits a byte stream into blocks of whole records.
pub(crate) struct ChunkReader<R> {
    reader: R,
    chunk_size: usize,
    /// Bytes after the last whole record of the previous block.
    carry: Vec<u8>,
    reads: u64,
    eof: bool,
}

impl<R: Read> ChunkReader<R> {
    pub(crate) fn new(reader: R, chunk_size: usize) -> Self {
        Self {
            reader,
            chunk_size,
            carry: Vec::new(),
            reads: 0,
            eof: false,
        }
    }

    /// Fill `buffer` with the next block. Returns `None` at the end of the stream.
    pub(crate) fn next_chunk(&mut self, mut buffer: Vec<u8>) -> Result<Option<Chunk>> {
        buffer.clear();
        buffer.append(&mut self.carry);
        let mut target = self.chunk_size;
        loop {
            if !self.eof {
                self.fill(&mut buffer, target)?;
            }
            if self.eof {
                if buffer.is_empty() {
                    return Ok(None);
                }
                return Ok(Some(self.take(buffer, None)));
            }
            // A record longer than the block keeps the block growing.
            if let Some(cut) = record_boundary(&buffer) {
                return Ok(Some(self.take(buffer, Some(cut))));
            }
            target = buffer.len() + self.chunk_size;
        }
    }

    /// Send blocks to `sender` until the stream ends, reusing returned buffers.
    pub(crate) fn send_all(
        &mut self,
        sender: SyncSender<Chunk>,
        recycled: Receiver<Vec<u8>>,
    ) -> Result<()> {
        loop {
            let buffer = recycled
                .try_recv()
                .unwrap_or_else(|_| Vec::with_capacity(self.chunk_size));
            match self.next_chunk(buffer)? {
                Some(chunk) => {
                    // Workers stop receiving after an error.
                    if sender.send(chunk).is_err() {
                        return Ok(());
                    }
                }
                None => return Ok(()),
            }
        }
    }

    /// Read until `buffer` holds `target` bytes or the stream ends.
    fn fill(&mut self, buffer: &mut Vec<u8>, target: usize) -> Result<()> {
        // Decoders return small reads, so zero the space once, not per read.
        let mut filled = buffer.len();
        buffer.resize(target.max(filled), 0);
        while filled < target {
            match self.reader.read(&mut buffer[filled..]) {
                Ok(0) => {
                    self.eof = true;
                    break;
                }
                Ok(n) => filled += n,
                Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(error) => return Err(error).context("Failed reading FASTQ data"),
            }
        }
        buffer.truncate(filled);
        Ok(())
    }

    fn take(&mut self, mut buffer: Vec<u8>, cut: Option<usize>) -> Chunk {
        if let Some(cut) = cut {
            self.carry.extend_from_slice(&buffer[cut..]);
            buffer.truncate(cut);
        }
        let first_read = self.reads;
        self.reads += count_records(&buffer);
        Chunk {
            first_read,
            data: buffer,
        }
    }
}

/// Position just after the last line feed that ends a whole record.
fn record_boundary(data: &[u8]) -> Option<usize> {
    let mut boundary = None;
    memchr::memchr_iter(LINE_FEED, data)
        .enumerate()
        .filter(|(line, _)| (line + 1) % LINES_PER_RECORD == 0)
        .for_each(|(_, pos)| boundary = Some(pos + 1));
    boundary
}

/// Number of records in a block, counting a last line without a line feed.
pub(crate) fn count_records(data: &[u8]) -> u64 {
    let mut lines = bytecount::count(data, LINE_FEED);
    if data.last().is_some_and(|b| *b != LINE_FEED) {
        lines += 1;
    }
    lines.div_ceil(LINES_PER_RECORD) as u64
}

/// Return a processed buffer to the reader for reuse.
pub(crate) fn recycle(sender: &Sender<Vec<u8>>, buffer: Vec<u8>) {
    // The reader may have finished already.
    let _ = sender.send(buffer);
}

#[cfg(test)]
mod test {
    use super::*;

    const RECORDS: &[u8] = b"@r1\nACGT\n+\nIIII\n@r2\nAC\n+\nII\n@r3\nACG\n+\nIII\n";

    #[test]
    fn test_detect_compression() {
        assert_eq!(Compression::from_header(b"@r1"), Compression::Plain);
        assert_eq!(
            Compression::from_header(&[0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0]),
            Compression::Gzip
        );
        let bgzf = [0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C'];
        assert_eq!(Compression::from_header(&bgzf), Compression::Bgzf);
    }

    #[test]
    fn test_chunks_end_on_records() {
        // A tiny block size forces records to span reads.
        let mut reader = ChunkReader::new(RECORDS, 5);
        let mut chunks = Vec::new();
        while let Some(chunk) = reader.next_chunk(Vec::new()).unwrap() {
            chunks.push(chunk);
        }
        let joined: Vec<u8> = chunks.iter().flat_map(|c| c.data.clone()).collect();
        assert_eq!(joined, RECORDS);
        assert!(chunks.iter().all(|c| c.data.starts_with(b"@")));
        let firsts: Vec<u64> = chunks.iter().map(|c| c.first_read).collect();
        assert_eq!(firsts, vec![0, 1, 2]);
    }

    #[test]
    fn test_count_records_without_final_line_feed() {
        assert_eq!(count_records(b"@r1\nACGT\n+\nIIII"), 1);
        assert_eq!(count_records(RECORDS), 3);
        assert_eq!(count_records(b""), 0);
    }
}
