//! A simple fasta parser that reads a fasta file and returns a sequence matrix.

use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use std::path::Path;

use indexmap::IndexSet;

use crate::helper::alphabet;
use crate::helper::sequence::SeqCheck;
use crate::helper::types::{DataType, Header, SeqMatrix};
use crate::parser;

/// Parse a fasta file sequence IDs.
/// return a set of sequence IDs.
pub fn parse_only_id(input: &Path) -> IndexSet<String> {
    let file = File::open(input).expect("Failed opening a fasta file.");
    let mut buff = BufReader::with_capacity(parser::IO_BUF_SIZE, file);
    let mut ids = IndexSet::new();
    // Reuse the line buffer and only decode header lines.
    let mut line = Vec::new();
    while matches!(buff.read_until(b'\n', &mut line), Ok(n) if n > 0) {
        if let Some(id) = line.strip_prefix(b">") {
            match std::str::from_utf8(id) {
                Ok(id) => ids.insert(id.trim().to_string()),
                Err(_) => break,
            };
        }
        line.clear();
    }
    ids
}

/// Parse a fasta file's unique IDs and the length of the first sequence.
#[cfg(feature = "align")]
pub(crate) fn parse_id_and_nchar(input: &Path) -> (IndexSet<String>, usize) {
    let file = File::open(input).expect("Failed opening a fasta file.");
    let mut buff = BufReader::with_capacity(parser::IO_BUF_SIZE, file);
    let mut ids = IndexSet::new();
    let mut line = Vec::new();
    let mut first_len = 0usize;
    let mut in_first = false;
    while matches!(buff.read_until(b'\n', &mut line), Ok(n) if n > 0) {
        if let Some(id) = line.strip_prefix(b">") {
            in_first = ids.is_empty();
            match std::str::from_utf8(id) {
                Ok(id) => {
                    ids.insert(id.trim().to_string());
                }
                Err(_) => break,
            };
        } else if in_first {
            first_len += line.trim_ascii().len();
        }
        line.clear();
    }
    (ids, first_len)
}

/// Stream records from an already-buffered FASTA reader without per-record allocations.
pub fn for_each_fasta_record<R: BufRead, F: FnMut(&[u8], &[u8])>(reader: &mut R, mut f: F) {
    let mut line = Vec::new();
    let mut id = Vec::new();
    let mut seq = Vec::new();
    while matches!(reader.read_until(b'\n', &mut line), Ok(n) if n > 0) {
        let trimmed = line.trim_ascii_end();
        if let Some(hdr) = trimmed.strip_prefix(b">") {
            if !id.is_empty() {
                f(&id, &seq);
                id.clear();
                seq.clear();
            }
            id.extend_from_slice(hdr.trim_ascii());
        } else {
            seq.extend_from_slice(trimmed.trim_ascii_start());
        }
        line.clear();
    }
    if !id.is_empty() {
        f(&id, &seq);
    }
}

/// Parse a fasta file and return a sequence matrix.
/// It does not allow duplicate IDs or spaces in the sequence IDs.
pub struct Fasta<'a> {
    input: &'a Path,
    datatype: &'a DataType,
    pub matrix: SeqMatrix,
    pub header: Header,
}

impl<'a> Fasta<'a> {
    /// Create a new Fasta instance.
    pub fn new(input: &'a Path, datatype: &'a DataType) -> Self {
        Self {
            input,
            datatype,
            matrix: SeqMatrix::new(),
            header: Header::new(),
        }
    }

    /// Parse a fasta file.
    pub fn parse(&mut self) {
        let file = File::open(self.input).expect("Failed opening a fasta file");
        let file_len = file.metadata().map(|m| m.len() as usize).unwrap_or(0);
        self.matrix.reserve_bytes(file_len.min(parser::MAX_RESERVE));
        let mut buff = BufReader::with_capacity(parser::IO_BUF_SIZE, file);
        self.parse_matrix(&mut buff);
        let mut seq_info = SeqCheck::new();
        assert!(
            !self.matrix.is_empty(),
            "{} is empty. \
        Make sure the file format is fasta or it is not an empty file!",
            self.input.display()
        );
        seq_info.check(&self.matrix);
        self.header.aligned = seq_info.is_alignment;
        self.header.nchar = seq_info.longest;
        self.header.ntax = self.matrix.len();
        self.match_header_datatype();
    }

    fn parse_matrix<R: BufRead>(&mut self, buff: &mut R) {
        let mut line = Vec::new();
        let mut current_id = String::new();
        while matches!(buff.read_until(b'\n', &mut line), Ok(n) if n > 0) {
            let trimmed = line.trim_ascii_end();
            if let Some(hdr) = trimmed.strip_prefix(b">") {
                if !current_id.is_empty() {
                    self.commit_record(std::mem::take(&mut current_id));
                }
                let id_str = std::str::from_utf8(hdr.trim_ascii()).unwrap_or_else(|_| {
                    panic!("Failed reading FASTA header in {}", self.input.display())
                });
                current_id.push_str(id_str);
            } else if !current_id.is_empty() {
                self.matrix
                    .data_mut()
                    .extend_from_slice(trimmed.trim_ascii_start());
            }
            line.clear();
        }
        if !current_id.is_empty() {
            self.commit_record(current_id);
        }
    }

    fn commit_record(&mut self, id: String) {
        if let Some(original_seq) = self.matrix.get(&id) {
            let same = original_seq == self.matrix.staged();
            panic!(
                "\nFound duplicate IDs for file {}. \
                First duplicate ID found: {}. \
                Both sequences are the same: {}.\n",
                self.input.display(),
                id,
                same
            );
        }
        alphabet::check_valid_seq_bytes(self.input, self.datatype, &id, self.matrix.staged());
        self.matrix.commit_staged(id);
    }

    #[inline]
    fn match_header_datatype(&mut self) {
        if let DataType::Aa = self.datatype {
            self.header.datatype = String::from("protein")
        };
    }
}

pub struct Records {
    pub id: String,
    pub seq: Vec<u8>,
}

impl Records {
    pub fn new(id: &str, seq: &[u8]) -> Self {
        Self {
            id: String::from(id),
            seq: seq.to_vec(),
        }
    }
}

pub struct FastaReader<R> {
    reader: BufReader<R>,
    id: String,
    seq: Vec<u8>,
    line: Vec<u8>,
}

impl<R: Read> FastaReader<R> {
    pub fn new(file: R) -> Self {
        Self {
            reader: BufReader::with_capacity(parser::IO_BUF_SIZE, file),
            id: String::new(),
            seq: Vec::new(),
            line: Vec::new(),
        }
    }

    /// Stream every record without allocating a `Records` struct per sequence.
    pub fn for_each_record<F: FnMut(&[u8], &[u8])>(&mut self, f: F) {
        for_each_fasta_record(&mut self.reader, f);
    }

    pub fn next_seq(&mut self) -> Option<Records> {
        // Reuse the line buffer instead of allocating a new buffer per line.
        loop {
            self.line.clear();
            match self.reader.read_until(b'\n', &mut self.line) {
                Ok(n) if n > 0 => (),
                _ => break,
            }
            let line = self.line.trim_ascii_end();
            if let Some(id_bytes) = line.strip_prefix(b">") {
                let id_str = std::str::from_utf8(id_bytes.trim_ascii()).ok()?;
                if self.id.is_empty() {
                    self.id.push_str(id_str);
                    self.seq.clear();
                } else {
                    let prev_id = std::mem::replace(&mut self.id, id_str.to_string());
                    let prev_seq = std::mem::take(&mut self.seq);
                    return Some(Records {
                        id: prev_id,
                        seq: prev_seq,
                    });
                }
            } else {
                self.seq.extend_from_slice(line.trim_ascii_start());
            }
        }

        // Return the last found record.
        if !self.id.is_empty() {
            let id = std::mem::take(&mut self.id);
            let seq = std::mem::take(&mut self.seq);
            Some(Records { id, seq })
        } else {
            None
        }
    }
}

impl<R: Read> Iterator for FastaReader<R> {
    type Item = Records;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_seq()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    const DNA: DataType = DataType::Dna;
    const AA: DataType = DataType::Aa;

    #[test]
    fn read_fasta_simple_test() {
        let path = Path::new("tests/files/simple.fas");
        let mut fasta = Fasta::new(path, &DNA);
        fasta.parse();

        assert_eq!(2, fasta.matrix.len());
    }

    #[test]
    fn check_is_alignment_test() {
        let path = Path::new("tests/files/simple.fas");
        let mut fasta = Fasta::new(path, &DNA);
        fasta.parse();

        assert!(fasta.header.aligned);
    }

    #[test]
    fn check_isnot_alignment_test() {
        let path = Path::new("tests/files/unaligned.fas");
        let mut fasta = Fasta::new(path, &DNA);
        fasta.parse();

        assert!(!fasta.header.aligned);
    }

    #[test]
    fn interleaved_fas_test() {
        let path = Path::new("tests/files/interleave.fas");
        let file = File::open(path).unwrap();
        let rec = FastaReader::new(file);
        let mut seq = SeqMatrix::new();
        rec.into_iter().for_each(|r| {
            seq.insert(r.id, &r.seq);
        });

        let res = b"AGTATGATGTATATGTAT".as_slice();
        let res_2 = b"AGTATGATGTATAAAAAA".as_slice();
        assert_eq!(Some(res), seq.get("ABCD"));
        assert_eq!(Some(res_2), seq.get("ABCE"));
    }

    #[test]
    fn simple_aa_fas_test() {
        let sample = Path::new("tests/files/simple_aa.fas");
        let mut fas = Fasta::new(sample, &AA);
        fas.parse();
        let res = b"MAYPMQLGFQDATSPI".as_slice();
        assert_eq!(Some(res), fas.matrix.get("ABCE"));
        assert_eq!(String::from("protein"), fas.header.datatype);
    }

    #[test]
    fn simple_trailing_space_fas_test() {
        let sample = Path::new("tests/files/simple_id_whitespaces.fas");
        let mut fas = Fasta::new(sample, &DNA);
        fas.parse();
        let res = b"ATGTAT".as_slice();
        assert_eq!(Some(res), fas.matrix.get("CDEF"));
    }
}
