//! FASTQ read statistics.
//!
//! Files are read in large blocks of whole records on one thread while
//! rayon workers parse and count the blocks. Counts from each worker are
//! merged at the end, so a single large file uses every core.
//!
//! Plain, gzip, and BGZF files are detected from their first bytes.
//! BGZF files are also decompressed in parallel.
mod accumulate;
mod chunk;
mod parse;
pub mod report;
mod sketch;

use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::mpsc::{channel, sync_channel};
use std::thread;

use anyhow::{Context, Result};
use rayon::prelude::*;

use crate::helper::types::{PhredOffset, SeqReadFmt};
use crate::stats::common::StreamStats;
use crate::stats::qscores::ReadQScore;
use crate::stats::read::{ReadRecord, ReadSummary};

use accumulate::{ReadAccumulator, detect_offset};
use chunk::{ChunkReader, chunk_size_for, open_fastq, recycle};
use report::{ReadMetrics, SampleReadReport};
pub use sketch::DEFAULT_DUP_SKETCH_SIZE;
use sketch::DuplicationSketch;

// Blocks waiting for a worker, per worker thread.
const QUEUED_CHUNKS_PER_THREAD: usize = 2;

/// Computes read statistics for one FASTQ file.
pub struct FastqAnalyzer<'a> {
    path: &'a Path,
    phred: PhredOffset,
    /// Adds per-position, adapter, and duplication statistics.
    complete: bool,
    /// Threads available for this file. Used to decompress BGZF files
    /// and to bound the blocks held in memory.
    threads: usize,
    /// Distinct sequences kept to estimate duplication. Zero keeps all.
    dup_sketch_size: usize,
    chunk_size: usize,
}

impl<'a> FastqAnalyzer<'a> {
    pub fn new(path: &'a Path, phred: PhredOffset, complete: bool, threads: usize) -> Self {
        Self {
            path,
            phred,
            complete,
            threads,
            dup_sketch_size: DEFAULT_DUP_SKETCH_SIZE,
            chunk_size: chunk_size_for(threads),
        }
    }

    /// Keep this many distinct sequences to estimate duplication, or zero
    /// to count every distinct sequence. Each kept sequence uses about
    /// 30 bytes, shared by every thread reading the file.
    pub fn with_dup_sketch_size(mut self, size: usize) -> Self {
        self.dup_sketch_size = size;
        self
    }

    /// Read the file once and return its statistics.
    pub fn analyze(&self) -> Result<SampleReadReport> {
        let duplicates = self
            .complete
            .then(|| DuplicationSketch::new(self.dup_sketch_size));
        let acc = self
            .accumulate(duplicates.as_ref())
            .with_context(|| format!("Failed summarizing {}", self.path.display()))?;
        Ok(SampleReadReport::from_accumulator(
            self.path,
            &acc,
            duplicates.as_ref(),
        ))
    }

    fn accumulate(&self, duplicates: Option<&DuplicationSketch>) -> Result<ReadAccumulator> {
        let reader = open_fastq(self.path, self.threads)?;
        let mut chunks = ChunkReader::new(reader, self.chunk_size);
        let Some(first) = chunks.next_chunk(Vec::with_capacity(self.chunk_size))? else {
            return Ok(ReadAccumulator::new(self.complete, self.default_offset()));
        };
        let offset = match self.phred {
            PhredOffset::Auto => detect_offset(&first.data),
            PhredOffset::Phred33 => 33,
            PhredOffset::Phred64 => 64,
        };
        let complete = self.complete;
        // Files summarized together share the memory budget.
        let capacity = self.threads.max(1) * QUEUED_CHUNKS_PER_THREAD;
        let (sender, receiver) = sync_channel(capacity);
        let (recycle_sender, recycled) = channel();
        thread::scope(|scope| {
            let reader = scope.spawn(move || -> Result<()> {
                if sender.send(first).is_err() {
                    return Ok(());
                }
                chunks.send_all(sender, recycled)
            });
            let counts = receiver
                .into_iter()
                .par_bridge()
                .try_fold(
                    || ReadAccumulator::new(complete, offset),
                    |mut acc, chunk| {
                        acc.add_chunk(&chunk.data, chunk.first_read, duplicates)?;
                        recycle(&recycle_sender, chunk.data);
                        Ok::<_, anyhow::Error>(acc)
                    },
                )
                .try_reduce(
                    || ReadAccumulator::new(complete, offset),
                    |a, b| Ok(a.merge(b)),
                );
            let read_result = reader
                .join()
                .expect("Failed joining the FASTQ reader thread");
            let mut acc = counts?;
            read_result?;
            acc.finish();
            Ok(acc)
        })
    }

    fn default_offset(&self) -> u8 {
        match self.phred {
            PhredOffset::Phred64 => 64,
            _ => 33,
        }
    }
}

/// Count the reads in a FASTQ file without parsing them.
pub fn count_reads(path: &Path, decode_threads: usize) -> Result<usize> {
    let mut reader = open_fastq(path, decode_threads)?;
    let mut buffer = vec![0u8; chunk_size_for(decode_threads)];
    let mut lines = 0;
    let mut last = None;
    loop {
        let n = match reader.read(&mut buffer) {
            Ok(0) => break,
            Ok(n) => n,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(error) => {
                return Err(error).with_context(|| format!("Failed reading {}", path.display()));
            }
        };
        lines += bytecount::count(&buffer[..n], b'\n');
        last = Some(buffer[n - 1]);
    }
    // A final line without a line feed still counts.
    if last.is_some_and(|b| b != b'\n') {
        lines += 1;
    }
    Ok(lines.div_ceil(4))
}

pub struct FastqSummaryMin {
    pub path: PathBuf,
    pub file_name: String,
    pub read_count: usize,
}

impl FastqSummaryMin {
    pub fn new(path: &Path) -> Self {
        Self {
            read_count: 0,
            file_name: file_name(path),
            path: path.to_path_buf(),
        }
    }

    /// Count reads. The compression is detected from the file content.
    pub fn summarize(&mut self, _file_fmt: &SeqReadFmt) {
        self.read_count = count_reads(&self.path, 1).expect("Failed counting reads");
    }
}

pub struct FastqSummary {
    /// Input file path
    pub path: PathBuf,
    /// File name with extension
    pub file_name: String,
    pub reads: ReadRecord,
    pub read_summary: ReadSummary,
    pub qscores: ReadQScore,
    /// All statistics for the file.
    pub report: SampleReadReport,
}

impl FastqSummary {
    pub fn new(path: &Path) -> Self {
        Self {
            path: path.to_path_buf(),
            file_name: file_name(path),
            reads: ReadRecord::new(),
            read_summary: ReadSummary::new(),
            qscores: ReadQScore::new(),
            report: SampleReadReport::default(),
        }
    }

    /// Wrap a report, filling the summary fields kept for older callers.
    pub fn from_report(path: &Path, report: SampleReadReport) -> Self {
        let mut summary = Self::new(path);
        summary.set_report(report);
        summary
    }

    /// Compute default statistics. The compression is detected from the file content.
    pub fn summarize(&mut self, _file_fmt: &SeqReadFmt) {
        let report = FastqAnalyzer::new(&self.path, PhredOffset::Auto, false, 1)
            .analyze()
            .expect("Failed summarizing reads");
        self.set_report(report);
    }

    fn set_report(&mut self, report: SampleReadReport) {
        let m = &report.metrics;
        self.reads = ReadRecord {
            stats: stream_stats(
                m.read_count,
                m.mean_read_length,
                m.min_read_length,
                m.max_read_length,
            ),
            len: m.base_count as usize,
            g_count: m.g_count as usize,
            c_count: m.c_count as usize,
            a_count: m.a_count as usize,
            t_count: m.t_count as usize,
            n_count: m.n_count as usize,
        };
        // Empty files keep zero content, as before.
        if m.base_count > 0 {
            self.read_summary.summarize(&self.reads);
        }
        self.qscores = ReadQScore {
            len: m.base_count as usize,
            low_q: m.low_qscore_bases as usize,
            q_sum: (m.mean_qscore * m.base_count as f64).round() as usize,
            stats: stream_stats(m.base_count, m.mean_qscore, m.min_qscore, m.max_qscore),
        };
        self.report = report;
    }

    pub fn metrics(&self) -> &ReadMetrics {
        &self.report.metrics
    }
}

fn stream_stats(count: u64, mean: f64, min: u64, max: u64) -> StreamStats {
    let has_values = count > 0;
    StreamStats {
        count: count as usize,
        mean,
        min: has_values.then_some(min as usize),
        max: has_values.then_some(max as usize),
    }
}

fn file_name(path: &Path) -> String {
    path.file_name()
        .expect("No file name")
        .to_str()
        .expect("File name not valid UTF-8")
        .to_string()
}

#[cfg(test)]
mod test {
    use std::path::Path;

    use super::*;

    const READ_1: &str = "tests/files/raw/read_1.fastq";

    #[test]
    fn test_read_count() {
        let path = Path::new(READ_1);
        let file_fmt = SeqReadFmt::Fastq;
        let mut summary = FastqSummaryMin::new(path);
        summary.summarize(&file_fmt);
        assert_eq!(summary.path, PathBuf::from(path));
        assert_eq!(summary.read_count, 2);
    }

    #[test]
    fn test_summary_fields() {
        let mut summary = FastqSummary::new(Path::new(READ_1));
        summary.summarize(&SeqReadFmt::Auto);
        assert_eq!(summary.reads.stats.count, 2);
        assert_eq!(summary.reads.len, 72);
        assert_eq!(summary.reads.g_count, 21);
        assert_eq!(summary.metrics().median_read_length, 36.0);
        assert!(summary.report.positions.is_empty());
    }

    #[test]
    fn test_complete_report() {
        let report = FastqAnalyzer::new(Path::new(READ_1), PhredOffset::Auto, true, 1)
            .analyze()
            .unwrap();
        assert_eq!(report.positions.len(), 36);
        assert_eq!(report.metrics.phred_offset, 33);
        assert!(report.duplication.is_some());
        assert_eq!(report.adapters.len(), 5);
        let json = report::ReadReport {
            samples: vec![report.clone()],
        }
        .to_json()
        .unwrap();
        // serde_json may round the last digit of a float, so compare counts.
        let parsed = report::ReadReport::from_json(&json).unwrap();
        assert_eq!(parsed.samples[0].positions.len(), report.positions.len());
        assert_eq!(parsed.samples[0].gc_distribution.len(), 101);
        assert_eq!(parsed.samples[0].metrics.q30_bases, 69);
    }

    #[test]
    fn test_empty_file() {
        let report = FastqAnalyzer::new(
            Path::new("tests/files/raw/read2_1.fq"),
            PhredOffset::Auto,
            true,
            1,
        )
        .analyze()
        .unwrap();
        assert_eq!(report.metrics.read_count, 0);
        assert_eq!(
            count_reads(Path::new("tests/files/raw/read2_1.fq"), 1).unwrap(),
            0
        );
    }

    const QC: &str = "tests/files/raw_qc/qc_adapter.fastq";

    fn analyze(path: &str, phred: PhredOffset) -> SampleReadReport {
        FastqAnalyzer::new(Path::new(path), phred, true, 2)
            .analyze()
            .unwrap()
    }

    #[test]
    fn test_qc_statistics() {
        let report = analyze(QC, PhredOffset::Auto);
        let m = &report.metrics;
        assert_eq!(m.read_count, 6);
        assert_eq!(m.base_count, 127);
        assert_eq!(m.median_read_length, 24.0);
        assert_eq!(m.n50_read_length, 24);
        assert_eq!(m.min_read_length, 5);
        assert_eq!(m.max_read_length, 30);
        assert_eq!(m.adapter_reads_pct, Some(50.0));
        assert!((m.poly_g_reads_pct.unwrap() - 100.0 / 6.0).abs() < 1e-9);
        let dup = report.duplication.as_ref().unwrap();
        assert_eq!(dup.reads, 6);
        assert_eq!(dup.distinct_sequences, 4.0);
        assert!(dup.exact);
        assert_eq!(report.overrepresented.len(), 1);
        assert_eq!(report.overrepresented[0].count, 3);
        assert_eq!(report.overrepresented[0].source, "Illumina Universal");
        // Illumina Universal starts at position 7 in three reads.
        let universal = &report.adapters[0].cumulative_pct;
        assert_eq!(universal[5], 0.0);
        assert_eq!(universal[6], 50.0);
        assert_eq!(universal[29], 50.0);
        // Only the short read has Q2 at its last position.
        assert_eq!(report.positions[4].min_qscore, 2);
        assert_eq!(m.low_quality_reads, 1);
        assert_eq!(report.read_quality_distribution[16], 1);
    }

    #[test]
    fn test_phred64_matches_phred33() {
        let phred33 = analyze(QC, PhredOffset::Auto);
        let phred64 = analyze("tests/files/raw_qc/qc_phred64.fastq", PhredOffset::Auto);
        assert_eq!(phred64.metrics.phred_offset, 64);
        assert_eq!(phred64.positions, phred33.positions);
        assert_eq!(phred64.metrics.q30_bases, phred33.metrics.q30_bases);
        // Forcing the wrong offset reports an error.
        let wrong = FastqAnalyzer::new(Path::new(QC), PhredOffset::Phred64, false, 1).analyze();
        assert!(wrong.is_err());
    }

    #[test]
    fn test_compressed_inputs_match_plain() {
        let plain = analyze(QC, PhredOffset::Auto);
        for path in [
            "tests/files/raw_qc/qc_gzip.fastq.gz",
            "tests/files/raw_qc/qc_bgzf.fastq.gz",
        ] {
            let report = analyze(path, PhredOffset::Auto);
            assert_eq!(report.metrics, plain.metrics);
            assert_eq!(report.positions, plain.positions);
        }
        assert_eq!(
            count_reads(Path::new("tests/files/raw_qc/qc_bgzf.fastq.gz"), 2).unwrap(),
            6
        );
    }

    #[test]
    fn test_small_chunks_match_one_chunk() {
        let whole = analyze(QC, PhredOffset::Phred33);
        // The offset is fixed because detection only sees the first block.
        let mut analyzer = FastqAnalyzer::new(Path::new(QC), PhredOffset::Phred33, true, 1);
        // Smaller than one record, so every block is cut and carried over.
        analyzer.chunk_size = 7;
        let split = analyzer.analyze().unwrap();
        assert_eq!(split, whole);
    }

    #[test]
    fn test_malformed_record() {
        let error = FastqAnalyzer::new(
            Path::new("tests/files/raw_invalid/malformed.fastq"),
            PhredOffset::Auto,
            false,
            1,
        )
        .analyze()
        .unwrap_err();
        let message = format!("{:#}", error);
        assert!(message.contains("malformed.fastq"));
        assert!(message.contains("Read 2"));
    }
}
