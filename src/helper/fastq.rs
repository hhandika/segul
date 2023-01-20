//! Methods for processing FASTQ files.

use std::default::Default;
use std::path::{Path, PathBuf};

/// Data types for all FASTQ records
#[derive(Debug, Clone, PartialEq)]
pub struct FastqRecords {
    /// File path
    pub path: PathBuf,
    /// Number of reads
    pub num_reads: usize,
    /// Number of bases
    pub num_bases: usize,
    /// Average read length
    pub min_read_len: usize,
    /// Average read length
    pub mean_read_len: usize,
    /// Maximum read length
    pub max_read_len: usize,
}

impl FastqRecords {
    pub fn new(path: &Path) -> Self {
        Self {
            path: PathBuf::from(path),
            num_reads: 0,
            num_bases: 0,
            min_read_len: 0,
            mean_read_len: 0,
            max_read_len: 0,
        }
    }

    pub fn summarize(&mut self, read_records: &[ReadRecord]) {
        self.num_reads = read_records.len();
        self.num_bases = read_records.iter().map(|x| x.len).sum();
        self.min_read_len = read_records.iter().map(|x| x.len).min().unwrap();
        self.mean_read_len = self.num_bases / self.num_reads;
        self.max_read_len = read_records.iter().map(|x| x.len).max().unwrap();
    }
}

/// Data types for all Q-Score records
#[derive(Debug, Clone, PartialEq)]
pub struct QScoreRecords {
    /// Number of bases
    pub len: usize,
    /// Number of bases with quality score < 20
    pub low_q: usize,
    /// Sum of quality scores
    pub sum: usize,
    /// Mean quality score
    pub mean: f64,
    /// Minimum quality score
    pub min: u8,
    /// Maximum quality score
    pub max: u8,
}

impl Default for QScoreRecords {
    fn default() -> Self {
        Self::new()
    }
}

impl QScoreRecords {
    pub fn new() -> Self {
        Self {
            len: 0,
            low_q: 0,
            sum: 0,
            mean: 0.0,
            min: 0,
            max: 0,
        }
    }

    pub fn summarize(&mut self, qread: &[ReadQScore]) {
        self.len = qread.len();
        self.low_q = qread.iter().map(|x| x.low_q).sum();
        self.sum = qread.iter().map(|x| x.sum).sum();
        self.mean = self.sum as f64 / self.len as f64;
        self.min = qread.iter().map(|x| x.min).min().unwrap_or(0);
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReadRecord {
    /// Read length
    pub len: usize,
    /// GC count in read
    pub gc_count: usize,
    /// AT count in read
    pub at_count: usize,
    /// N count in read
    pub n_count: usize,
}

impl Default for ReadRecord {
    fn default() -> Self {
        Self::new()
    }
}

impl ReadRecord {
    pub fn new() -> Self {
        Self {
            len: 0,
            gc_count: 0,
            at_count: 0,
            n_count: 0,
        }
    }

    pub fn summarize(&mut self, read: &[u8]) {
        self.len = read.len();
        self.gc_count = read
            .iter()
            .map(|x| usize::from(*x == b'G' || *x == b'C'))
            .sum();
        self.at_count = read
            .iter()
            .map(|x| usize::from(*x == b'A' || *x == b'T'))
            .sum();
        self.n_count = read.iter().map(|x| usize::from(*x == b'N')).sum();
    }
}

/// Q-Score per read
#[derive(Debug, Clone, PartialEq)]
pub struct ReadQScore {
    /// Q-Score length
    pub len: usize,
    pub low_q: usize,
    pub sum: usize,
    pub min: u8,
    pub max: u8,
}

impl Default for ReadQScore {
    fn default() -> Self {
        Self::new()
    }
}

impl ReadQScore {
    pub fn new() -> Self {
        Self {
            len: 0,
            low_q: 0,
            sum: 0,
            min: 0,
            max: 0,
        }
    }

    pub fn summarize(&mut self, qread: &[u8]) {
        self.len = qread.len();
        self.low_q = qread.iter().map(|x| usize::from(*x < 20)).sum();
        self.min = qread.iter().copied().min().unwrap_or(0);
        self.max = qread.iter().copied().max().unwrap_or(0);
        self.sum = qread.iter().map(|x| *x as usize).sum();
    }
}
