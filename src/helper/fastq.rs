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
        self.max = qread.iter().map(|x| x.max).max().unwrap_or(0);
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
        read.iter().for_each(|r| match r {
            b'G' | b'g' | b'C' | b'c' => self.gc_count += 1,
            b'A' | b'T' | b'a' | b't' => self.at_count += 1,
            b'N' | b'n' => self.n_count += 1,
            _ => (),
        });
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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_base_counts() {
        let mut read = ReadRecord::new();
        read.summarize(b"ATGC");
        assert_eq!(read.len, 4);
        assert_eq!(read.gc_count, 2);
        assert_eq!(read.at_count, 2);
        assert_eq!(read.n_count, 0);
    }

    #[test]
    fn test_qscore_counts() {
        let mut qread = ReadQScore::new();
        qread.summarize(&[31, 34, 35, 35, 15, 20]);
        assert_eq!(qread.len, 6);
        assert_eq!(qread.low_q, 1);
        assert_eq!(qread.sum, 170);
        assert_eq!(qread.min, 15);
        assert_eq!(qread.max, 35);
    }

    #[test]
    fn test_fastq_records() {
        let mut fq = FastqRecords::new(Path::new("test.fq"));
        let read_records = vec![
            ReadRecord {
                len: 100,
                gc_count: 50,
                at_count: 50,
                n_count: 0,
            },
            ReadRecord {
                len: 100,
                gc_count: 50,
                at_count: 50,
                n_count: 0,
            },
        ];
        fq.summarize(&read_records);
        assert_eq!(fq.num_reads, 2);
        assert_eq!(fq.num_bases, 200);
        assert_eq!(fq.min_read_len, 100);
        assert_eq!(fq.mean_read_len, 100);
        assert_eq!(fq.max_read_len, 100);
    }

    #[test]
    fn test_qscore_records() {
        let mut qscore = QScoreRecords::new();
        let qscore_records = vec![
            ReadQScore {
                len: 100,
                low_q: 0,
                sum: 100,
                min: 30,
                max: 40,
            },
            ReadQScore {
                len: 100,
                low_q: 1,
                sum: 100,
                min: 10,
                max: 33,
            },
        ];
        qscore.summarize(&qscore_records);
        assert_eq!(qscore.len, 2);
        assert_eq!(qscore.low_q, 1);
        assert_eq!(qscore.sum, 200);
        assert_eq!(qscore.mean, 100.0);
        assert_eq!(qscore.min, 10);
        assert_eq!(qscore.max, 40);
    }
}
