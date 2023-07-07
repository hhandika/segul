//! Data types for all FASTQ records
use std::collections::BTreeMap;

use super::common::StreamStats;

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
        self.len = qread.iter().map(|x| x.len).sum();
        self.low_q = qread.iter().map(|x| x.low_q).sum();
        self.sum = qread.iter().map(|x| x.sum).sum();
        self.mean = self.sum as f64 / self.len as f64;
        self.min = qread.iter().map(|x| x.min).min().unwrap_or(0);
        self.max = qread.iter().map(|x| x.max).max().unwrap_or(0);
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReadSummary {
    /// GC count in read
    pub gc_count: usize,
    /// GC content in read
    pub gc_content: f64,
    /// AT count in read
    pub at_count: usize,
    /// AT content in read
    pub at_content: f64,
    /// N content in read
    pub n_content: f64,
}

impl Default for ReadSummary {
    fn default() -> Self {
        Self::new()
    }
}

impl ReadSummary {
    pub fn new() -> Self {
        Self {
            gc_count: 0,
            gc_content: 0.0,
            at_count: 0,
            at_content: 0.0,
            n_content: 0.0,
        }
    }

    pub fn summarize(&mut self, record: &ReadRecord) {
        self.gc_count = record.g_count + record.c_count;
        self.gc_content = self.gc_count as f64 / record.len as f64;
        self.at_count = record.a_count + record.t_count;
        self.at_content = self.at_count as f64 / record.len as f64;
        self.n_content = record.n_count as f64 / record.len as f64;
    }
}

/// Statistics for a single file read records
#[derive(Debug, Clone, PartialEq)]
pub struct ReadRecord {
    /// read common stats
    pub stats: StreamStats,
    /// Read length
    pub len: usize,
    /// G count in read
    pub g_count: usize,
    /// C count in read
    pub c_count: usize,
    /// A count in read
    pub a_count: usize,
    /// T count in read
    pub t_count: usize,
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
            stats: StreamStats::new(),
            len: 0,
            g_count: 0,
            c_count: 0,
            a_count: 0,
            t_count: 0,
            n_count: 0,
        }
    }

    pub fn summarize(&mut self, read: &[u8]) {
        self.len += read.len();
        read.iter().for_each(|r| match r {
            b'G' | b'g' => self.g_count += 1,
            b'C' | b'c' => self.c_count += 1,
            b'A' | b'a' => self.a_count += 1,
            b'T' | b't' => self.t_count += 1,
            b'N' | b'n' => self.n_count += 1,
            _ => (),
        });
        self.stats.update(self.len, &read.len());
    }

    pub fn add(&mut self, base: &u8) {
        match base {
            b'G' | b'g' => self.g_count += 1,
            b'C' | b'c' => self.c_count += 1,
            b'A' | b'a' => self.a_count += 1,
            b'T' | b't' => self.t_count += 1,
            b'N' | b'n' => self.n_count += 1,
            _ => (),
        }
    }
}

pub struct FastqMappedRead {
    pub reads: BTreeMap<i32, ReadRecord>,
    pub qscores: BTreeMap<i32, QScoreStream>,
}

impl Default for FastqMappedRead {
    fn default() -> Self {
        Self::new()
    }
}

impl FastqMappedRead {
    pub fn new() -> Self {
        Self {
            reads: BTreeMap::new(),
            qscores: BTreeMap::new(),
        }
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
        self.len += qread.len();
        self.low_q = qread.iter().map(|x| usize::from(*x < 20)).sum();
        self.min = qread.iter().copied().min().unwrap_or(0);
        self.max = qread.iter().copied().max().unwrap_or(0);
        self.sum = qread.iter().map(|x| *x as usize).sum();
    }
}

pub struct QScoreStream {
    pub mean: f64,
    pub min: Option<u8>,
    pub max: Option<u8>,
    sum: usize,
    count: usize,
}

impl Default for QScoreStream {
    fn default() -> Self {
        Self::new()
    }
}

impl QScoreStream {
    pub fn new() -> Self {
        Self {
            mean: 0.0,
            min: None,
            max: None,
            sum: 0,
            count: 0,
        }
    }

    pub fn update(&mut self, score: &u8) {
        self.calculate_mean(score);
        self.calculate_min(score);
        self.calculate_max(score);
    }

    fn calculate_mean(&mut self, score: &u8) {
        self.sum += *score as usize;
        self.count += 1;
        self.mean = self.sum as f64 / self.count as f64;
    }

    fn calculate_min(&mut self, score: &u8) {
        if let Some(min) = self.min {
            if min > *score {
                self.min = Some(*score);
            }
        } else {
            self.min = Some(*score);
        }
    }

    fn calculate_max(&mut self, score: &u8) {
        if let Some(max) = self.max {
            if max < *score {
                self.max = Some(*score);
            }
        } else {
            self.max = Some(*score);
        }
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
        assert_eq!(read.g_count, 1);
        assert_eq!(read.c_count, 1);
        assert_eq!(read.a_count, 1);
        assert_eq!(read.t_count, 1);
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
    fn test_qscore_records() {
        let mut qscore = QScoreRecords::new();
        let qscore_records = vec![
            ReadQScore {
                len: 100,
                low_q: 0,
                sum: 3000,
                min: 30,
                max: 40,
            },
            ReadQScore {
                len: 100,
                low_q: 1,
                sum: 1000,
                min: 10,
                max: 33,
            },
        ];
        qscore.summarize(&qscore_records);
        assert_eq!(qscore.len, 200);
        assert_eq!(qscore.low_q, 1);
        assert_eq!(qscore.sum, 4000);
        assert_eq!(qscore.mean, 20.0);
        assert_eq!(qscore.min, 10);
        assert_eq!(qscore.max, 40);
    }

    #[test]
    fn test_streaming_qscore() {
        let qscore = vec![40, 40, 10, 30, 30];
        let mut score = QScoreStream::new();
        for s in qscore {
            score.update(&s);
        }
        assert_eq!(score.mean, 30.0);
        assert_eq!(score.min, Some(10));
        assert_eq!(score.max, Some(40));
    }
}
