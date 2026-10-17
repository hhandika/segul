//! Data types for all FASTQ records
use super::common::StreamStats;

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

/// Read summary statistics for a single file
impl ReadSummary {
    /// Create a new ReadSummary instance
    pub fn new() -> Self {
        Self {
            gc_count: 0,
            gc_content: 0.0,
            at_count: 0,
            at_content: 0.0,
            n_content: 0.0,
        }
    }
    /// Summarize a read record and update the stats
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
    /// Create a new ReadRecord instance
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

    /// Summarize a read record and update the stats
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

    /// Add a base to the read record
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
    fn test_multi_read_stats() {
        let reads = [b"ATGCN", b"ATGCC"];
        let mut read = ReadRecord::new();
        reads.iter().for_each(|&r| read.summarize(r));
        let mut summary = ReadSummary::new();
        summary.summarize(&read);
        assert_eq!(read.len, 10);
        assert_eq!(read.g_count, 2);
        assert_eq!(read.c_count, 3);
        assert_eq!(read.a_count, 2);
        assert_eq!(read.t_count, 2);
        assert_eq!(read.n_count, 1);
        assert_eq!(read.stats.count, 2);
        assert_eq!(read.stats.min.unwrap_or(0), 5);
        assert_eq!(read.stats.max.unwrap_or(0), 5);
        assert_eq!(read.stats.mean, 5.0);
        assert_eq!(summary.gc_count, 5);
        assert_eq!(summary.gc_content, 0.5);
        assert_eq!(summary.at_count, 4);
        assert_eq!(summary.at_content, 0.4);
        assert_eq!(summary.n_content, 0.1);
    }
}
