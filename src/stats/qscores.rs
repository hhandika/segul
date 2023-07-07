//! Q-Score per read
use super::common::StreamStats;

const LOW_QSCORE: u8 = 20;

#[derive(Debug, Clone, PartialEq)]
pub struct ReadQScore {
    /// Q-Score length
    pub len: usize,
    pub low_q: u8,
    pub q_sum: usize,
    pub stats: StreamStats,
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
            q_sum: 0,
            stats: StreamStats::new(),
        }
    }

    /// Summarize Q-Score per read
    pub fn summarize(&mut self, qscore: &[u8]) {
        self.len += qscore.len();
        qscore.iter().for_each(|r| self.update(r));
    }

    /// Update Q-Score per read
    /// Q-Score less than 20 is considered low quality
    pub fn update(&mut self, qscore: &u8) {
        self.q_sum += *qscore as usize;
        self.low_q += if *qscore < LOW_QSCORE { 1 } else { 0 };
        self.stats.update(self.q_sum, &usize::from(*qscore));
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn test_qscore_counts() {
        let mut qread = ReadQScore::new();
        qread.summarize(&[31, 34, 35, 35, 15, 30]);
        assert_eq!(qread.len, 6);
        assert_eq!(qread.low_q, 1);
        assert_eq!(qread.stats.min.unwrap_or(0), 15);
        assert_eq!(qread.stats.max.unwrap_or(0), 35);
        assert_approx_eq!(qread.stats.mean, 29.2, 1f64);
    }

    #[test]
    fn test_multi_read_qscore_count() {
        let read_qscores = vec![
            vec![31, 34, 35, 35, 15, 25],
            vec![31, 34, 35, 35, 15, 25],
            vec![31, 34, 35, 35, 15, 25],
        ];
        let mut qread = ReadQScore::new();
        read_qscores.iter().for_each(|x| qread.summarize(x));
        assert_eq!(qread.len, 18);
        assert_eq!(qread.low_q, 3);
        assert_eq!(qread.stats.min.unwrap_or(0), 15);
        assert_eq!(qread.stats.max.unwrap_or(0), 35);
        assert_approx_eq!(qread.stats.mean, 29.2, 1f64);
    }
}
