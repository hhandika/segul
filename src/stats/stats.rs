//! Statistics for a vector of values and genomic sequences
use std::cmp::Reverse;

/// Common statistics for a vector of values
/// Mean, median, min, max, and stdev
///
/// # Example
/// ```
/// use assert_approx_eq::assert_approx_eq;
///
/// use segul::stats::stats::CommonStats;
///
/// let values = vec![1, 2, 3, 4, 5];
/// let mut stats = CommonStats::new();
/// stats.calculate(&values);
///
/// assert_eq!(stats.sum, 15);
/// assert_eq!(stats.mean, 3.0);
/// assert_eq!(stats.median, 3.0);
/// assert_eq!(stats.min, 1);
/// assert_eq!(stats.max, 5);
/// assert_approx_eq!(stats.stdev, 1.58, 2f64);
/// ```
pub struct CommonStats {
    pub sum: usize,
    pub mean: f64,
    pub median: f64,
    pub min: usize,
    pub max: usize,
    pub stdev: f64,
}

impl Default for CommonStats {
    fn default() -> Self {
        Self::new()
    }
}

impl CommonStats {
    /// Create a new instance of CommonStats
    pub fn new() -> Self {
        Self {
            sum: 0,
            mean: 0.0,
            median: 0.0,
            min: 0,
            max: 0,
            stdev: 0.0,
        }
    }

    /// Calculate the common statistics
    pub fn calculate(&mut self, values: &[usize]) {
        self.sum(values);
        self.mean(values);
        self.median(values);
        self.min(values);
        self.max(values);
        self.stdev(values);
    }

    fn sum(&mut self, values: &[usize]) {
        self.sum = values.iter().sum();
    }

    fn min(&mut self, values: &[usize]) {
        self.min = *values.iter().min().unwrap_or(&0);
    }

    fn max(&mut self, values: &[usize]) {
        self.max = *values.iter().max().unwrap_or(&0);
    }

    fn mean(&mut self, vec: &[usize]) {
        let n = vec.len() as f64;
        self.mean = self.sum as f64 / n;
    }

    fn median(&mut self, values: &[usize]) {
        let sorted_vec = self.sort_vector_asc(values);
        let n = sorted_vec.len();
        let midpoint = n / 2;

        if n % 2 == 0 {
            self.median = (sorted_vec[midpoint - 1] + sorted_vec[midpoint]) as f64 / 2.0;
        } else {
            self.median = sorted_vec[midpoint] as f64;
        }
    }

    fn stdev(&mut self, values: &[usize]) {
        let var = self.variance(values);
        self.stdev = var.sqrt();
    }

    fn variance(&self, values: &[usize]) -> f64 {
        let d_mean = self.dev_mean(values);
        let n = values.len() as f64 - 1.0;
        self.sum_of_square(&d_mean) / n
    }

    #[inline]
    fn sort_vector_asc(&self, values: &[usize]) -> Vec<usize> {
        let mut sorted_vec = values.to_vec();
        sorted_vec.sort_unstable();

        sorted_vec
    }

    #[inline(always)]
    fn sum_of_square(&self, vec: &[f64]) -> f64 {
        let d: f64 = vec.iter().map(|val| val.powf(2.0)).sum();
        d
    }

    #[inline(always)]
    fn dev_mean(&self, vec: &[usize]) -> Vec<f64> {
        vec.iter().map(|&val| val as f64 - self.mean).collect()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct StreamStats {
    pub count: usize,
    pub mean: f64,
    pub min: Option<usize>,
    pub max: Option<usize>,
}

impl Default for StreamStats {
    fn default() -> Self {
        Self::new()
    }
}

impl StreamStats {
    pub fn new() -> Self {
        Self {
            count: 0,
            mean: 0.0,
            min: None,
            max: None,
        }
    }

    pub fn update(&mut self, sum: usize, values: &usize) {
        self.calculate_mean(sum);
        self.calculate_min(values);
        self.calculate_max(values);
    }

    fn calculate_mean(&mut self, sum: usize) {
        self.count += 1;
        self.mean = sum as f64 / self.count as f64;
    }

    fn calculate_min(&mut self, values: &usize) {
        if let Some(min) = self.min {
            if min > *values {
                self.min = Some(*values);
            }
        } else {
            self.min = Some(*values);
        }
    }

    fn calculate_max(&mut self, values: &usize) {
        if let Some(max) = self.max {
            if max < *values {
                self.max = Some(*values);
            }
        } else {
            self.max = Some(*values);
        }
    }
}
/// N50, N75, and N90 statistics for a vector of values
/// N50 is the length of the shortest contig at 50% of the total length
/// N75 is the length of the shortest contig at 75% of the total length
/// N90 is the length of the shortest contig at 90% of the total length
///
/// # Example
/// ```
/// use segul::stats::stats::NStats;
///
/// // Contig lengths stored in a vector
/// let contigs = vec![2, 3, 4, 5, 6, 7, 8, 9, 10];
/// let total_len = contigs.iter().sum::<usize>();
/// let mut nstats = NStats::new();
/// nstats.calculate(&contigs, total_len);
///
/// assert_eq!(8, nstats.n50);
/// assert_eq!(6, nstats.n75);
/// assert_eq!(4, nstats.n90);
/// ```
pub struct NStats {
    pub n50: usize,
    pub n75: usize,
    pub n90: usize,
    sum: usize,
}

impl NStats {
    pub fn new() -> Self {
        Self {
            n50: 0,
            n75: 0,
            n90: 0,
            sum: 0,
        }
    }

    pub fn calculate(&mut self, contigs: &[usize], total_len: usize) {
        self.sum = total_len;
        let sorted_contig = self.sort_vec_desc(contigs);
        let csum = self.cumsum(&sorted_contig);
        self.n50(&sorted_contig, &csum);
        self.n75(&sorted_contig, &csum);
        self.n90(&sorted_contig, &csum);
    }

    fn n50(&mut self, sorted_contigs: &[usize], csum: &[usize]) {
        let n50_len = self.n_len(0.5);
        let idx = self.n_idx(n50_len, csum);
        self.n50 = sorted_contigs[idx];
    }

    fn n75(&mut self, sorted_contigs: &[usize], csum: &[usize]) {
        let n75_len = self.n_len(0.75);
        let idx = self.n_idx(n75_len, csum);
        self.n75 = sorted_contigs[idx];
    }

    fn n90(&mut self, sorted_contigs: &[usize], csum: &[usize]) {
        let n90_len = self.n_len(0.9);
        let idx = self.n_idx(n90_len, csum);
        self.n90 = sorted_contigs[idx];
    }

    fn sort_vec_desc(&self, vec: &[usize]) -> Vec<usize> {
        let mut sorted_vec = vec.to_vec();
        sorted_vec.sort_by_key(|v| Reverse(*v));

        sorted_vec
    }

    fn cumsum(&self, vec: &[usize]) -> Vec<usize> {
        let mut csum = Vec::new();
        let mut sum = 0;
        vec.iter().for_each(|v| {
            sum += v;
            csum.push(sum);
        });
        csum
    }

    fn n_idx(&mut self, n: usize, csum: &[usize]) -> usize {
        csum.iter().position(|i| *i >= n).unwrap()
    }

    fn n_len(&mut self, i: f64) -> usize {
        let n = self.sum as f64 * i;

        n as usize
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn median_test() {
        let odd: Vec<usize> = vec![1, 4, 3, 5, 6];
        let even: Vec<usize> = vec![1, 4, 3, 5, 6, 6, 8, 10];
        let mut stat_odd = CommonStats::new();
        stat_odd.calculate(&odd);
        let mut stat_even = CommonStats::new();
        stat_even.calculate(&even);
        assert_eq!(4.0, stat_odd.median);
        assert_eq!(5.5, stat_even.median);
    }
    #[test]
    fn var_test() {
        let data: Vec<usize> = vec![1, 4, 3, 5, 6, 6, 8, 10];
        let exp = 7.982143;
        let mut stat = CommonStats::new();
        stat.calculate(&data);
        let res = stat.variance(&data);
        assert_approx_eq!(exp, res, 6f64);
    }

    #[test]
    fn stdev_test() {
        let data: Vec<usize> = vec![1, 4, 3, 5, 6, 6, 8, 10];
        let exp = 2.825269;
        let mut stat = CommonStats::new();
        stat.calculate(&data);
        assert_approx_eq!(exp, stat.stdev, 6f64);
    }

    #[test]
    fn csum_test() {
        let a = vec![1, 2, 3];
        let res = vec![1, 3, 6];
        let nstat = NStats::new();
        assert_eq!(res, nstat.cumsum(&a));
    }

    #[test]
    fn sorted_vec_desc_test() {
        let a = vec![1, 2, 3];
        let res = vec![3, 2, 1];
        let nstat = NStats::new();
        assert_eq!(res, nstat.sort_vec_desc(&a));
    }

    #[test]
    fn n50_stats_test() {
        let contigs = vec![2, 3, 4, 5, 6, 7, 8, 9, 10];
        let sum = contigs.iter().sum::<usize>();
        let mut seq = NStats::new();
        seq.calculate(&contigs, sum);
        assert_eq!(8, seq.n50);
        assert_eq!(6, seq.n75);
        assert_eq!(4, seq.n90);
    }
}
