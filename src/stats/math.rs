//! Heru Handika
//! Module for statistics
use std::cmp::Reverse;

pub fn mean(vec: &[usize], sum: usize) -> f64 {
    let n = vec.len() as f64;
    sum as f64 / n
}

pub fn median(vec: &[usize]) -> f64 {
    let sorted_vec = sort_vector_asc(&vec);
    let n = sorted_vec.len();
    let midpoint = n / 2;

    let med;
    if n % 2 == 0 {
        med = (sorted_vec[midpoint - 1] + sorted_vec[midpoint]) as f64 / 2.0
    } else {
        med = sorted_vec[midpoint] as f64
    }

    med
}

pub fn stdev(vec: &[usize], mean: &f64) -> f64 {
    let var = variance(vec, mean);
    var.sqrt()
}

pub struct NStats {
    pub n50: usize,
    pub n75: usize,
    pub n90: usize,
    sum: usize,
}

impl NStats {
    pub fn new(sum: usize) -> Self {
        Self {
            n50: 0,
            n75: 0,
            n90: 0,
            sum,
        }
    }

    pub fn calculate(&mut self, contigs: &[usize]) {
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

#[inline]
fn sort_vector_asc(vec: &[usize]) -> Vec<usize> {
    let mut sorted_vec = vec.to_vec();
    sorted_vec.sort_unstable();

    sorted_vec
}

#[inline(always)]
fn sum_of_square(vec: &[f64]) -> f64 {
    let d: f64 = vec.iter().map(|val| val.powf(2.0)).sum();
    d
}

#[inline(always)]
fn dev_mean(vec: &[usize], mean: &f64) -> Vec<f64> {
    vec.iter().map(|&val| val as f64 - *mean).collect()
}

fn variance(vec: &[usize], mean: &f64) -> f64 {
    let d_mean = dev_mean(vec, mean);
    let n = vec.len() as f64 - 1.0;
    sum_of_square(&d_mean) / n
}

#[cfg(test)]
mod test {
    use super::*;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn median_test() {
        let odd: Vec<usize> = vec![1, 4, 3, 5, 6];
        let even: Vec<usize> = vec![1, 4, 3, 5, 6, 6, 8, 10];
        assert_eq!(4.0, median(&odd));
        assert_eq!(5.5, median(&even));
    }
    #[test]
    fn var_test() {
        let data: Vec<usize> = vec![1, 4, 3, 5, 6, 6, 8, 10];
        let mean = 5.375;
        let exp = 7.982143;
        let res = variance(&data, &mean);
        assert_approx_eq!(exp, res, 6f64);
    }

    #[test]
    fn stdev_test() {
        let data: Vec<usize> = vec![1, 4, 3, 5, 6, 6, 8, 10];
        let mean = 5.375;

        let exp = 2.825269;
        let res = stdev(&data, &mean);
        assert_approx_eq!(exp, res, 6f64);
    }

    #[test]
    fn csum_test() {
        let a = vec![1, 2, 3];
        let res = vec![1, 3, 6];
        let nstat = NStats::new(6);
        assert_eq!(res, nstat.cumsum(&a));
    }

    #[test]
    fn sorted_vec_desc_test() {
        let a = vec![1, 2, 3];
        let res = vec![3, 2, 1];
        let nstat = NStats::new(6);
        assert_eq!(res, nstat.sort_vec_desc(&a));
    }

    #[test]
    fn n50_stats_test() {
        let contigs = vec![2, 3, 4, 5, 6, 7, 8, 9, 10];
        let sum = contigs.iter().sum::<usize>();
        let mut seq = NStats::new(sum);
        seq.calculate(&contigs);
        assert_eq!(8, seq.n50);
        assert_eq!(6, seq.n75);
        assert_eq!(4, seq.n90);
    }
}
