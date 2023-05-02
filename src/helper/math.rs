//! Heru Handika
//! Module for statistics

use std::cmp::Reverse;

#[inline]
fn sort_vector_asc(vec: &[u32]) -> Vec<u32> {
    let mut sorted_vec = vec.to_vec();
    sorted_vec.sort_unstable();

    sorted_vec
}

pub fn median(vec: &[u32]) -> f64 {
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

#[inline(always)]
fn sum_of_square(vec: &[f64]) -> f64 {
    let d: f64 = vec.iter().map(|val| val.powf(2.0)).sum();
    d
}

#[inline(always)]
fn dev_mean(vec: &[u32], mean: &f64) -> Vec<f64> {
    vec.iter().map(|&val| val as f64 - *mean).collect()
}

fn variance(vec: &[u32], mean: &f64) -> f64 {
    let d_mean = dev_mean(vec, mean);
    let n = vec.len() as f64 - 1.0;
    sum_of_square(&d_mean) / n
}

pub fn stdev(vec: &[u32], mean: &f64) -> f64 {
    let var = variance(vec, mean);
    var.sqrt()
}

fn sort_vec_desc(vec: &[u32]) -> Vec<u32> {
    let mut sorted_vec = vec.to_vec();
    sorted_vec.sort_by_key(|v| Reverse(*v));

    sorted_vec
}

fn cumsum(vec: &[u32]) -> Vec<u32> {
    let mut csum = Vec::new();
    let mut sum = 0;
    vec.iter().for_each(|v| {
        sum += v;
        csum.push(sum);
    });
    csum
}

pub struct NStats {
    sorted_contigs: Vec<u32>,
    csum_contigs: Vec<u32>,
    sum_contigs: u32,
    pub n50: u32,
    pub n75: u32,
    pub n90: u32,
}

impl NStats {
    pub fn new(contigs: &[u32]) -> Self {
        let mut nstats = Self {
            sorted_contigs: sort_vec_desc(contigs),
            csum_contigs: Vec::new(),
            sum_contigs: contigs.iter().sum::<u32>(),
            n50: 0,
            n75: 0,
            n90: 0,
        };

        nstats.csum_contigs = cumsum(&nstats.sorted_contigs);

        nstats
    }

    pub fn get_n50(&mut self) {
        let n50_len = self.n_len(0.5);
        let idx = self.get_n_idx(n50_len);
        self.n50 = self.sorted_contigs[idx];
    }

    pub fn get_n75(&mut self) {
        let n75_len = self.n_len(0.75);
        let idx = self.get_n_idx(n75_len);
        self.n75 = self.sorted_contigs[idx];
    }

    pub fn get_n90(&mut self) {
        let n90_len = self.n_len(0.9);
        let idx = self.get_n_idx(n90_len);
        self.n90 = self.sorted_contigs[idx];
    }

    fn get_n_idx(&mut self, n: u32) -> usize {
        self.csum_contigs.iter().position(|i| *i >= n).unwrap()
    }

    fn n_len(&mut self, i: f64) -> u32 {
        let n = self.sum_contigs as f64 * i;

        n as u32
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use assert_approx_eq::assert_approx_eq;

    #[test]
    fn median_test() {
        let odd: Vec<u32> = vec![1, 4, 3, 5, 6];
        let even: Vec<u32> = vec![1, 4, 3, 5, 6, 6, 8, 10];
        assert_eq!(4.0, median(&odd));
        assert_eq!(5.5, median(&even));
    }
    #[test]
    fn var_test() {
        let data: Vec<u32> = vec![1, 4, 3, 5, 6, 6, 8, 10];
        let mean = 5.375;
        let exp = 7.982143;
        let res = variance(&data, &mean);
        assert_approx_eq!(exp, res, 6f64);
    }

    #[test]
    fn stdev_test() {
        let data: Vec<u32> = vec![1, 4, 3, 5, 6, 6, 8, 10];
        let mean = 5.375;

        let exp = 2.825269;
        let res = stdev(&data, &mean);
        assert_approx_eq!(exp, res, 6f64);
    }

    #[test]
    fn csum_test() {
        let a = vec![1, 2, 3];
        let res = vec![1, 3, 6];
        assert_eq!(res, cumsum(&a));
    }

    #[test]
    fn sorted_vec_desc_test() {
        let a = vec![1, 2, 3];
        let res = vec![3, 2, 1];

        assert_eq!(res, sort_vec_desc(&a));
    }

    #[test]
    fn n50_stats_test() {
        let contigs = vec![2, 3, 4, 5, 6, 7, 8, 9, 10];
        let mut seq = NStats::new(&contigs);
        seq.get_n50();
        seq.get_n90();
        seq.get_n75();

        assert_eq!(8, seq.n50);
        assert_eq!(6, seq.n75);
        assert_eq!(4, seq.n90);
    }
}
