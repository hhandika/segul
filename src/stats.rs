//! A module for sequence statistics.

use std::collections::HashMap;
use std::path::Path;

use indexmap::IndexMap;

use crate::alignment::Alignment;
use crate::common::SeqFormat;

#[allow(dead_code)]
pub struct AlnStats {
    parsimony_inf: usize,
    site_matrix: HashMap<usize, String>,
}

#[allow(dead_code)]
impl AlnStats {
    pub fn new() -> Self {
        Self {
            parsimony_inf: 0,
            site_matrix: HashMap::new(),
        }
    }

    pub fn get_stats(&mut self, path: &Path, input_format: &SeqFormat, interleave: bool) {
        let mut aln = Alignment::new();
        aln.get_aln_any(path, input_format, interleave);
        self.index_sites(&aln.alignment);
        self.parsimony_inf = self.count_parsimony_informative();
        println!("Parsimony informative sites: {}", self.parsimony_inf);
    }

    fn index_sites(&mut self, matrix: &IndexMap<String, String>) {
        matrix.values().for_each(|seq| {
            seq.chars()
                .enumerate()
                .for_each(|(idx, dna)| match self.site_matrix.get_mut(&idx) {
                    Some(value) => match dna {
                        '-' | 'N' | '?' | '.' => (),
                        _ => value.push(dna),
                    },
                    None => match dna {
                        '-' | 'N' | '?' | '.' => (),
                        _ => {
                            self.site_matrix.insert(idx, dna.to_string());
                        }
                    },
                })
        });
    }

    fn count_parsimony_informative(&mut self) -> usize {
        let mut parsim: usize = 0;
        self.site_matrix.values().for_each(|site| {
            let n_patterns = self.get_pattern(&site);
            if n_patterns >= 2 {
                parsim += 1
            }
        });

        parsim
    }

    fn get_pattern(&self, site: &str) -> usize {
        let mut uniques: Vec<char> = site.chars().collect();
        uniques.sort_unstable();
        uniques.dedup();
        let mut n_patterns = 0;
        uniques.iter().for_each(|c| {
            let patterns = site.matches(&c.to_string()).count();
            if patterns >= 2 {
                n_patterns += 1;
            }
        });

        n_patterns
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn get_matrix(id: &[&str], seq: &[&str]) -> IndexMap<String, String> {
        let mut matrix = IndexMap::new();
        id.iter().zip(seq.iter()).for_each(|(i, s)| {
            matrix.insert(i.to_string(), s.to_string());
        });

        matrix
    }

    #[test]
    fn pattern_count_test() {
        let site = "AATT";
        let site_2 = "AATTGG";
        let pattern = AlnStats::new().get_pattern(&site);
        assert_eq!(2, pattern);
        assert_eq!(3, AlnStats::new().get_pattern(site_2));
    }

    #[test]
    fn pattern_count_all_test() {
        let site = "AAAA";
        let pattern = AlnStats::new().get_pattern(&site);
        assert_eq!(1, pattern);
    }

    #[test]
    fn count_parsimony_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT", "ATTA", "ATGC", "ATGA"];
        let mat = get_matrix(&id, &seq);
        let mut dna = AlnStats::new();
        dna.index_sites(&mat);
        let parsim = dna.count_parsimony_informative();
        assert_eq!(1, parsim);
    }

    #[test]
    fn count_parsimony_gap_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT---", "ATTA---", "ATGC---", "ATGA---"];
        let mat = get_matrix(&id, &seq);
        let mut dna = AlnStats::new();
        dna.index_sites(&mat);
        let parsim = dna.count_parsimony_informative();
        assert_eq!(1, parsim);
    }
}
