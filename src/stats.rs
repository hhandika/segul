//! A module for sequence statistics.

use std::collections::HashMap;
use std::path::Path;

use indexmap::IndexMap;

use crate::alignment::Alignment;
use crate::common::SeqFormat;
use crate::utils;

pub struct SiteStats {
    site_matrix: HashMap<usize, Vec<u8>>,
}

impl SiteStats {
    pub fn new() -> Self {
        Self {
            site_matrix: HashMap::new(),
        }
    }

    pub fn get_stats(&mut self, path: &Path, input_format: &SeqFormat) {
        let mut aln = Alignment::new();
        aln.get_aln_any(path, input_format);
        self.index_sites(&aln.alignment);
        aln.alignment.clear();
        let mut dna = DnaStats::new();
        dna.count_chars(&self.site_matrix);
        let (conserved, var_sites, parsim) = self.get_site_stats();
        println!("Sites: {}", utils::fmt_num(&(conserved + var_sites)));
        println!("Conserved sites: {}", utils::fmt_num(&conserved));
        println!("Variable sites: {}", utils::fmt_num(&var_sites));
        println!("Parsimony informative sites: {}", parsim);
        let all_chars: usize = aln.alignment.values().map(|seq| seq.len()).sum();
        println!("\nAll chars: {}", utils::fmt_num(&all_chars));
        println!("A: {}", utils::fmt_num(&dna.a_count));
        println!("C: {}", utils::fmt_num(&dna.c_count));
        println!("G: {}", utils::fmt_num(&dna.g_count));
        println!("T: {}", utils::fmt_num(&dna.t_count));
        println!("N: {}", utils::fmt_num(&dna.n_count));
        println!("?: {}", utils::fmt_num(&dna.missings));
        println!("-: {}", utils::fmt_num(&dna.gaps));
    }

    fn index_sites(&mut self, matrix: &IndexMap<String, String>) {
        matrix.values().for_each(|seq| {
            seq.bytes()
                .enumerate()
                .for_each(|(idx, dna)| match self.site_matrix.get_mut(&idx) {
                    Some(value) => match dna {
                        b'a' | b'g' | b't' | b'c' | b'A' | b'G' | b'T' | b'C' => value.push(dna),
                        _ => (), // ignore ambigous characters
                    },
                    None => match dna {
                        b'a' | b'g' | b't' | b'c' | b'A' | b'G' | b'T' | b'C' => {
                            self.site_matrix.insert(idx, vec![dna]);
                        }
                        _ => (),
                    },
                })
        });
    }

    fn get_site_stats(&mut self) -> (usize, usize, usize) {
        let mut parsim: usize = 0;
        let mut var_sites = 0;
        let mut conserved = 0;
        self.site_matrix.values().for_each(|site| {
            let (cons, var, n_patterns) = self.get_patterns(site);
            if n_patterns >= 2 {
                parsim += 1;
            }
            var_sites += var;
            conserved += cons;
        });

        (conserved, var_sites, parsim)
    }

    fn get_patterns(&self, site: &[u8]) -> (usize, usize, usize) {
        let mut uniques: Vec<u8> = site.to_vec();
        uniques.sort_unstable();
        uniques.dedup();

        // We consider variable sites
        // when the characters not all the same
        let mut var_sites = 0;
        let mut conserved = 0;
        let mut n_patterns = 0;
        if uniques.len() != 1 {
            var_sites += 1;
            uniques.iter().for_each(|ch| {
                let patterns = site.iter().filter(|&site| site == ch).count();
                if patterns >= 2 {
                    n_patterns += 1;
                }
            });
        } else {
            conserved += 1;
        }

        (conserved, var_sites, n_patterns)
    }
}

struct DnaStats {
    a_count: usize,
    c_count: usize,
    g_count: usize,
    t_count: usize,
    n_count: usize,
    missings: usize,
    gaps: usize,
    others: usize,
}

impl DnaStats {
    fn new() -> Self {
        Self {
            a_count: 0,
            c_count: 0,
            g_count: 0,
            t_count: 0,
            n_count: 0,
            missings: 0,
            gaps: 0,
            others: 0,
        }
    }

    fn count_chars(&mut self, site: &HashMap<usize, Vec<u8>>) {
        site.values().for_each(|seqs| {
            seqs.iter().for_each(|ch| match ch {
                b'a' | b'A' => self.a_count += 1,
                b'c' | b'C' => self.c_count += 1,
                b'g' | b'G' => self.g_count += 1,
                b't' | b'T' => self.t_count += 1,
                b'n' | b'N' => self.n_count += 1,
                b'?' | b'.' => self.missings += 1,
                b'-' => self.gaps += 1,
                _ => self.others += 1,
            })
        })
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
        let site = b"AATT";
        let site_2 = b"AATTGG";
        let (_, _, pattern) = SiteStats::new().get_patterns(site);
        let (_, _, pattern_2) = SiteStats::new().get_patterns(site_2);
        assert_eq!(2, pattern);
        assert_eq!(3, pattern_2);
    }

    #[test]
    fn count_parsimony_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT", "ATTA", "ATGC", "ATGA"];
        let mat = get_matrix(&id, &seq);
        let mut dna = SiteStats::new();
        dna.index_sites(&mat);
        let (_, _, parsim) = dna.get_site_stats();
        assert_eq!(1, parsim);
    }

    #[test]
    fn count_variable_sites_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT", "ATTA", "ATGC", "ATGA"];
        let mat = get_matrix(&id, &seq);
        let mut dna = SiteStats::new();
        dna.index_sites(&mat);
        let (_, _, parsim) = dna.get_site_stats();
        assert_eq!(1, parsim);
    }

    #[test]
    fn count_parsimony_gap_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT---", "ATTA---", "ATGC---", "ATGA---"];
        let mat = get_matrix(&id, &seq);
        let mut dna = SiteStats::new();
        dna.index_sites(&mat);
        let (_, var_sites, parsim) = dna.get_site_stats();
        assert_eq!(1, parsim);
        assert_eq!(3, var_sites);
    }
}
