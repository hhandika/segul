//! A module for sequence statistics.
use std::collections::HashMap;
use std::io::{self, BufWriter, Result, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;

use indexmap::IndexMap;
use rayon::prelude::*;

use crate::alignment::Alignment;
use crate::common::SeqFormat;
use crate::utils;

pub fn get_seq_stats(path: &Path, input_format: &SeqFormat) {
    let spin = utils::set_spinner();
    spin.set_message("Getting alignments...");
    let mut aln = Alignment::new();
    aln.get_aln_any(path, input_format);
    spin.set_message("Counting characters...");
    let mut dna = Dna::new();
    dna.count_chars(&aln.alignment);
    spin.set_message("Getting summary stats...");
    let mut sites = Sites::new();
    sites.get_stats(&aln.alignment);
    spin.finish_with_message("DONE!\n");
    display_stats(&sites, &dna).unwrap();
}

pub fn get_stats_dir(files: &[PathBuf], input_format: &SeqFormat) {
    let (send, rec) = channel();

    files.par_iter().for_each_with(send, |s, file| {
        s.send(get_stats(file, input_format)).unwrap();
    });

    let stats: Vec<(Dna, Sites)> = rec.iter().collect();

    display_summary(&stats);
}

fn get_stats(path: &Path, input_format: &SeqFormat) -> (Dna, Sites) {
    let mut aln = Alignment::new();
    aln.get_aln_any(path, input_format);
    let mut dna = Dna::new();
    dna.count_chars(&aln.alignment);
    let mut sites = Sites::new();
    sites.get_stats(&aln.alignment);

    (dna, sites)
}

fn display_summary(stats: &[(Dna, Sites)]) {
    stats.iter().for_each(|(dna, site)| {
        display_stats(site, dna).unwrap();
    })
}

fn display_stats(site: &Sites, dna: &Dna) -> Result<()> {
    let io = io::stdout();
    let mut writer = BufWriter::new(io);

    writeln!(writer, "\x1b[0;33mSites\x1b[0m")?;
    writeln!(writer, "Count\t\t: {}", utils::fmt_num(&site.counts))?;
    writeln!(writer, "Conserved\t: {}", utils::fmt_num(&site.conserved))?;
    writeln!(writer, "Variable\t: {}", utils::fmt_num(&site.variable))?;
    writeln!(writer, "Parsimony inf.\t: {}\n", site.pars_inf)?;

    writeln!(writer, "\x1b[0;33mCharacters\x1b[0m")?;
    writeln!(writer, "Total\t: {}", utils::fmt_num(&dna.nchars))?;
    writeln!(writer, "A\t: {}", utils::fmt_num(&dna.a_count))?;
    writeln!(writer, "C\t: {}", utils::fmt_num(&dna.c_count))?;
    writeln!(writer, "G\t: {}", utils::fmt_num(&dna.g_count))?;
    writeln!(writer, "T\t: {}", utils::fmt_num(&dna.t_count))?;
    writeln!(writer, "N\t: {}", utils::fmt_num(&dna.n_count))?;
    writeln!(writer, "?\t: {}", utils::fmt_num(&dna.missings))?;
    writeln!(writer, "-\t: {}", utils::fmt_num(&dna.gaps))?;

    Ok(())
}

pub struct Sites {
    conserved: usize,
    variable: usize,
    pars_inf: usize,
    counts: usize,
}

impl Sites {
    pub fn new() -> Self {
        Self {
            conserved: 0,
            variable: 0,
            pars_inf: 0,
            counts: 0,
        }
    }

    pub fn get_stats(&mut self, matrix: &IndexMap<String, String>) {
        let site_matrix = self.index_sites(matrix);
        self.get_site_stats(&site_matrix);
        self.count_sites();
    }

    fn index_sites(&mut self, matrix: &IndexMap<String, String>) -> HashMap<usize, Vec<u8>> {
        let mut site_matrix: HashMap<usize, Vec<u8>> = HashMap::new();
        matrix.values().for_each(|seq| {
            seq.bytes().enumerate().for_each(|(idx, dna)| {
                match site_matrix.get_mut(&idx) {
                    Some(value) => match dna {
                        b'a' | b'g' | b't' | b'c' | b'A' | b'G' | b'T' | b'C' => value.push(dna),
                        _ => (), // ignore ambigous characters
                    },
                    None => match dna {
                        b'a' | b'g' | b't' | b'c' | b'A' | b'G' | b'T' | b'C' => {
                            site_matrix.insert(idx, vec![dna]);
                        }
                        _ => (),
                    },
                }
            })
        });

        site_matrix
    }

    fn get_site_stats(&mut self, site_matrix: &HashMap<usize, Vec<u8>>) {
        site_matrix.values().for_each(|site| {
            let n_patterns = self.get_patterns(site);
            if n_patterns >= 2 {
                self.pars_inf += 1;
            }
        });
    }

    fn get_patterns(&mut self, site: &[u8]) -> usize {
        let mut uniques: Vec<u8> = site.to_vec();
        uniques.sort_unstable();
        uniques.dedup();

        // We consider variable sites
        // when the characters not all the same
        let mut n_patterns = 0;
        if uniques.len() >= 2 {
            self.variable += 1;
            uniques.iter().for_each(|ch| {
                let patterns = site.iter().filter(|&site| site == ch).count();
                if patterns >= 2 {
                    n_patterns += 1;
                }
            });
        } else {
            self.conserved += 1;
        }

        n_patterns
    }

    fn count_sites(&mut self) {
        self.counts = self.conserved + self.variable;
    }
}

struct Dna {
    a_count: usize,
    c_count: usize,
    g_count: usize,
    t_count: usize,
    n_count: usize,
    missings: usize,
    gaps: usize,
    others: usize,
    nchars: usize,
}

impl Dna {
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
            nchars: 0,
        }
    }

    fn count_chars(&mut self, matrix: &IndexMap<String, String>) {
        matrix.values().for_each(|seqs| {
            seqs.bytes().for_each(|ch| {
                match ch {
                    b'a' | b'A' => self.a_count += 1,
                    b'c' | b'C' => self.c_count += 1,
                    b'g' | b'G' => self.g_count += 1,
                    b't' | b'T' => self.t_count += 1,
                    b'n' | b'N' => self.n_count += 1,
                    b'?' | b'.' => self.missings += 1,
                    b'-' => self.gaps += 1,
                    _ => self.others += 1,
                }
                self.nchars += 1;
            })
        })
    }
}

#[cfg(test)]
mod test {
    // use indexmap::IndexMap;
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
        let pattern = Sites::new().get_patterns(site);
        let pattern_2 = Sites::new().get_patterns(site_2);
        assert_eq!(2, pattern);
        assert_eq!(3, pattern_2);
    }

    #[test]
    fn count_parsimony_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT", "ATTA", "ATGC", "ATGA"];
        let mat = get_matrix(&id, &seq);
        let mut site = Sites::new();
        let smat = site.index_sites(&mat);
        site.get_site_stats(&smat);
        assert_eq!(1, site.pars_inf);
    }

    #[test]
    fn count_variable_sites_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT", "ATTA", "ATGC", "ATGA"];
        let mat = get_matrix(&id, &seq);
        let mut site = Sites::new();
        let smat = site.index_sites(&mat);
        site.get_site_stats(&smat);
        assert_eq!(1, site.pars_inf);
    }

    #[test]
    fn count_parsimony_gap_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT---", "ATTA---", "ATGC---", "ATGA---"];
        let mat = get_matrix(&id, &seq);
        let mut site = Sites::new();
        let smat = site.index_sites(&mat);
        site.get_site_stats(&smat);
        assert_eq!(1, site.pars_inf);
        assert_eq!(3, site.variable);
    }

    #[test]
    fn get_site_stats_test() {
        let path = Path::new("test_files/concat.fasta");
        let input_format = SeqFormat::Fasta;
        let mut aln = Alignment::new();
        aln.get_aln_any(path, &input_format);
        let mut site = Sites::new();
        let smat = site.index_sites(&aln.alignment);
        site.get_site_stats(&smat);
        assert_eq!(18, site.conserved);
        assert_eq!(8, site.variable);
        assert_eq!(2, site.pars_inf);
    }
}
