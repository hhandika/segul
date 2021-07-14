//! A module for sequence statistics.
use std::collections::HashMap;
use std::fmt::Debug;

use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;

use indexmap::IndexMap;
use rayon::prelude::*;

use crate::helper::alignment::Alignment;
use crate::helper::common::SeqFormat;
use crate::helper::finder::IDs;
use crate::helper::utils;
use crate::writer::sumwriter;

pub fn get_pars_inf(matrix: &IndexMap<String, String>) -> usize {
    Sites::new().get_pars_inf_only(matrix)
}

pub struct SeqStats<'a> {
    input_format: &'a SeqFormat,
    output: &'a str,
    ntax: usize,
    decrement: usize,
}

impl<'a> SeqStats<'a> {
    pub fn new(input_format: &'a SeqFormat, output: &'a str, decrement: usize) -> Self {
        Self {
            input_format,
            output,
            ntax: 0,
            decrement,
        }
    }

    pub fn get_seq_stats_file(&self, path: &Path) {
        let spin = utils::set_spinner();
        spin.set_message("Getting alignments...");
        let (site, dna) = self.get_stats(path);
        spin.finish_with_message("DONE!\n");
        sumwriter::CsvWriter::new(self.output)
            .write_summary_file(&site, &dna)
            .expect("CANNOT WRITE PER LOCUS SUMMARY STATS");
        sumwriter::print_stats(&site, &dna).unwrap();
    }

    pub fn get_stats_dir(&mut self, files: &[PathBuf]) {
        let spin = utils::set_spinner();
        spin.set_message("Indexing alignments...");
        self.get_ntax(files);
        spin.set_message("Getting summary stats...");
        let mut stats: Vec<(Sites, Dna)> = self.par_get_stats(files);
        stats.sort_by(|a, b| alphanumeric_sort::compare_path(&a.0.path, &b.0.path));
        let (sites, dna, complete) = self.get_summary_dna(&stats);
        sumwriter::CsvWriter::new(self.output)
            .write_summary_dir(&stats)
            .expect("CANNOT WRITE PER LOCUS SUMMARY STATS");
        let sum = sumwriter::SummaryWriter::new(&sites, &dna, &complete);
        sum.write_sum_to_file(self.output)
            .expect("CANNOT CREATE FILE FOR SUMMARY OUPUT");
        spin.finish_with_message("DONE!\n");
        sum.print_summary().expect("CANNOT WRITE SUMMARY TO STDOUT");
    }

    fn get_ntax(&mut self, files: &[PathBuf]) {
        self.ntax = IDs::new(files, self.input_format).get_id_all().len();
    }

    fn par_get_stats(&self, files: &[PathBuf]) -> Vec<(Sites, Dna)> {
        let (send, rec) = channel();
        files.par_iter().for_each_with(send, |s, file| {
            s.send(self.get_stats(file)).unwrap();
        });
        rec.iter().collect()
    }

    fn get_stats(&self, path: &Path) -> (Sites, Dna) {
        let mut aln = Alignment::new();
        aln.get_aln_any(path, self.input_format);
        let mut dna = Dna::new();
        dna.count_chars(&aln);
        let mut sites = Sites::new();
        sites.get_stats(path, &aln.alignment);

        (sites, dna)
    }

    fn get_summary_dna(&self, stats: &[(Sites, Dna)]) -> (SiteSummary, DnaSummary, Completeness) {
        let (sites, dna): (Vec<Sites>, Vec<Dna>) =
            stats.par_iter().map(|p| (p.0.clone(), p.1.clone())).unzip();
        let mut sum_sites = SiteSummary::new();
        sum_sites.get_summary(&sites);
        let mut sum_dna = DnaSummary::new();
        sum_dna.get_summary(&dna);
        let mut ntax_comp = Completeness::new(&self.ntax, self.decrement);
        ntax_comp.get_ntax_completeness(&dna);
        (sum_sites, sum_dna, ntax_comp)
    }
}

pub struct SiteSummary {
    // General site summary
    pub total_loci: usize,
    pub total_sites: usize,
    pub min_sites: usize,
    pub max_sites: usize,
    pub mean_sites: f64,

    // Conserved site summary
    pub cons_loci: usize,
    pub prop_cons_loci: f64,
    pub total_cons_site: usize,
    pub prop_cons_site: f64,
    pub min_cons_site: usize,
    pub max_cons_site: usize,
    pub mean_cons_site: f64,

    // Variable site summary
    pub var_loci: usize,
    pub prop_var_loci: f64,
    pub total_var_site: usize,
    pub prop_var_site: f64,
    pub min_var_site: usize,
    pub max_var_site: usize,
    pub mean_var_site: f64,

    // Parsimony inf site summary
    pub inf_loci: usize,
    pub prop_inf_loci: f64,
    pub total_inf_site: usize,
    pub prop_inf_site: f64,
    pub min_inf_site: usize,
    pub max_inf_site: usize,
    pub mean_inf_site: f64,
}

impl SiteSummary {
    fn new() -> Self {
        Self {
            total_sites: 0,
            total_loci: 0,
            min_sites: 0,
            max_sites: 0,
            mean_sites: 0.0,
            cons_loci: 0,
            prop_cons_loci: 0.0,
            total_cons_site: 0,
            prop_cons_site: 0.0,
            min_cons_site: 0,
            max_cons_site: 0,
            mean_cons_site: 0.0,
            var_loci: 0,
            prop_var_loci: 0.0,
            total_var_site: 0,
            prop_var_site: 0.0,
            min_var_site: 0,
            max_var_site: 0,
            mean_var_site: 0.0,
            inf_loci: 0,
            prop_inf_loci: 0.0,
            total_inf_site: 0,
            prop_inf_site: 0.0,
            min_inf_site: 0,
            max_inf_site: 0,
            mean_inf_site: 0.0,
        }
    }

    fn get_summary(&mut self, sites: &[Sites]) {
        self.total_loci = sites.len();
        self.total_sites = sites.iter().map(|s| s.counts).sum();
        self.min_sites = sites.iter().map(|s| s.counts).min().unwrap();
        self.max_sites = sites.iter().map(|s| s.counts).max().unwrap();
        self.mean_sites = self.total_sites as f64 / self.total_loci as f64;
        self.count_cons_sites(&sites);
        self.count_var_sites(&sites);
        self.count_inf_sites(&sites);
    }

    fn count_cons_sites(&mut self, sites: &[Sites]) {
        self.cons_loci = sites.iter().filter(|s| s.variable == 0).count();
        self.prop_cons_loci = self.cons_loci as f64 / self.total_loci as f64;
        self.total_cons_site = sites.iter().map(|s| s.conserved).sum();
        self.prop_cons_site = self.total_cons_site as f64 / self.total_sites as f64;
        self.min_cons_site = sites.iter().map(|s| s.conserved).min().unwrap();
        self.max_cons_site = sites.iter().map(|s| s.conserved).max().unwrap();
        self.mean_cons_site = self.total_cons_site as f64 / self.total_sites as f64;
    }

    fn count_var_sites(&mut self, sites: &[Sites]) {
        self.var_loci = sites.iter().filter(|s| s.variable > 0).count();
        self.prop_var_loci = self.var_loci as f64 / self.total_loci as f64;
        self.total_var_site = sites.iter().map(|s| s.variable).sum();
        self.prop_var_site = self.total_var_site as f64 / self.total_sites as f64;
        self.min_var_site = sites.iter().map(|s| s.variable).min().unwrap();
        self.max_var_site = sites.iter().map(|s| s.variable).max().unwrap();
        self.mean_var_site = self.total_var_site as f64 / self.total_sites as f64;
    }

    fn count_inf_sites(&mut self, sites: &[Sites]) {
        self.inf_loci = sites.iter().filter(|s| s.pars_inf > 0).count();
        self.prop_inf_loci = self.inf_loci as f64 / self.total_loci as f64;
        self.total_inf_site = sites.iter().map(|s| s.pars_inf).sum();
        self.prop_inf_site = self.total_inf_site as f64 / self.total_sites as f64;
        self.min_inf_site = sites.iter().map(|s| s.pars_inf).min().unwrap();
        self.max_inf_site = sites.iter().map(|s| s.pars_inf).max().unwrap();
        self.mean_inf_site = self.total_inf_site as f64 / self.total_sites as f64;
    }
}

pub struct DnaSummary {
    pub min_tax: usize,
    pub max_tax: usize,
    pub mean_tax: f64,
    pub gc_content: f64,
    pub at_content: f64,
    pub missing_data: usize,
    pub prop_missing_data: f64,
    pub total_chars: usize,
    pub total_nucleotides: usize,
    pub total_a: usize,
    pub total_c: usize,
    pub total_g: usize,
    pub total_t: usize,
    pub total_n: usize,
    pub total_missings: usize,
    pub total_gaps: usize,
    pub total_undetermined: usize,
}

impl DnaSummary {
    fn new() -> Self {
        Self {
            total_chars: 0,
            min_tax: 0,
            max_tax: 0,
            mean_tax: 0.0,
            gc_content: 0.0,
            at_content: 0.0,
            missing_data: 0,
            prop_missing_data: 0.0,
            total_nucleotides: 0,
            total_a: 0,
            total_c: 0,
            total_g: 0,
            total_t: 0,
            total_n: 0,
            total_missings: 0,
            total_gaps: 0,
            total_undetermined: 0,
        }
    }

    fn get_summary(&mut self, dna: &[Dna]) {
        self.min_tax = dna.iter().map(|d| d.ntax).min().unwrap();
        self.max_tax = dna.iter().map(|d| d.ntax).max().unwrap();
        let sum_tax: usize = dna.iter().map(|d| d.ntax).sum();
        self.mean_tax = sum_tax as f64 / dna.len() as f64;
        self.total_chars = dna.iter().map(|d| d.total_chars).sum();
        self.count_chars(dna);
        self.get_total_nucleotides();
        self.count_gc_at_content();
        self.count_missing_data();
    }

    fn count_chars(&mut self, dna: &[Dna]) {
        self.total_a = dna.iter().map(|d| d.a_count).sum();
        self.total_t = dna.iter().map(|d| d.t_count).sum();
        self.total_g = dna.iter().map(|d| d.g_count).sum();
        self.total_c = dna.iter().map(|d| d.c_count).sum();
        self.total_n = dna.iter().map(|d| d.n_count).sum();
        self.total_missings = dna.iter().map(|d| d.missings).sum();
        self.total_gaps = dna.iter().map(|d| d.gaps).sum();
        self.total_undetermined = dna.iter().map(|d| d.undetermined).sum();
    }

    fn get_total_nucleotides(&mut self) {
        self.total_nucleotides = self.total_a + self.total_t + self.total_g + self.total_c
    }

    fn count_gc_at_content(&mut self) {
        self.gc_content = (self.total_g + self.total_c) as f64 / self.total_chars as f64;
        self.at_content = (self.total_g + self.total_c) as f64 / self.total_chars as f64;
    }

    fn count_missing_data(&mut self) {
        self.missing_data = self.total_missings + self.total_gaps + self.total_n;
        self.prop_missing_data = self.missing_data as f64 / self.total_chars as f64;
    }
}

pub struct Completeness {
    pub completeness: Vec<(usize, usize)>,
    pub total_tax: usize,
    decrement: usize,
}

impl Completeness {
    fn new(total_tax: &usize, decrement: usize) -> Self {
        Self {
            completeness: Vec::new(),
            total_tax: *total_tax,
            decrement,
        }
    }

    fn get_ntax_completeness(&mut self, dna: &[Dna]) {
        let ntax: Vec<usize> = dna.iter().map(|d| d.ntax).collect();
        let mut values: usize = 100;

        while values > 0 {
            let percent = values as f64 / 100.0;
            let ntax_comp = self.count_min_tax(&ntax, percent);
            self.completeness.push((values, ntax_comp));
            if ntax_comp == ntax.len() {
                break;
            } else {
                values -= self.decrement;
            }
        }
    }

    fn count_min_tax(&self, ntax: &[usize], percent: f64) -> usize {
        ntax.iter()
            .filter(|&n| n >= &self.compute_min_taxa(percent))
            .count()
    }

    fn compute_min_taxa(&self, percent: f64) -> usize {
        (self.total_tax as f64 * percent).floor() as usize
    }
}

#[derive(Debug, Clone)]
pub struct Sites {
    pub path: PathBuf,
    pub conserved: usize,
    pub variable: usize,
    pub pars_inf: usize,
    pub counts: usize,
    pub prop_var: f64,
    pub prop_cons: f64,
    pub prop_pinf: f64,
}

impl Sites {
    fn new() -> Self {
        Self {
            path: PathBuf::new(),
            conserved: 0,
            variable: 0,
            pars_inf: 0,
            counts: 0,
            prop_var: 0.0,
            prop_cons: 0.0,
            prop_pinf: 0.0,
        }
    }

    fn get_stats(&mut self, path: &Path, matrix: &IndexMap<String, String>) {
        self.path = path.to_path_buf();
        let site_matrix = self.index_sites(matrix);
        self.get_site_stats(&site_matrix);
        self.count_sites();
        self.get_proportion();
    }

    fn get_pars_inf_only(&mut self, matrix: &IndexMap<String, String>) -> usize {
        let site_matrix = self.index_sites(matrix);
        self.get_site_stats(&site_matrix);
        self.pars_inf
    }

    fn index_sites(&mut self, matrix: &IndexMap<String, String>) -> HashMap<usize, Vec<u8>> {
        let mut site_matrix: HashMap<usize, Vec<u8>> = HashMap::new();
        matrix.values().for_each(|seq| {
            seq.bytes().enumerate().for_each(|(idx, dna)| {
                match site_matrix.get_mut(&idx) {
                    Some(value) => match dna {
                        b'a' | b'g' | b't' | b'c' | b'A' | b'G' | b'T' | b'C' => value.push(dna),
                        _ => (), // ignore ambiguous characters
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
        if uniques.len() > 1 {
            self.variable += 1;
            uniques.iter().for_each(|&ch| {
                let patterns = bytecount::count(site, ch);
                if patterns > 1 {
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

    fn get_proportion(&mut self) {
        self.prop_cons = self.conserved as f64 / self.counts as f64;
        self.prop_var = self.variable as f64 / self.counts as f64;
        self.prop_pinf = self.pars_inf as f64 / self.counts as f64;
    }
}

#[derive(Debug, Clone)]
pub struct Dna {
    pub a_count: usize,
    pub c_count: usize,
    pub g_count: usize,
    pub t_count: usize,
    pub n_count: usize,
    pub missings: usize,
    pub gaps: usize,
    pub undetermined: usize,
    pub total_chars: usize,
    pub ntax: usize,
    pub missing_data: usize,
    pub prop_missing_data: f64,
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
            undetermined: 0,
            total_chars: 0,
            ntax: 0,
            missing_data: 0,
            prop_missing_data: 0.0,
        }
    }

    fn count_chars(&mut self, aln: &Alignment) {
        self.ntax = aln.header.ntax;
        self.total_chars = aln.header.nchar * self.ntax;
        aln.alignment.values().for_each(|seqs| {
            seqs.bytes().for_each(|ch| match ch {
                b'a' | b'A' => self.a_count += 1,
                b'c' | b'C' => self.c_count += 1,
                b'g' | b'G' => self.g_count += 1,
                b't' | b'T' => self.t_count += 1,
                b'n' | b'N' => self.n_count += 1,
                b'?' | b'.' | b'~' => self.missings += 1,
                b'O' | b'o' | b'X' | b'x' => self.missings += 1, // Following iqtree treatments
                b'-' => self.gaps += 1,
                _ => self.undetermined += 1,
            })
        });

        self.count_missing_data();
    }

    fn count_missing_data(&mut self) {
        self.missing_data = self.missings + self.gaps + self.n_count;
        self.prop_missing_data = self.missing_data as f64 / self.total_chars as f64;
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

    #[test]
    fn filter_min_tax_test() {
        let ntax = vec![10, 8, 20, 30, 60];
        let comp = Completeness::new(&60, 2);
        assert_eq!(2, comp.count_min_tax(&ntax, 0.5))
    }
}
