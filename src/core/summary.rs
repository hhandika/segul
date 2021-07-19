//! A module for sequence statistics.
use std::collections::{BTreeMap, HashMap};
use std::fmt::Debug;

use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;

use indexmap::IndexMap;
use rayon::prelude::*;

use crate::helper::common::{DataType, InputFmt};
use crate::helper::finder::IDs;
use crate::helper::sequence::Sequence;
use crate::helper::utils;
use crate::writer::sumwriter;

pub fn get_pars_inf(matrix: &IndexMap<String, String>, datatype: &DataType) -> usize {
    Sites::new().get_pars_inf_only(matrix, datatype)
}

pub struct SeqStats<'a> {
    input_format: &'a InputFmt,
    output: &'a str,
    datatype: &'a DataType,
    ntax: usize,
    interval: usize,
}

impl<'a> SeqStats<'a> {
    pub fn new(
        input_format: &'a InputFmt,
        output: &'a str,
        interval: usize,
        datatype: &'a DataType,
    ) -> Self {
        Self {
            input_format,
            output,
            ntax: 0,
            interval,
            datatype,
        }
    }

    pub fn get_seq_stats_file(&self, path: &Path) {
        let spin = utils::set_spinner();
        spin.set_message("Getting alignments...");
        let (site, dna) = self.get_stats(path);
        spin.finish_with_message("DONE!\n");
        sumwriter::CsvWriter::new(self.output, self.datatype)
            .write_summary_file(&site, &dna)
            .expect("CANNOT WRITE PER LOCUS SUMMARY STATS");
        sumwriter::print_stats(&site, &dna).unwrap();
    }

    pub fn get_stats_dir(&mut self, files: &[PathBuf]) {
        let spin = utils::set_spinner();
        spin.set_message("Indexing alignments...");
        self.get_ntax(files);

        spin.set_message("Getting summary stats...");
        let mut stats: Vec<(Sites, Chars)> = self.par_get_stats(files);
        stats.sort_by(|a, b| alphanumeric_sort::compare_path(&a.0.path, &b.0.path));
        let (sites, dna, complete) = self.get_summary_dna(&stats);

        spin.set_message("Writing output files...");
        sumwriter::CsvWriter::new(self.output, self.datatype)
            .write_summary_dir(&stats)
            .expect("Failed writing a per locus csv file");

        let sum = sumwriter::SummaryWriter::new(&sites, &dna, &complete, self.datatype);
        sum.write_sum_to_file(self.output)
            .expect("File writing a summary statistics file");
        spin.finish_with_message("DONE!\n");
        sum.print_summary().expect("Failed writing to stdout");
    }

    fn get_ntax(&mut self, files: &[PathBuf]) {
        self.ntax = IDs::new(files, self.input_format, self.datatype)
            .get_id_all()
            .len();
    }

    fn par_get_stats(&self, files: &[PathBuf]) -> Vec<(Sites, Chars)> {
        let (send, rec) = channel();
        files.par_iter().for_each_with(send, |s, file| {
            s.send(self.get_stats(file)).unwrap();
        });
        rec.iter().collect()
    }

    fn get_stats(&self, path: &Path) -> (Sites, Chars) {
        let mut aln = Sequence::new();
        aln.get_alignment(path, self.input_format, self.datatype);
        let mut dna = Chars::new();
        dna.count_chars(&aln);
        let mut sites = Sites::new();
        sites.get_stats(path, &aln.matrix, self.datatype);

        (sites, dna)
    }

    fn get_summary_dna(
        &self,
        stats: &[(Sites, Chars)],
    ) -> (SiteSummary, CharSummary, Completeness) {
        let (sites, dna): (Vec<Sites>, Vec<Chars>) =
            stats.par_iter().map(|p| (p.0.clone(), p.1.clone())).unzip();
        let mut sum_sites = SiteSummary::new();
        sum_sites.get_summary(&sites);
        let mut sum_dna = CharSummary::new();
        sum_dna.get_summary(&dna);
        let mut mat_comp = Completeness::new(&self.ntax, self.interval);
        mat_comp.matrix_completeness(&dna);
        (sum_sites, sum_dna, mat_comp)
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

pub struct CharSummary {
    pub min_tax: usize,
    pub max_tax: usize,
    pub mean_tax: f64,
    pub gc_content: f64,
    pub at_content: f64,
    pub missing_data: usize,
    pub prop_missing_data: f64,
    pub total_chars: usize,
    pub total_nucleotides: usize,
    pub chars: BTreeMap<char, usize>,
}

impl CharSummary {
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
            chars: BTreeMap::new(),
        }
    }

    fn get_summary(&mut self, chars: &[Chars]) {
        self.min_tax = chars.iter().map(|d| d.ntax).min().unwrap();
        self.max_tax = chars.iter().map(|d| d.ntax).max().unwrap();
        let sum_tax: usize = chars.iter().map(|d| d.ntax).sum();
        self.mean_tax = sum_tax as f64 / chars.len() as f64;
        self.total_chars = chars.iter().map(|d| d.total_chars).sum();
        self.count_chars(chars);
        self.count_nucleotides();
        self.compute_gc_content(chars);
        self.compute_at_content(chars);
        self.count_missing_data();
    }

    fn count_chars(&mut self, chars: &[Chars]) {
        chars
            .iter()
            .flat_map(|ch| ch.chars.iter())
            .for_each(|(ch, count)| {
                *self.chars.entry(ch.to_ascii_uppercase()).or_insert(0) += count;
            });
    }

    fn count_nucleotides(&mut self) {
        self.total_nucleotides = self
            .chars
            .iter()
            .filter(|&(k, _)| *k == 'A' || *k == 'C' || *k == 'G' || *k == 'T')
            .map(|(_, count)| count)
            .sum();
    }

    fn compute_gc_content(&mut self, chars: &[Chars]) {
        let gc_count: usize = chars.iter().map(|c| c.gc_count).sum();
        self.gc_content = gc_count as f64 / self.total_chars as f64;
    }

    fn compute_at_content(&mut self, chars: &[Chars]) {
        let at_count: usize = chars.iter().map(|c| c.at_count).sum();
        self.at_content = at_count as f64 / self.total_chars as f64;
    }

    fn count_missing_data(&mut self) {
        self.missing_data = self
            .chars
            .iter()
            .filter(|&(k, _)| *k == '-' || *k == '?' || *k == 'N')
            .map(|(_, count)| count)
            .sum();
        self.prop_missing_data = self.missing_data as f64 / self.total_chars as f64;
    }
}

pub struct Completeness {
    pub completeness: Vec<(usize, usize)>,
    pub total_tax: usize,
    interval: usize,
}

impl Completeness {
    fn new(total_tax: &usize, interval: usize) -> Self {
        Self {
            completeness: Vec::new(),
            total_tax: *total_tax,
            interval,
        }
    }

    fn matrix_completeness(&mut self, chars: &[Chars]) {
        let ntax: Vec<usize> = chars.iter().map(|d| d.ntax).collect();
        let mut values: usize = 100;

        while values > 0 {
            let percent = values as f64 / 100.0;
            let mat_comp = self.count_mat_completeness(&ntax, percent);
            self.completeness.push((values, mat_comp));
            if mat_comp == ntax.len() {
                break;
            } else {
                values -= self.interval;
            }
        }
    }

    fn count_mat_completeness(&self, ntax: &[usize], percent: f64) -> usize {
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

    fn get_stats(&mut self, path: &Path, matrix: &IndexMap<String, String>, datatype: &DataType) {
        self.path = path.to_path_buf();
        let site_matrix = self.index_sites(matrix, datatype);
        self.get_site_stats(&site_matrix);
        self.count_sites();
        self.get_proportion();
    }

    fn get_pars_inf_only(
        &mut self,
        matrix: &IndexMap<String, String>,
        datatype: &DataType,
    ) -> usize {
        let site_matrix = self.index_sites(matrix, datatype);
        self.get_site_stats(&site_matrix);
        self.pars_inf
    }

    fn index_sites(
        &self,
        matrix: &IndexMap<String, String>,
        datatype: &DataType,
    ) -> HashMap<usize, Vec<u8>> {
        match datatype {
            DataType::Dna => self.index_site_dna(matrix),
            DataType::Aa => self.index_site_aa(matrix),
            _ => unreachable!(),
        }
    }

    fn index_site_dna(&self, matrix: &IndexMap<String, String>) -> HashMap<usize, Vec<u8>> {
        let mut site_matrix: HashMap<usize, Vec<u8>> = HashMap::new();
        matrix.values().for_each(|seq| {
            seq.bytes()
                .enumerate()
                .for_each(|(idx, dna)| match site_matrix.get_mut(&idx) {
                    Some(value) => {
                        // ignore ambiguous characters
                        if self.is_non_ambiguous_dna(&dna) {
                            value.push(dna);
                        }
                    }
                    None => {
                        if self.is_non_ambiguous_dna(&dna) {
                            site_matrix.insert(idx, vec![dna]);
                        }
                    }
                })
        });

        site_matrix
    }

    fn index_site_aa(&self, matrix: &IndexMap<String, String>) -> HashMap<usize, Vec<u8>> {
        let mut site_matrix: HashMap<usize, Vec<u8>> = HashMap::new();
        matrix.values().for_each(|seq| {
            seq.bytes().enumerate().for_each(|(idx, aa)| {
                match site_matrix.get_mut(&idx) {
                    Some(value) => {
                        // ignore ambiguous characters
                        if !self.is_ambiguous_aa(&aa) {
                            value.push(aa)
                        }
                    }
                    None => {
                        if !self.is_ambiguous_aa(&aa) {
                            site_matrix.insert(idx, vec![aa]);
                        }
                    }
                }
            })
        });

        site_matrix
    }

    fn is_non_ambiguous_dna(&self, ch: &u8) -> bool {
        let non_ambiguous_dna = b"acgtACGT";
        non_ambiguous_dna.contains(ch)
    }

    // We match ambiguous for aa because it is shorter
    // then matching non-ambiguous aa
    fn is_ambiguous_aa(&self, ch: &u8) -> bool {
        let ambiguous_aa = b"XBZJU?-.~*";
        ambiguous_aa.contains(ch)
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
pub struct Chars {
    pub total_chars: usize,
    pub ntax: usize,
    pub chars: BTreeMap<char, usize>,
    pub gc_count: usize,
    pub at_count: usize,
    pub missing_data: usize,
    pub prop_missing_data: f64,
}

impl Chars {
    fn new() -> Self {
        Self {
            total_chars: 0,
            ntax: 0,
            chars: BTreeMap::new(),
            gc_count: 0,
            at_count: 0,
            missing_data: 0,
            prop_missing_data: 0.0,
        }
    }

    fn count_chars(&mut self, aln: &Sequence) {
        self.ntax = aln.header.ntax;
        self.total_chars = aln.header.nchar * self.ntax;
        aln.matrix
            .values()
            .flat_map(|seqs| seqs.chars())
            .for_each(|ch| {
                *self.chars.entry(ch.to_ascii_uppercase()).or_insert(0) += 1;
            });
        self.count_gc();
        self.count_at();
        self.count_missing_data();
    }

    fn count_gc(&mut self) {
        self.gc_count = self
            .chars
            .iter()
            .filter(|&(k, _)| *k == 'G' || *k == 'C')
            .map(|(_, v)| v)
            .sum();
    }

    fn count_at(&mut self) {
        self.at_count = self
            .chars
            .iter()
            .filter(|&(k, _)| *k == 'A' || *k == 'T')
            .map(|(_, v)| v)
            .sum();
    }

    fn count_missing_data(&mut self) {
        self.missing_data = self
            .chars
            .iter()
            .filter(|&(ch, _)| *ch == '-' || *ch == '?' || *ch == 'N' || *ch == 'n')
            .map(|(_, count)| count)
            .sum();
        self.prop_missing_data = self.missing_data as f64 / self.total_chars as f64;
    }
}

#[cfg(test)]
mod test {
    use super::*;

    const DNA: DataType = DataType::Dna;

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
        let smat = site.index_sites(&mat, &DNA);
        site.get_site_stats(&smat);
        assert_eq!(1, site.pars_inf);
    }

    #[test]
    fn count_variable_sites_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT", "ATTA", "ATGC", "ATGA"];
        let mat = get_matrix(&id, &seq);
        let mut site = Sites::new();
        let smat = site.index_sites(&mat, &DNA);
        site.get_site_stats(&smat);
        assert_eq!(1, site.pars_inf);
    }

    #[test]
    fn count_parsimony_gap_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT---", "ATTA---", "ATGC---", "ATGA---"];
        let mat = get_matrix(&id, &seq);
        let mut site = Sites::new();
        let smat = site.index_sites(&mat, &DNA);
        site.get_site_stats(&smat);
        assert_eq!(1, site.pars_inf);
        assert_eq!(3, site.variable);
    }

    #[test]
    fn get_site_stats_test() {
        let path = Path::new("test_files/concat.fasta");
        let input_format = InputFmt::Fasta;
        let mut aln = Sequence::new();
        aln.get_alignment(path, &input_format, &DNA);
        let mut site = Sites::new();
        let smat = site.index_sites(&aln.matrix, &DNA);
        site.get_site_stats(&smat);
        assert_eq!(18, site.conserved);
        assert_eq!(8, site.variable);
        assert_eq!(2, site.pars_inf);
    }

    #[test]
    fn filter_min_tax_test() {
        let ntax = vec![10, 8, 20, 30, 60];
        let comp = Completeness::new(&60, 2);
        assert_eq!(2, comp.count_mat_completeness(&ntax, 0.5))
    }

    #[test]
    fn dna_count_test() {
        let path = Path::new("test_files/concat.fasta");
        let input_format = InputFmt::Fasta;
        let mut aln = Sequence::new();
        aln.get_alignment(path, &input_format, &DNA);
        let mut dna = Chars::new();
        dna.count_chars(&aln);
        assert_eq!(4, dna.ntax);
        assert_eq!(104, dna.total_chars);
        assert_eq!(Some(&48), dna.chars.get(&'A'));
        assert_eq!(Some(&22), dna.chars.get(&'T'));
        assert_eq!(Some(&10), dna.chars.get(&'G'));
        assert_eq!(None, dna.chars.get(&'C'));
        assert_eq!(Some(&24), dna.chars.get(&'?'));
        assert_eq!(24, dna.missing_data);
        assert_eq!(None, dna.chars.get(&'-'));
    }
}
