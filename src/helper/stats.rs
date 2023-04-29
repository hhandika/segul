//! Statistics for the alignment and sequences.
use std::fmt::Debug;
use std::path::{Path, PathBuf};

use ahash::AHashMap as HashMap;

use crate::helper::types::{DataType, Header, SeqMatrix};

/// Get parsimony informative sites.
///
/// # Example
/// ```
/// use std::path::Path;
/// use segul::helper::types::{DataType, InputFmt};
/// use segul::helper::sequence::{SeqParser, SeqCheck};
/// use segul::helper::stats;
///
/// let file = Path::new("tests/files/concat.fasta");
/// let datatype = &DataType::Dna;
/// let input_fmt = &InputFmt::Fasta;
///
/// let seq = SeqParser::new(&file, datatype);
/// let (matrix, _) = seq.parse(input_fmt);
/// let pars_inf = stats::get_pars_inf(&matrix, datatype);
/// assert_eq!(pars_inf, 2);
/// ```
pub fn get_pars_inf(matrix: &SeqMatrix, datatype: &DataType) -> usize {
    Sites::default().get_pars_inf_only(matrix, datatype)
}

/// Get site summary statistics from a collection of alignments.
pub struct SiteSummary {
    /// General site summary
    pub total_loci: usize,
    pub total_sites: usize,
    pub min_sites: usize,
    pub max_sites: usize,
    pub mean_sites: f64,

    /// Conserved site summary
    pub cons_loci: usize,
    pub prop_cons_loci: f64,
    pub total_cons_site: usize,
    pub prop_cons_site: f64,
    pub min_cons_site: usize,
    pub max_cons_site: usize,
    pub mean_cons_site: f64,

    /// Variable site summary
    pub var_loci: usize,
    pub prop_var_loci: f64,
    pub total_var_site: usize,
    pub prop_var_site: f64,
    pub min_var_site: usize,
    pub max_var_site: usize,
    pub mean_var_site: f64,

    /// Parsimony inf site summary
    pub inf_loci: usize,
    pub prop_inf_loci: f64,
    pub total_inf_site: usize,
    pub prop_inf_site: f64,
    pub min_inf_site: usize,
    pub max_inf_site: usize,
    pub mean_inf_site: f64,
}

impl Default for SiteSummary {
    fn default() -> Self {
        Self::new()
    }
}

impl SiteSummary {
    /// Create a new `SiteSummary` instance.
    pub fn new() -> Self {
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

    /// Summarize all the sites from a collection of alignments.
    /// Returns a `SiteSummary` instance.
    pub fn summarize(&mut self, sites: &[Sites]) {
        self.total_loci = sites.len();
        self.total_sites = sites.iter().map(|s| s.counts).sum();
        self.min_sites = sites.iter().map(|s| s.counts).min().unwrap();
        self.max_sites = sites.iter().map(|s| s.counts).max().unwrap();
        self.mean_sites = self.total_sites as f64 / self.total_loci as f64;
        self.count_cons_sites(sites);
        self.count_var_sites(sites);
        self.count_inf_sites(sites);
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

/// A summary of the characters in a collection of alignments.
pub struct CharSummary {
    /// The minimum number of taxa in a collection of alignments.
    pub min_tax: usize,
    /// The maximum number of taxa in a collection of alignments.
    pub max_tax: usize,
    /// The mean number of taxa in a collection of alignments.
    pub mean_tax: f64,
    /// The GC content of the collection of alignments.
    pub gc_content: f64,
    /// The AT content of the collection of alignments.
    pub at_content: f64,
    /// The total number of missing data characters in the collection of alignments.
    pub missing_data: usize,
    /// The proportion of missing data characters in the collection of alignments.
    pub prop_missing_data: f64,
    /// The total number of characters in the collection of alignments.
    pub total_chars: usize,
    /// The total number of nucleotides in the collection of alignments.
    pub total_nucleotides: usize,
    /// A map of the characters in the collection of alignments.
    pub chars: HashMap<char, usize>,
}

impl Default for CharSummary {
    fn default() -> Self {
        Self::new()
    }
}

impl CharSummary {
    /// Create a new `CharSummary` instance.
    pub fn new() -> Self {
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
            chars: HashMap::new(),
        }
    }

    /// Summarize all the characters from a collection of alignments.
    /// Returns a `CharSummary` instance.
    /// # Arguments
    /// * `chars` - A collection of `CharMatrix` instances.
    /// * `datatype` - The datatype of the characters.
    pub fn summarize(&mut self, chars: &[CharMatrix], datatype: &DataType) {
        self.min_tax = chars.iter().map(|d| d.ntax).min().unwrap();
        self.max_tax = chars.iter().map(|d| d.ntax).max().unwrap();
        let sum_tax: usize = chars.iter().map(|d| d.ntax).sum();
        self.mean_tax = sum_tax as f64 / chars.len() as f64;
        self.total_chars = chars.iter().map(|d| d.chars.total_chars).sum();
        self.missing_data = chars.iter().map(|d| d.chars.missing_data).sum();
        self.count_chars(chars);
        if DataType::Dna == *datatype {
            self.total_nucleotides = chars.iter().map(|d| d.chars.nucleotides).sum();
            self.compute_gc_content(chars);
            self.compute_at_content(chars);
        }
        self.count_prop_missing_data();
    }

    fn count_chars(&mut self, chars: &[CharMatrix]) {
        chars
            .iter()
            .flat_map(|ch| ch.chars.chars.iter())
            .for_each(|(ch, count)| {
                *self.chars.entry(*ch).or_insert(0) += count;
            });
    }

    fn compute_gc_content(&mut self, chars: &[CharMatrix]) {
        let gc_count: usize = chars.iter().map(|c| c.chars.gc_count).sum();
        self.gc_content = gc_count as f64 / self.total_chars as f64;
    }

    fn compute_at_content(&mut self, chars: &[CharMatrix]) {
        let at_count: usize = chars.iter().map(|c| c.chars.at_count).sum();
        self.at_content = at_count as f64 / self.total_chars as f64;
    }

    #[inline]
    fn count_prop_missing_data(&mut self) {
        self.prop_missing_data = self.missing_data as f64 / self.total_chars as f64;
    }
}

/// Data completeness summary from a collection of alignments.
pub struct Completeness {
    pub completeness: Vec<(usize, usize)>,
    pub total_tax: usize,
    interval: usize,
}

impl Completeness {
    pub fn new(total_tax: &usize, interval: usize) -> Self {
        Self {
            completeness: Vec::new(),
            total_tax: *total_tax,
            interval,
        }
    }

    /// Compute the matrix completeness for a collection of alignments.
    /// Completeness is defined as the number of alignments with at least
    /// a given percentage of taxa.
    /// It will stop computing the completeness if
    /// it is equal to the total number of taxa.
    pub fn matrix_completeness(&mut self, chars: &[CharMatrix]) {
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

/// A site summary from an alignment.
#[derive(Debug, Clone)]
pub struct Sites {
    /// The path to the alignment.
    /// This is used to retrieve the alignment name.
    /// Used by the writer to sort the output.
    pub path: PathBuf,
    /// The number of conserved sites.
    pub conserved: usize,
    /// The number of variable sites.
    pub variable: usize,
    /// The number of parsimony informative sites.
    pub pars_inf: usize,
    /// The total number of sites.
    pub counts: usize,
    /// The proportion of variable sites.
    pub prop_var: f64,
    /// The proportion of conserved sites.
    pub prop_cons: f64,
    /// The proportion of parsimony informative sites.
    pub prop_pinf: f64,
}

impl Default for Sites {
    /// Create a new `Sites` instance.
    /// The path is set to an empty `PathBuf`.
    fn default() -> Self {
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
}

impl Sites {
    /// Create a new `Sites` instance.
    /// The path is set to the given `Path`.
    pub fn new(path: &Path) -> Self {
        Self {
            path: PathBuf::from(path),
            conserved: 0,
            variable: 0,
            pars_inf: 0,
            counts: 0,
            prop_var: 0.0,
            prop_cons: 0.0,
            prop_pinf: 0.0,
        }
    }

    /// Get the site statistics from an alignment.
    pub fn get_stats(&mut self, matrix: &SeqMatrix, datatype: &DataType) {
        let site_matrix = self.index_sites(matrix, datatype);
        self.get_site_stats(&site_matrix);
        self.count_sites();
        self.get_proportion();
    }

    fn get_pars_inf_only(&mut self, matrix: &SeqMatrix, datatype: &DataType) -> usize {
        let site_matrix = self.index_sites(matrix, datatype);
        self.get_site_stats(&site_matrix);
        self.pars_inf
    }

    fn index_sites(&self, matrix: &SeqMatrix, datatype: &DataType) -> HashMap<usize, Vec<u8>> {
        match datatype {
            DataType::Dna => self.index_site_dna(matrix),
            DataType::Aa => self.index_site_aa(matrix),
            _ => unreachable!(),
        }
    }

    fn index_site_dna(&self, matrix: &SeqMatrix) -> HashMap<usize, Vec<u8>> {
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

    fn index_site_aa(&self, matrix: &SeqMatrix) -> HashMap<usize, Vec<u8>> {
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

/// A struct to hold the character counts.
#[derive(Debug, Clone)]
pub struct CharMatrix {
    /// The number of taxa in the alignment.
    pub ntax: usize,
    /// The character counts.
    pub chars: Chars,
}

impl Default for CharMatrix {
    fn default() -> Self {
        Self::new()
    }
}

impl CharMatrix {
    /// Create a new `CharMatrix`.
    pub fn new() -> Self {
        Self {
            ntax: 0,
            chars: Chars::new(),
        }
    }

    /// Count the characters in an alignment.
    pub fn count_chars(&mut self, matrix: &SeqMatrix, header: &Header, datatype: &DataType) {
        self.ntax = header.ntax;
        self.chars.total_chars = header.nchar * self.ntax;
        self.parse_chars(matrix);
        if DataType::Dna == *datatype {
            self.chars.count_gc();
            self.chars.count_at();
            self.chars.count_nucleotides();
        }
        self.chars.count_missing_data();
        self.chars.calculate_prop_missing_data();
    }

    fn parse_chars(&mut self, matrix: &SeqMatrix) {
        matrix
            .values()
            .flat_map(|seqs| seqs.chars())
            .for_each(|ch| {
                *self.chars.chars.entry(ch.to_ascii_uppercase()).or_insert(0) += 1;
            });
    }
}

/// A struct to hold the character counts per taxon.
pub struct Taxa {
    /// The character counts per taxon.
    /// The key is the taxon name and the value is a `Chars` struct.
    pub records: HashMap<String, Chars>,
}

impl Default for Taxa {
    fn default() -> Self {
        Self::new()
    }
}

impl Taxa {
    /// Create a new `Taxa`.
    pub fn new() -> Self {
        Self {
            records: HashMap::new(),
        }
    }

    /// Summarize the character counts per taxon.
    ///
    /// # Arguments
    /// * `aln` - A sequence matrix.
    /// * `datatype` - The datatype of the alignment.
    pub fn summarize_taxa(&mut self, aln: &SeqMatrix, datatype: &DataType) {
        aln.iter().for_each(|(id, seq)| {
            let mut chars = Chars::new();
            // insert character to matrix
            seq.chars().for_each(|ch| {
                *chars.chars.entry(ch.to_ascii_uppercase()).or_insert(0) += 1;
            });
            chars.total_chars = seq.len();
            if DataType::Dna == *datatype {
                chars.count_gc();
                chars.count_at();
                chars.count_nucleotides();
            }
            chars.count_missing_data();
            self.records.insert(id.to_string(), chars);
        });
    }
}

/// A struct to hold the character counts per site.
/// A struct to hold the character counts.
#[derive(Debug, Clone)]
pub struct Chars {
    /// The character counts.
    /// The key is the character and the value is the count.
    pub chars: HashMap<char, usize>,
    /// The total number of characters in the alignment.
    pub total_chars: usize,
    /// The number of G and C characters in the alignment.
    pub gc_count: usize,
    /// The number of A and T characters in the alignment.
    pub at_count: usize,
    /// The number of nucleotides in the alignment.
    pub nucleotides: usize,
    /// The number of missing data characters in the alignment.
    /// This includes gaps (-) and missing data (?) characters.
    pub missing_data: usize,
    /// The proportion of missing data characters in the alignment.
    pub prop_missing_data: f64,
}

impl Default for Chars {
    fn default() -> Self {
        Self::new()
    }
}

impl Chars {
    /// Create a new `Chars`.
    pub fn new() -> Self {
        Self {
            chars: HashMap::new(),
            total_chars: 0,
            gc_count: 0,
            at_count: 0,
            nucleotides: 0,
            missing_data: 0,
            prop_missing_data: 0.0,
        }
    }

    /// Calculate GC count.
    pub fn count_gc(&mut self) {
        self.gc_count = self
            .chars
            .iter()
            .filter(|&(k, _)| *k == 'G' || *k == 'C')
            .map(|(_, v)| v)
            .sum();
    }

    /// Calculate AT count.
    pub fn count_at(&mut self) {
        self.at_count = self
            .chars
            .iter()
            .filter(|&(k, _)| *k == 'A' || *k == 'T')
            .map(|(_, v)| v)
            .sum();
    }

    /// Calculate nucleotide count.
    pub fn count_nucleotides(&mut self) {
        self.nucleotides = self.gc_count + self.at_count;
    }

    /// Calculate missing data count.
    pub fn count_missing_data(&mut self) {
        self.missing_data = self
            .chars
            .iter()
            .filter(|&(ch, _)| *ch == '-' || *ch == '?')
            .map(|(_, count)| count)
            .sum();
    }

    /// Calculate proportion of missing data.
    pub fn calculate_prop_missing_data(&mut self) {
        self.prop_missing_data = self.missing_data as f64 / self.total_chars as f64;
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::helper::sequence::SeqParser;
    use crate::helper::types::{DataType, InputFmt};
    use indexmap::IndexMap;

    const DNA: DataType = DataType::Dna;

    fn get_matrix(id: &[&str], seq: &[&str]) -> SeqMatrix {
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
        let pattern = Sites::default().get_patterns(site);
        let pattern_2 = Sites::default().get_patterns(site_2);
        assert_eq!(2, pattern);
        assert_eq!(3, pattern_2);
    }

    #[test]
    fn count_parsimony_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT", "ATTA", "ATGC", "ATGA"];
        let mat = get_matrix(&id, &seq);
        let mut site = Sites::default();
        let smat = site.index_sites(&mat, &DNA);
        site.get_site_stats(&smat);
        assert_eq!(1, site.pars_inf);
    }

    #[test]
    fn count_variable_sites_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT", "ATTA", "ATGC", "ATGA"];
        let mat = get_matrix(&id, &seq);
        let mut site = Sites::default();
        let site_mat = site.index_sites(&mat, &DNA);
        site.get_site_stats(&site_mat);
        assert_eq!(1, site.pars_inf);
    }

    #[test]
    fn count_parsimony_gap_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT---", "ATTA---", "ATGC---", "ATGA---"];
        let mat = get_matrix(&id, &seq);
        let mut site = Sites::default();
        let smat = site.index_sites(&mat, &DNA);
        site.get_site_stats(&smat);
        assert_eq!(1, site.pars_inf);
        assert_eq!(3, site.variable);
    }

    #[test]
    fn get_site_stats_test() {
        let path = Path::new("tests/files/concat.fasta");
        let input_format = InputFmt::Fasta;
        let aln = SeqParser::new(path, &DNA);
        let (matrix, _) = aln.get_alignment(&input_format);
        let mut site = Sites::default();
        let smat = site.index_sites(&matrix, &DNA);
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
    fn char_count_test() {
        let path = Path::new("tests/files/concat.fasta");
        let input_format = InputFmt::Fasta;
        let aln = SeqParser::new(path, &DNA);
        let (matrix, header) = aln.get_alignment(&input_format);
        let mut dna = CharMatrix::new();
        dna.count_chars(&matrix, &header, &DataType::Dna);
        assert_eq!(4, dna.ntax);
        assert_eq!(104, dna.chars.total_chars);
        assert_eq!(Some(&48), dna.chars.chars.get(&'A'));
        assert_eq!(Some(&22), dna.chars.chars.get(&'T'));
        assert_eq!(Some(&10), dna.chars.chars.get(&'G'));
        assert_eq!(None, dna.chars.chars.get(&'C'));
        assert_eq!(Some(&24), dna.chars.chars.get(&'?'));
        assert_eq!(24, dna.chars.missing_data);
        assert_eq!(None, dna.chars.chars.get(&'-'));
    }
}
