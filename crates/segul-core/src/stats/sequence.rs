//! Statistics for the alignment and sequences.
use std::fmt::Debug;
use std::path::{Path, PathBuf};

use ahash::AHashMap as HashMap;

use crate::helper::types::{DataType, Header, SeqMatrix};

pub type SiteMap = HashMap<usize, Vec<u8>>;

const SITE_BLOCK: usize = 4096;

const DNA_STATE_BIT: [u32; 256] = build_dna_state_bits();
const AA_STATE_BIT: [u32; 256] = build_aa_state_bits();

const fn build_dna_state_bits() -> [u32; 256] {
    let mut table = [0u32; 256];
    table[b'A' as usize] = 1 << 0;
    table[b'a' as usize] = 1 << 0;
    table[b'C' as usize] = 1 << 1;
    table[b'c' as usize] = 1 << 1;
    table[b'G' as usize] = 1 << 2;
    table[b'g' as usize] = 1 << 2;
    table[b'T' as usize] = 1 << 3;
    table[b't' as usize] = 1 << 3;
    table
}

const fn build_aa_state_bits() -> [u32; 256] {
    let mut table = [0u32; 256];
    let mut b = b'A';
    while b <= b'Z' {
        let ambiguous = matches!(b, b'X' | b'B' | b'Z' | b'J' | b'U');
        if !ambiguous {
            let bit = 1u32 << (b - b'A');
            table[b as usize] = bit;
            table[b.to_ascii_lowercase() as usize] = bit;
        }
        b += 1;
    }
    table
}

fn state_bit_table(datatype: &DataType) -> &'static [u32; 256] {
    match datatype {
        DataType::Dna => &DNA_STATE_BIT,
        DataType::Aa => &AA_STATE_BIT,
        _ => unreachable!(),
    }
}

#[inline]
fn scan_site_block(
    rows: &[&[u8]],
    start: usize,
    width: usize,
    table: &[u32; 256],
    seen1: &mut [u32; SITE_BLOCK],
    seen2: &mut [u32; SITE_BLOCK],
) {
    seen1[..width].fill(0);
    seen2[..width].fill(0);
    for row in rows {
        let n = row.len().saturating_sub(start).min(width);
        for (i, &b) in row[start..start + n].iter().enumerate() {
            let bit = table[b as usize];
            seen2[i] |= seen1[i] & bit;
            seen1[i] |= bit;
        }
    }
}

/// Get parsimony informative sites.
///
/// # Example
/// ```
/// use std::path::Path;
/// use segul_core::helper::types::{DataType, InputFmt};
/// use segul_core::helper::sequence::{SeqParser, SeqCheck};
/// use segul_core::stats::sequence as stats;
///
/// let file = Path::new("../../tests/files/concat.fasta");
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
        self.scan_sites(matrix, datatype, |_, _| {});
        self.count_sites();
        self.get_proportion();
    }

    /// Get parsimony informative for each site
    /// Return the site location and the number of parsimony informative sites.
    pub fn get_pars_inf_per_site(
        &mut self,
        matrix: &SeqMatrix,
        datatype: &DataType,
    ) -> Vec<(usize, usize)> {
        let mut pars_inf_per_site: Vec<(usize, usize)> = Vec::new();
        self.scan_sites(matrix, datatype, |site, n_patterns| {
            pars_inf_per_site.push((site, n_patterns));
        });
        pars_inf_per_site
    }

    /// Get sites with missing data less than a threshold.
    /// Return the site location
    pub fn get_site_without_missing_data(&self, matrix: &SeqMatrix, threshold: f64) -> Vec<usize> {
        let rows: Vec<&[u8]> = matrix.values().collect();
        let ntax = rows.len();
        if ntax == 0 {
            return Vec::new();
        }
        let nchar = rows.iter().map(|r| r.len()).max().unwrap_or(0);
        let mut sites = Vec::with_capacity(nchar);
        let mut missing = [0u32; SITE_BLOCK];
        for start in (0..nchar).step_by(SITE_BLOCK) {
            let width = (nchar - start).min(SITE_BLOCK);
            missing[..width].fill(0);
            for row in &rows {
                let n = row.len().saturating_sub(start).min(width);
                for (i, &b) in row[start..start + n].iter().enumerate() {
                    missing[i] += u32::from(b == b'-' || b == b'?');
                }
            }
            for (i, &count) in missing[..width].iter().enumerate() {
                if (count as f64 / ntax as f64) < threshold {
                    sites.push(start + i);
                }
            }
        }
        sites
    }

    /// Get sites with parsimony informative greater than a threshold.
    /// Return the site location and the number of parsimony informative sites.
    pub fn get_site_with_pars_informative(
        &mut self,
        matrix: &SeqMatrix,
        datatype: &DataType,
        threshold: usize,
    ) -> Vec<(usize, usize)> {
        let mut sites = Vec::new();
        self.scan_sites(matrix, datatype, |site, n_patterns| {
            if n_patterns >= threshold {
                sites.push((site, n_patterns));
            }
        });
        sites
    }

    /// Get missing data per site
    /// Return the site location and the proportion of missing data.
    pub fn get_missing_data_per_site(
        &self,
        indexed_sites: &SiteMap,
        site_count: usize,
    ) -> Vec<(usize, f64)> {
        let mut missing_data_per_site: Vec<(usize, f64)> = Vec::new();
        indexed_sites.iter().for_each(|(site, seq)| {
            let n_missing = seq.iter().filter(|&ch| *ch == b'-' || *ch == b'?').count();
            missing_data_per_site.push((*site, n_missing as f64 / site_count as f64));
        });
        missing_data_per_site
    }

    /// Get the sites as is
    /// Return the site location including missing data characters.
    pub fn index_site_with_missing_data(&self, matrix: &SeqMatrix) -> SiteMap {
        let mut site_matrix: SiteMap = HashMap::new();
        matrix.values().for_each(|seq| {
            seq.iter()
                .copied()
                .enumerate()
                .for_each(|(idx, dna)| match site_matrix.get_mut(&idx) {
                    Some(value) => value.push(dna),
                    None => {
                        site_matrix.insert(idx, vec![dna]);
                    }
                })
        });

        site_matrix
    }

    pub fn index_site_without_ambiguity(&self, matrix: &SeqMatrix, datatype: &DataType) -> SiteMap {
        match datatype {
            DataType::Dna => self.index_site_dna(matrix),
            DataType::Aa => self.index_site_aa(matrix),
            _ => unreachable!(),
        }
    }

    fn get_pars_inf_only(&mut self, matrix: &SeqMatrix, datatype: &DataType) -> usize {
        self.scan_sites(matrix, datatype, |_, _| {});
        self.pars_inf
    }

    fn scan_sites<F>(&mut self, matrix: &SeqMatrix, datatype: &DataType, mut on_pis: F)
    where
        F: FnMut(usize, usize),
    {
        let table = state_bit_table(datatype);
        let rows: Vec<&[u8]> = matrix.values().collect();
        let nchar = rows.iter().map(|r| r.len()).max().unwrap_or(0);
        let mut seen1 = [0u32; SITE_BLOCK];
        let mut seen2 = [0u32; SITE_BLOCK];
        for start in (0..nchar).step_by(SITE_BLOCK) {
            let width = (nchar - start).min(SITE_BLOCK);
            scan_site_block(&rows, start, width, table, &mut seen1, &mut seen2);
            for i in 0..width {
                let s1 = seen1[i];
                let n_states = s1.count_ones();
                if n_states == 0 {
                    continue;
                }
                if n_states == 1 {
                    self.conserved += 1;
                } else {
                    self.variable += 1;
                    let n_patterns = seen2[i].count_ones() as usize;
                    if n_patterns >= 2 {
                        self.pars_inf += 1;
                        on_pis(start + i, n_patterns);
                    }
                }
            }
        }
    }

    fn index_site_dna(&self, matrix: &SeqMatrix) -> SiteMap {
        let mut site_matrix: SiteMap = HashMap::new();
        matrix.values().for_each(|seq| {
            seq.iter().copied().enumerate().for_each(|(idx, dna)| {
                if self.is_non_ambiguous_dna(&dna) {
                    let base = dna.to_ascii_uppercase();
                    match site_matrix.get_mut(&idx) {
                        Some(value) => value.push(base),
                        None => {
                            site_matrix.insert(idx, vec![base]);
                        }
                    }
                }
            })
        });

        site_matrix
    }

    fn index_site_aa(&self, matrix: &SeqMatrix) -> SiteMap {
        let mut site_matrix: SiteMap = HashMap::new();
        matrix.values().for_each(|seq| {
            seq.iter().copied().enumerate().for_each(|(idx, aa)| {
                if !self.is_ambiguous_aa(&aa) {
                    let res = aa.to_ascii_uppercase();
                    match site_matrix.get_mut(&idx) {
                        Some(value) => value.push(res),
                        None => {
                            site_matrix.insert(idx, vec![res]);
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
        ambiguous_aa.contains(&ch.to_ascii_uppercase())
    }

    #[cfg(test)]
    fn compute_parsimony_informative(&mut self, site_matrix: &SiteMap) {
        site_matrix.values().for_each(|site| {
            let n_patterns = self.get_patterns(site);
            if n_patterns >= 2 {
                self.pars_inf += 1;
            }
        });
    }

    #[cfg(test)]
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

    /// Derive alignment character counts from already-computed per-taxon counts.
    pub fn count_from_taxa(&mut self, taxa: &Taxa, header: &Header, datatype: &DataType) {
        self.ntax = header.ntax;
        self.chars.total_chars = header.nchar * self.ntax;
        for t_chars in taxa.records.values() {
            for (ch, count) in &t_chars.chars {
                *self.chars.chars.entry(*ch).or_insert(0) += count;
            }
            self.chars.missing_data += t_chars.missing_data;
            if DataType::Dna == *datatype {
                self.chars.gc_count += t_chars.gc_count;
                self.chars.at_count += t_chars.at_count;
                self.chars.nucleotides += t_chars.nucleotides;
            }
        }
        self.chars.calculate_prop_missing_data();
    }

    fn parse_chars(&mut self, matrix: &SeqMatrix) {
        let mut counts = [[0usize; 256]; 4];
        count_bytes_4lane(matrix.as_bytes(), &mut counts);
        self.chars.chars = histogram_to_chars(&counts);
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
        let mut counts = [[0usize; 256]; 4];
        aln.iter().for_each(|(id, seq)| {
            counts.iter_mut().for_each(|lane| lane.fill(0));
            count_bytes_4lane(seq, &mut counts);
            let mut chars = Chars::new();
            chars.chars = histogram_to_chars(&counts);
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

#[inline]
fn count_bytes_4lane(seq: &[u8], counts: &mut [[usize; 256]; 4]) {
    let (chunks, remainder) = seq.as_chunks::<4>();
    for c in chunks {
        counts[0][c[0] as usize] += 1;
        counts[1][c[1] as usize] += 1;
        counts[2][c[2] as usize] += 1;
        counts[3][c[3] as usize] += 1;
    }
    for &b in remainder {
        counts[0][b as usize] += 1;
    }
}

fn histogram_to_chars(counts: &[[usize; 256]; 4]) -> HashMap<char, usize> {
    let mut map = HashMap::new();
    for b in 0..=255u8 {
        let total: usize = counts.iter().map(|lane| lane[b as usize]).sum();
        if total > 0 {
            let ch = (b as char).to_ascii_uppercase();
            *map.entry(ch).or_insert(0) += total;
        }
    }
    map
}

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
    use core::str;

    use super::*;
    use crate::helper::sequence::SeqParser;
    use crate::helper::types::{DataType, InputFmt};

    const DNA: DataType = DataType::Dna;

    fn get_matrix(id: &[&str], seq: &[&str]) -> SeqMatrix {
        let mut matrix = SeqMatrix::new();
        id.iter().zip(seq.iter()).for_each(|(i, s)| {
            matrix.insert(i.to_string(), s.as_bytes());
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
        let smat = site.index_site_without_ambiguity(&mat, &DNA);
        site.compute_parsimony_informative(&smat);
        assert_eq!(1, site.pars_inf);
    }

    #[test]
    fn count_variable_sites_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT", "ATTA", "ATGC", "ATGA"];
        let mat = get_matrix(&id, &seq);
        let mut site = Sites::default();
        let site_mat = site.index_site_without_ambiguity(&mat, &DNA);
        site.compute_parsimony_informative(&site_mat);
        assert_eq!(1, site.pars_inf);
    }

    #[test]
    fn count_parsimony_gap_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT---", "ATTA---", "ATGC---", "ATGA---"];
        let mat = get_matrix(&id, &seq);
        let mut site = Sites::default();
        let smat = site.index_site_without_ambiguity(&mat, &DNA);
        site.compute_parsimony_informative(&smat);
        assert_eq!(1, site.pars_inf);
        assert_eq!(3, site.variable);
    }

    #[test]
    fn get_site_stats_test() {
        let path = Path::new("../../tests/files/concat.fasta");
        let input_format = InputFmt::Fasta;
        let aln = SeqParser::new(path, &DNA);
        let (matrix, _) = aln.get_alignment(&input_format);
        let mut site = Sites::default();
        let smat = site.index_site_without_ambiguity(&matrix, &DNA);
        site.compute_parsimony_informative(&smat);
        assert_eq!(18, site.conserved);
        assert_eq!(8, site.variable);
        assert_eq!(2, site.pars_inf);
    }

    #[test]
    fn test_site_missing_data() {
        let mut index_site = HashMap::new();
        // Create dummy hash map data with sites
        index_site.insert(0, vec![b'A', b'T', b'?', b'-']);
        let site = Sites::default();
        let data = site.get_missing_data_per_site(&index_site, 4);
        assert_eq!(1, data.len());
        let site_length = index_site.get(&0).unwrap().len();
        let missing_data = 2;

        assert_eq!(missing_data as f64 / site_length as f64, data[0].1,);
    }

    #[test]
    fn test_site_indexing() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AA--", "AT--", "ATGC", "ATGA"];
        let mat = get_matrix(&id, &seq);
        let site = Sites::default();
        let site_mat = site.index_site_with_missing_data(&mat);
        assert_eq!(4, site_mat.len());
        assert_eq!(4, site_mat.get(&0).unwrap().len());
        let missing = site.get_missing_data_per_site(&site_mat, 4);
        let site_dna = "--GG";
        let matrix_seq = site_mat.get(&2).unwrap();
        let matrix_dna = str::from_utf8(matrix_seq.as_slice()).unwrap();
        assert_eq!(site_dna, matrix_dna);
        assert_eq!(4, missing.len());
        let threshold = 0.4;
        assert!(0.5 > threshold);
        let sites = site.get_site_without_missing_data(&mat, threshold);
        assert_eq!(2, sites.len());
    }

    #[test]
    fn filter_min_tax_test() {
        let ntax = vec![10, 8, 20, 30, 60];
        let comp = Completeness::new(&60, 2);
        assert_eq!(2, comp.count_mat_completeness(&ntax, 0.5))
    }

    #[test]
    fn char_count_test() {
        let path = Path::new("../../tests/files/concat.fasta");
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

    #[test]
    fn bitmask_matches_legacy_sitemap_dna_and_aa() {
        let ids = ["t1", "t2", "t3", "t4"];
        // Site 0: conserved mixed-case a/A
        // Site 1: variable non-PIS (A, A, C, G)
        // Site 2: PIS (A, A, T, T)
        // Site 3: all gaps/ambiguous (skipped)
        // Site 4: PIS with trailing gap (C, C, g, -)
        let dna_seqs = ["aAC-C", "AAC?C", "ACtNg", "AGt-g"];
        let dna_mat = get_matrix(&ids, &dna_seqs);

        let mut fast_dna = Sites::default();
        fast_dna.get_stats(&dna_mat, &DataType::Dna);

        let mut legacy_dna = Sites::default();
        let smat = legacy_dna.index_site_without_ambiguity(&dna_mat, &DataType::Dna);
        legacy_dna.compute_parsimony_informative(&smat);
        legacy_dna.count_sites();

        assert_eq!(fast_dna.conserved, legacy_dna.conserved);
        assert_eq!(fast_dna.variable, legacy_dna.variable);
        assert_eq!(fast_dna.pars_inf, legacy_dna.pars_inf);
        assert_eq!(fast_dna.counts, legacy_dna.counts);
        assert_eq!(fast_dna.conserved, 1);
        assert_eq!(fast_dna.variable, 3);
        assert_eq!(fast_dna.pars_inf, 2);

        let pis_sites =
            Sites::default().get_site_with_pars_informative(&dna_mat, &DataType::Dna, 2);
        assert_eq!(pis_sites, vec![(2, 2), (4, 2)]);

        let aa_seqs = ["mARK-", "MaSK?", "MASRx", "MASK*"];
        let aa_mat = get_matrix(&ids, &aa_seqs);
        let mut fast_aa = Sites::default();
        fast_aa.get_stats(&aa_mat, &DataType::Aa);

        let mut legacy_aa = Sites::default();
        let aa_smat = legacy_aa.index_site_without_ambiguity(&aa_mat, &DataType::Aa);
        legacy_aa.compute_parsimony_informative(&aa_smat);
        legacy_aa.count_sites();

        assert_eq!(fast_aa.conserved, legacy_aa.conserved);
        assert_eq!(fast_aa.variable, legacy_aa.variable);
        assert_eq!(fast_aa.pars_inf, legacy_aa.pars_inf);
        assert_eq!(fast_aa.counts, legacy_aa.counts);
    }
}
