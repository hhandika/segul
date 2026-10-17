//! Statistics for contigs
use std::{io::BufRead, path::Path};

use crate::{
    helper::{
        files,
        types::{infer_contig_fmt_auto, ContigFmt},
    },
    stats::common::{CommonStats, NStats},
};

use crate::parser::fasta::FastaReader;

/// Summary statistics for contigs
pub struct ContigSummary {
    /// Path to the file
    pub file_path: String,
    /// Name of the contig file
    pub contig_name: String,
    /// Contig count
    pub contig_count: usize,
    /// Total number of bases
    pub base_count: usize,
    /// Total number of nucleotides
    pub nucleotide: usize,
    /// Total number of G nucleotides
    pub g_count: usize,
    /// Total number of C nucleotides
    pub c_count: usize,
    /// Total number of A nucleotides
    pub a_count: usize,
    /// Total number of T nucleotides
    pub t_count: usize,
    /// Total number of N nucleotides
    pub n_count: usize,
    /// Total number of unknown characters
    pub unknown: usize,
    /// GC content
    pub gc_content: f64,
    /// AT content
    pub at_content: f64,
    /// N50
    pub nstats: NStats,
    /// Common stat for contigs
    /// Mean, median, min, max, and stdev
    pub stats: CommonStats,
    /// Number of contigs > 750 bp
    pub contig750: usize,
    /// Number of contigs > 1000 bp
    pub contig1000: usize,
    /// Number of contigs > 1500 bp
    pub contig1500: usize,
}

impl Default for ContigSummary {
    fn default() -> Self {
        Self::new()
    }
}

impl ContigSummary {
    pub fn new() -> Self {
        Self {
            file_path: String::new(),
            contig_name: String::new(),
            contig_count: 0,
            base_count: 0,
            nucleotide: 0,
            g_count: 0,
            c_count: 0,
            a_count: 0,
            t_count: 0,
            n_count: 0,
            unknown: 0,
            gc_content: 0.0,
            at_content: 0.0,
            nstats: NStats::new(),
            stats: CommonStats::new(),
            contig750: 0,
            contig1000: 0,
            contig1500: 0,
        }
    }

    pub fn summarize(&mut self, path: &Path, file_fmt: &ContigFmt) {
        self.file_path = path.display().to_string();
        self.contig_name = path
            .file_stem()
            .expect("No file name")
            .to_str()
            .expect("File name not valid UTF-8")
            .to_string();
        let contigs = self.parse_file(path, file_fmt);
        self.summarize_contigs(&contigs);
    }

    fn parse_file(&mut self, path: &Path, file_fmt: &ContigFmt) -> Vec<usize> {
        match file_fmt {
            ContigFmt::Fasta => {
                let mut buff = files::open_file(path);
                self.count(&mut buff)
            }
            ContigFmt::Gzip => {
                let mut buff = files::decode_gzip(path);
                self.count(&mut buff)
            }
            ContigFmt::Auto => {
                let file_fmt = infer_contig_fmt_auto(path);
                self.parse_file(path, &file_fmt)
            }
        }
    }

    fn count<R: BufRead>(&mut self, buff: &mut R) -> Vec<usize> {
        let reader = FastaReader::new(buff);
        let mut contigs = Vec::new();
        reader.into_iter().for_each(|r| {
            contigs.push(r.seq.len());
            r.seq.bytes().for_each(|s| match s {
                b'G' | b'g' => self.g_count += 1,
                b'C' | b'c' => self.c_count += 1,
                b'A' | b'a' => self.a_count += 1,
                b'T' | b't' => self.t_count += 1,
                b'N' | b'n' => self.n_count += 1,
                _ => self.unknown += 1,
            })
        });
        contigs
    }

    fn summarize_contigs(&mut self, contigs: &[usize]) {
        self.common_stats(contigs);
        self.calculate_nstat(contigs);
        self.count_contigs(contigs);
        self.count_bases();
        self.count_nucleotide();
        self.calculate_gc_content();
        self.calculate_at_content();
        self.count_contig750(contigs);
        self.count_contig1000(contigs);
        self.count_contig1500(contigs);
    }

    fn count_contigs(&mut self, contigs: &[usize]) {
        self.contig_count = contigs.len();
    }

    fn common_stats(&mut self, contigs: &[usize]) {
        self.stats.calculate(contigs);
    }

    fn calculate_nstat(&mut self, contigs: &[usize]) {
        self.nstats.calculate(contigs, self.stats.sum);
    }

    fn count_bases(&mut self) {
        self.base_count = self.g_count + self.c_count + self.a_count + self.t_count + self.n_count;
    }

    fn count_nucleotide(&mut self) {
        self.nucleotide = self.base_count - self.n_count;
    }

    fn calculate_gc_content(&mut self) {
        self.gc_content = (self.g_count + self.c_count) as f64 / self.base_count as f64;
    }

    fn calculate_at_content(&mut self) {
        self.at_content = (self.a_count + self.t_count) as f64 / self.base_count as f64;
    }

    fn count_contig750(&mut self, contigs: &[usize]) {
        self.contig750 = contigs.iter().filter(|&c| *c > 750).count();
    }

    fn count_contig1000(&mut self, contigs: &[usize]) {
        self.contig1000 = contigs.iter().filter(|&c| *c > 1000).count();
    }

    fn count_contig1500(&mut self, contigs: &[usize]) {
        self.contig1500 = contigs.iter().filter(|&c| *c > 1500).count();
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_parse_fasta() {
        let mut parser = ContigSummary::new();
        let path = Path::new("tests/files/contigs/contigs1.fa");
        let file_fmt = ContigFmt::Fasta;
        parser.summarize(path, &file_fmt);
        assert_eq!(parser.contig_name, "contigs1");
        assert_eq!(parser.file_path, "tests/files/contigs/contigs1.fa");
        assert_eq!(parser.g_count, 53);
        assert_eq!(parser.c_count, 121);
        assert_eq!(parser.a_count, 187);
        assert_eq!(parser.t_count, 124);
        assert_eq!(parser.n_count, 0);
        assert_eq!(parser.unknown, 0);
        assert_eq!(parser.base_count, 485);
        assert_eq!(parser.nucleotide, 485);
        assert_eq!(parser.contig_count, 2);
        assert_eq!(parser.stats.sum, 485);
        assert_eq!(parser.stats.min, 227);
        assert_eq!(parser.stats.max, 258);
        assert_eq!(parser.stats.mean, 242.5);
        assert_eq!(parser.stats.median, 242.5);
        assert_eq!(parser.contig1000, 0);
        assert_eq!(parser.contig1500, 0);
        assert_eq!(parser.contig750, 0);
    }
}
