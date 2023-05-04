//! Statistics for contigs
use std::{io::BufRead, path::Path};

use crate::{
    helper::{
        files,
        types::{infer_contig_fmt_auto, ContigFmt},
    },
    stats::math::{self, NStats},
};

use super::fasta::FastaReader;

/// Summary statistics for contigs
pub struct ContigSummaryParser {
    /// Path to the file
    pub file_path: String,
    /// Name of the file
    pub file_name: String,
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
    pub n50: usize,
    /// N75
    pub n75: usize,
    /// N90
    pub n90: usize,
    /// Mean
    pub mean: f64,
    /// Median
    pub median: f64,
    /// Sum of all contig lengths
    pub total_len: usize,
    /// Maximum contig length
    pub max_len: usize,
    /// Minimum contig length
    pub min_len: usize,
    /// Number of contigs > 750 bp
    pub contig750_count: usize,
    /// Number of contigs > 1000 bp
    pub contig1000_count: usize,
    /// Number of contigs > 1500 bp
    pub contig1500_count: usize,
}

impl ContigSummaryParser {
    pub fn new() -> Self {
        Self {
            file_path: String::new(),
            file_name: String::new(),
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
            n50: 0,
            n75: 0,
            n90: 0,
            mean: 0.0,
            median: 0.0,
            total_len: 0,
            max_len: 0,
            min_len: 0,
            contig750_count: 0,
            contig1000_count: 0,
            contig1500_count: 0,
        }
    }

    pub fn parse(&mut self, path: &Path, file_fmt: &ContigFmt) {
        self.file_path = path.display().to_string();
        self.file_name = path
            .file_name()
            .expect("No file name")
            .to_str()
            .expect("File name not valid UTF-8")
            .to_string();
        let contigs = self.parse_file(path, file_fmt);
        self.summarize(&contigs);
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

    fn summarize(&mut self, contigs: &[usize]) {
        self.total_len = self.total_len(contigs);
        let mut nstats = NStats::new(self.total_len);
        nstats.count(contigs);
        self.n50 = nstats.n50;
        self.n75 = nstats.n75;
        self.n90 = nstats.n90;
        self.contig_count = contigs.len();
        self.base_count = self.base_count();
        self.nucleotide = self.nucleotide();
        self.gc_content = self.gc_content();
        self.at_content = self.at_content();
        self.mean = self.mean(contigs);
        self.median = self.median(contigs);
        self.max_len = self.max(contigs);
        self.min_len = self.min(contigs);
        self.contig750_count = self.contig750(contigs);
        self.contig1000_count = self.contig1000(contigs);
        self.contig1500_count = self.contig1500(contigs);
    }

    fn total_len(&mut self, contigs: &[usize]) -> usize {
        contigs.iter().sum::<usize>()
    }

    fn mean(&mut self, contigs: &[usize]) -> f64 {
        math::mean(contigs, self.total_len)
    }

    fn median(&mut self, contigs: &[usize]) -> f64 {
        math::median(contigs)
    }

    fn min(&mut self, contigs: &[usize]) -> usize {
        *contigs.iter().min().unwrap_or(&0)
    }

    fn max(&mut self, contigs: &[usize]) -> usize {
        *contigs.iter().max().unwrap_or(&0)
    }

    fn base_count(&mut self) -> usize {
        self.g_count + self.c_count + self.a_count + self.t_count + self.n_count
    }

    fn nucleotide(&mut self) -> usize {
        self.base_count - self.n_count
    }

    fn gc_content(&mut self) -> f64 {
        (self.g_count + self.c_count) as f64 / self.base_count as f64
    }

    fn at_content(&mut self) -> f64 {
        (self.a_count + self.t_count) as f64 / self.base_count as f64
    }

    fn contig750(&mut self, contigs: &[usize]) -> usize {
        contigs.iter().filter(|&c| *c > 750).count()
    }

    fn contig1000(&mut self, contigs: &[usize]) -> usize {
        contigs.iter().filter(|&c| *c > 1000).count()
    }

    fn contig1500(&mut self, contigs: &[usize]) -> usize {
        contigs.iter().filter(|&c| *c > 1500).count()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_parse_fasta() {
        let mut parser = ContigSummaryParser::new();
        let path = Path::new("tests/files/contigs/contigs1.fa");
        let file_fmt = ContigFmt::Fasta;
        parser.parse(path, &file_fmt);
        assert_eq!(parser.file_name, "contigs1.fa");
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
    }
}
