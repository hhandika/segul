//! Statistics for contigs
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use flate2::read::MultiGzDecoder;

use crate::{
    helper::types::{ContigFmt, infer_contig_fmt_auto},
    stats::common::{CommonStats, NStats},
};

/// Lower bounds of the contig length bins, in bp. The last bin is open.
pub const LENGTH_BIN_EDGES: [usize; 10] = [
    0, 200, 500, 1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000,
];

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
    /// Contigs per length bin, aligned with [`LENGTH_BIN_EDGES`].
    pub length_distribution: Vec<usize>,
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
            length_distribution: vec![0; LENGTH_BIN_EDGES.len()],
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
        const BUF_SIZE: usize = 64 * 1024;
        let file = File::open(path).expect("Failed opening a contig file");
        match file_fmt {
            ContigFmt::Fasta => {
                let mut buff = BufReader::with_capacity(BUF_SIZE, file);
                self.count(&mut buff)
            }
            ContigFmt::Gzip => {
                let mut buff = BufReader::with_capacity(BUF_SIZE, MultiGzDecoder::new(file));
                self.count(&mut buff)
            }
            ContigFmt::Auto => {
                let file_fmt = infer_contig_fmt_auto(path);
                self.parse_file(path, &file_fmt)
            }
        }
    }

    fn count<R: BufRead>(&mut self, buff: &mut R) -> Vec<usize> {
        let mut contigs = Vec::new();
        // Byte histogram for the whole file, split into four lanes
        // so consecutive identical bases do not update the same counter.
        let mut counts = [[0usize; 256]; 4];
        let mut line = Vec::new();
        let mut current_len = 0usize;
        let mut in_record = false;
        while matches!(buff.read_until(b'\n', &mut line), Ok(n) if n > 0) {
            let trimmed = line.trim_ascii_end();
            if trimmed.starts_with(b">") {
                if in_record {
                    contigs.push(current_len);
                }
                in_record = true;
                current_len = 0;
            } else if in_record {
                let seq = trimmed.trim_ascii_start();
                current_len += seq.len();
                count_bytes(seq, &mut counts);
            }
            line.clear();
        }
        if in_record {
            contigs.push(current_len);
        }
        self.tally_bases(&counts);
        contigs
    }

    fn tally_bases(&mut self, counts: &[[usize; 256]; 4]) {
        let total = |b: u8| -> usize { counts.iter().map(|lane| lane[b as usize]).sum() };
        let all: usize = counts.iter().flatten().sum();
        self.g_count += total(b'G') + total(b'g');
        self.c_count += total(b'C') + total(b'c');
        self.a_count += total(b'A') + total(b'a');
        self.t_count += total(b'T') + total(b't');
        self.n_count += total(b'N') + total(b'n');
        let known = b"GgCcAaTtNn".iter().map(|&b| total(b)).sum::<usize>();
        self.unknown += all - known;
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
        self.count_length_bins(contigs);
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
        self.gc_content = proportion(self.g_count + self.c_count, self.base_count);
    }

    fn calculate_at_content(&mut self) {
        self.at_content = proportion(self.a_count + self.t_count, self.base_count);
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

    fn count_length_bins(&mut self, contigs: &[usize]) {
        contigs
            .iter()
            .for_each(|&len| self.length_distribution[length_bin(len)] += 1);
    }
}

/// Index of the length bin holding a contig of `len` bp.
pub fn length_bin(len: usize) -> usize {
    LENGTH_BIN_EDGES.partition_point(|&edge| edge <= len) - 1
}

/// Share of `count` in `total`, 0 when `total` is 0.
pub(crate) fn proportion(count: usize, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        count as f64 / total as f64
    }
}

fn count_bytes(seq: &[u8], counts: &mut [[usize; 256]; 4]) {
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
        assert_eq!(parser.length_distribution[1], 2);
        assert_eq!(parser.length_distribution.iter().sum::<usize>(), 2);
    }

    #[test]
    fn test_length_bin() {
        assert_eq!(length_bin(0), 0);
        assert_eq!(length_bin(199), 0);
        assert_eq!(length_bin(200), 1);
        assert_eq!(length_bin(1_000), 3);
        assert_eq!(length_bin(5_000_000), LENGTH_BIN_EDGES.len() - 1);
    }

    #[test]
    fn test_empty_contigs() {
        let mut parser = ContigSummary::new();
        parser.summarize_contigs(&[]);
        assert_eq!(parser.gc_content, 0.0);
        assert_eq!(parser.at_content, 0.0);
    }

    #[test]
    fn test_count_bytes_mixed_case_and_unknown() {
        let mut parser = ContigSummary::new();
        let mut counts = [[0usize; 256]; 4];
        count_bytes(b"GgCcAaTtNnRY-", &mut counts);
        parser.tally_bases(&counts);
        assert_eq!(
            (2, 2, 2, 2, 2, 3),
            (
                parser.g_count,
                parser.c_count,
                parser.a_count,
                parser.t_count,
                parser.n_count,
                parser.unknown
            )
        );
    }

    #[test]
    fn test_streaming_wrapped_crlf_and_gzip_contigs() {
        use flate2::{Compression, write::GzEncoder};
        use std::io::Write;
        use tempfile::TempDir;

        let dir = TempDir::with_prefix("contig_stream").unwrap();
        let fasta_path = dir.path().join("wrapped.fa");
        let gz_path = dir.path().join("wrapped.fa.gz");
        let empty_path = dir.path().join("empty.fa");
        let content = b">c1 desc\r\nACGT\r\nTGCA\r\n>c2\r\nGG\r\n";
        std::fs::write(&fasta_path, content).unwrap();
        std::fs::write(&empty_path, b"").unwrap();
        let mut gz = GzEncoder::new(Vec::new(), Compression::default());
        gz.write_all(content).unwrap();
        std::fs::write(&gz_path, gz.finish().unwrap()).unwrap();

        let mut plain = ContigSummary::new();
        plain.summarize(&fasta_path, &ContigFmt::Fasta);
        let mut gzipped = ContigSummary::new();
        gzipped.summarize(&gz_path, &ContigFmt::Gzip);
        assert_eq!(plain.contig_count, 2);
        assert_eq!(plain.base_count, 10);
        assert_eq!(plain.stats.min, 2);
        assert_eq!(plain.stats.max, 8);
        assert_eq!(gzipped.contig_count, plain.contig_count);
        assert_eq!(gzipped.base_count, plain.base_count);
        assert_eq!(gzipped.g_count, plain.g_count);

        let mut empty = ContigSummary::new();
        empty.summarize(&empty_path, &ContigFmt::Fasta);
        assert_eq!(empty.contig_count, 0);
        assert_eq!(empty.base_count, 0);
    }
}
