//! Split a concatenated alignment
//! into multiple alignments based on partitions
use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{Result, bail};
use colored::Colorize;
use rayon::prelude::*;

use crate::core::{OutputPrint, PartitionPrint};
use crate::helper::files;
use crate::helper::sequence::SeqParser;
use crate::helper::types::{
    DataType, Header, InputFmt, OutputFmt, Partition, PartitionFmt, SeqMatrix,
};
use crate::helper::utils;
use crate::parser::partition::PartitionParser;
use crate::report::spec::{AlignSplitParams, TaskKind, canonical};
use crate::report::stats::TaskCounts;
use crate::report::{DirSnapshot, ReportContext, RunReport, RunReportBuilder, report_path};
use crate::writer::sequences::SeqWriter;

impl OutputPrint for AlignmentSplitting<'_> {}
impl PartitionPrint for AlignmentSplitting<'_> {}

pub struct AlignmentSplitting<'a> {
    input: &'a Path,
    datatype: &'a DataType,
    input_fmt: &'a InputFmt,
    output: &'a Path,
    output_fmt: &'a OutputFmt,
    report_context: ReportContext,
}

impl<'a> AlignmentSplitting<'a> {
    pub fn new(
        input: &'a Path,
        datatype: &'a DataType,
        input_fmt: &'a InputFmt,
        output: &'a Path,
        output_fmt: &'a OutputFmt,
    ) -> Self {
        Self {
            input,
            datatype,
            input_fmt,
            output,
            output_fmt,
            report_context: ReportContext::default(),
        }
    }

    /// Attach frontend provenance and SHA-256 checksum settings for `report.json`.
    pub fn with_report_context(mut self, context: ReportContext) -> Self {
        self.report_context = context;
        self
    }

    /// Split the alignment by partitions, panicking on error.
    pub fn split(
        &self,
        part_path: &Path,
        partition_fmt: &PartitionFmt,
        prefix: &Option<String>,
        is_uncheck: bool,
    ) {
        self.try_split(part_path, partition_fmt, prefix, is_uncheck)
            .unwrap_or_else(|err| panic!("{err}"));
    }

    /// Split the alignment by partitions and write a unified `report.json`.
    pub fn try_split_with_report(
        &self,
        part_path: &Path,
        partition_fmt: &PartitionFmt,
        prefix: &Option<String>,
        is_uncheck: bool,
    ) -> Result<RunReport> {
        let params = AlignSplitParams {
            input_format: canonical::input_fmt(self.input_fmt).to_string(),
            datatype: canonical::datatype(self.datatype).to_string(),
            output_format: canonical::output_fmt(self.output_fmt).to_string(),
            partition_format: canonical::partition_fmt(partition_fmt).to_string(),
            skip_checking: is_uncheck,
            prefix: prefix.clone(),
        };
        let mut aux = Vec::new();
        if part_path != self.input {
            aux.push(("partition", part_path.to_path_buf()));
        }
        let snapshot = DirSnapshot::capture(self.output);
        let builder = RunReportBuilder::new(TaskKind::AlignSplit, &self.report_context)
            .primary_inputs(&[self.input.to_path_buf()])
            .auxiliary_inputs(aux)
            .parameters(&params)?;
        self.try_split(part_path, partition_fmt, prefix, is_uncheck)?;
        let written = snapshot.written_files(self.output);
        let counts = TaskCounts {
            files_read: 1,
            files_written: written.len(),
            items_processed: Some(written.len()),
            items_skipped: None,
        };
        let report_file = report_path(self.output, prefix.as_deref());
        builder
            .outputs(self.output, written.into_iter().map(|p| ("alignment", p)))
            .statistics(&counts)?
            .write_to(&report_file)
    }

    /// Split the alignment by partitions, returning an error on invalid partitions or I/O failure.
    pub fn try_split(
        &self,
        part_path: &Path,
        partition_fmt: &PartitionFmt,
        prefix: &Option<String>,
        is_uncheck: bool,
    ) -> Result<()> {
        let partitions = PartitionParser::new(part_path, partition_fmt, is_uncheck).try_parse()?;
        self.print_partition_info(part_path, &partitions.len());
        let spin = utils::set_spinner();
        spin.set_message("Parsing input sequence file...");
        let aln_matrix = self.parse_sequence();
        let aln_len = aln_matrix.values().next().map_or(0, <[u8]>::len);
        self.validate_split_targets(&partitions, aln_len, prefix)?;

        spin.set_message("Splitting alignment...");
        let file_counts = AtomicUsize::new(0);
        partitions.par_iter().try_for_each(|part| -> Result<()> {
            let runs = part.contiguous_site_runs();
            let nchar: usize = runs.iter().map(|(s, e)| e - s).sum();
            let kept_rows: Vec<usize> = (0..aln_matrix.len())
                .filter(|&row| {
                    let (_, seq) = aln_matrix.get_index(row).expect("Row index out of bounds");
                    !runs
                        .iter()
                        .all(|&(start, end)| self.is_empty_seq(&seq[start..end]))
                })
                .collect();
            if kept_rows.is_empty() {
                log::warn!(
                    "Skipping partition '{}': all sequences are gaps or missing.",
                    part.gene
                );
                return Ok(());
            }
            let mut header = Header::new();
            if let DataType::Aa = self.datatype {
                header.datatype = String::from("protein");
            }
            header.nchar = nchar;
            header.ntax = kept_rows.len();
            header.aligned = true;
            let filename = self.parse_filename(&part.gene, prefix);
            let output_path = files::create_output_fname(self.output, &filename, self.output_fmt);
            let mut out = SeqWriter::new(&output_path, &aln_matrix, &header);
            out.write_runs(self.output_fmt, &kept_rows, &runs)?;
            file_counts.fetch_add(1, Ordering::Relaxed);
            Ok(())
        })?;

        spin.finish_with_message("Finished splitting alignment!\n");
        self.print_output_info(file_counts.load(Ordering::Relaxed));
        Ok(())
    }

    fn validate_split_targets(
        &self,
        partitions: &[Partition],
        aln_len: usize,
        prefix: &Option<String>,
    ) -> Result<()> {
        let mut out_paths = HashSet::with_capacity(partitions.len());
        for part in partitions {
            let max_site = part.max_site();
            if max_site > aln_len {
                bail!(
                    "Partition '{}' references site {}, which exceeds alignment length ({}).",
                    part.gene,
                    max_site,
                    aln_len
                );
            }
            let fname = self.parse_filename(&part.gene, prefix);
            let out_path = files::create_output_fname(self.output, &fname, self.output_fmt);
            if !out_paths.insert(out_path.clone()) {
                bail!(
                    "Multiple partitions resolve to the same output file '{}'.",
                    out_path.display()
                );
            }
        }
        Ok(())
    }

    // Generate a filename for each locus based on the locus name
    // We get rid of any characters that are not alphanumeric, underscore, or dash
    fn parse_filename(&self, gene_name: &str, prefix: &Option<String>) -> PathBuf {
        let filename = match prefix {
            Some(prefix) => format!("{}_{}", prefix, gene_name),
            None => String::from(gene_name),
        };

        PathBuf::from(filename)
    }

    #[cfg_attr(not(test), allow(dead_code))]
    fn generate_matrix_from_runs(
        &self,
        matrix: &SeqMatrix,
        runs: &[(usize, usize)],
        total_len: usize,
    ) -> SeqMatrix {
        let mut seq_matrix = SeqMatrix::with_capacity(matrix.len(), matrix.len() * total_len);
        for (taxon, seq) in matrix.iter() {
            let part_seq = self.slice_runs(seq, runs, total_len);
            if !self.is_empty_seq(&part_seq) {
                seq_matrix.insert(taxon.clone(), &part_seq);
            }
        }
        seq_matrix
    }

    fn is_empty_seq(&self, seq: &[u8]) -> bool {
        seq.iter().all(|&b| b == b'-' || b == b'?')
    }

    #[cfg_attr(not(test), allow(dead_code))]
    fn slice_runs(&self, seq: &[u8], runs: &[(usize, usize)], total_len: usize) -> Vec<u8> {
        let mut out = Vec::with_capacity(total_len);
        for &(start, end) in runs {
            let slice = seq
                .get(start..end)
                .expect("Errors in slicing the alignment");
            out.extend_from_slice(slice);
        }
        out
    }

    fn parse_sequence(&self) -> SeqMatrix {
        let aln = SeqParser::new(self.input, self.datatype);
        let (matrix, _) = aln.get_alignment(self.input_fmt);
        matrix
    }

    fn print_output_info(&self, file_counts: usize) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&file_counts));
        log::info!("{:18}: {}", "Output dir", self.output.display());
        self.print_output_fmt(self.output_fmt);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::helper::types::SiteRange;
    use tempfile::TempDir;

    macro_rules! input_split {
        ($var:ident, $input:expr_2021) => {
            let $var = AlignmentSplitting::new(
                &Path::new($input),
                &DataType::Dna,
                &InputFmt::Fasta,
                &Path::new("../../tests/files/"),
                &OutputFmt::Fasta,
            );
        };
    }

    #[test]
    fn test_parse_filename() {
        input_split!(split, "../../tests/files/test.fasta");
        assert_eq!(
            split.parse_filename("test", &Some("test".to_string())),
            PathBuf::from("test_test")
        );
    }

    #[test]
    fn test_sequence_slicing_runs() {
        input_split!(split, "../../tests/files/test_aln.fasta");
        let seq = b"AAAAAAAAGGGGTTTTCCCC";

        let slice = split.slice_runs(seq, &[(0, 10)], 10);
        let slice_2 = split.slice_runs(seq, &[(10, 15)], 5);
        assert_eq!(slice, b"AAAAAAAAGG");
        assert_eq!(slice_2, b"GGTTT");
    }

    #[test]
    fn test_generate_matrix_merged_and_out_of_order_column_order() {
        input_split!(split, "../../tests/files/partition/concat_part.fas");
        let matrix = split.parse_sequence();
        // Sites 11-15 and 1-5 given out of order -> sorted ascending [0..5, 10..15]
        let part = Partition::from_ranges(
            "merged",
            vec![SiteRange::new(11, 15), SiteRange::new(1, 5)],
            None,
        );
        let runs = part.contiguous_site_runs();
        assert_eq!(runs, vec![(0, 5), (10, 15)]);
        let new_matrix = split.generate_matrix_from_runs(&matrix, &runs, 10);
        assert_eq!(new_matrix.len(), 4);
        assert_eq!(new_matrix.get("ABCD").unwrap(), b"aaaaattttt");
        assert_eq!(new_matrix.get("ABCG").unwrap(), b"aaaaaccccc");
    }

    #[test]
    fn test_generate_matrix_strided_codon_column_order() {
        input_split!(split, "../../tests/files/partition/concat_part.fas");
        let matrix = split.parse_sequence();
        // 1-based sites: 1,4,7,10,13 (from 1-15\3) and 2,5,8,11,14 (from 2-15\3)
        // Sorted 1-based: 1,2, 4,5, 7,8, 10,11, 13,14 -> runs: (0,2),(3,5),(6,8),(9,11),(12,14)
        let part = Partition::from_ranges(
            "pos12",
            vec![
                SiteRange::with_step(1, 15, 3),
                SiteRange::with_step(2, 15, 3),
            ],
            None,
        );
        let runs = part.contiguous_site_runs();
        assert_eq!(runs, vec![(0, 2), (3, 5), (6, 8), (9, 11), (12, 14)]);
        let new_matrix = split.generate_matrix_from_runs(&matrix, &runs, 10);
        // ABCG is aaaaagggggccccc -> positions 1,2='aa', 4,5='aa', 7,8='gg', 10,11='gc', 13,14='cc'
        assert_eq!(new_matrix.get("ABCG").unwrap(), b"aaaagggccc");
    }

    #[test]
    fn test_split_out_of_bounds_returns_error() {
        let tmp = TempDir::with_prefix("segul_split_oob").unwrap();
        let part_file = tmp.path().join("oob.txt");
        std::fs::write(&part_file, "DNA, p1 = 1-20\n").unwrap();
        let split = AlignmentSplitting::new(
            Path::new("../../tests/files/partition/concat_part.fas"),
            &DataType::Dna,
            &InputFmt::Fasta,
            tmp.path(),
            &OutputFmt::Fasta,
        );
        let err = split
            .try_split(&part_file, &PartitionFmt::Raxml, &None, true)
            .unwrap_err();
        assert!(err.to_string().contains("exceeds alignment length"));
    }
}
