//! Transpose all sequences from locus alignments to create individual-level reference sequences.
//!
//! Each output file contains all loci for a single sample.
//! Sequences are labeled with the locus (input file stem) they come from
//! and gaps (`-`) and missing data (`?`) are removed.
use std::collections::HashMap;
use std::path::{Path, PathBuf};

use colored::Colorize;
use rayon::prelude::*;

use crate::core::OutputPrint;
use crate::helper::files;
use crate::helper::sequence::SeqParser;
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix};
use crate::helper::utils;
use crate::writer::sequences::SeqWriter;

impl OutputPrint for SequenceTranspose<'_> {}

/// Transpose locus alignments into per-sample sequence files.
pub struct SequenceTranspose<'a> {
    /// The input format of the sequence alignment files.
    input_fmt: &'a InputFmt,
    /// The data type of the sequence alignment files.
    datatype: &'a DataType,
    /// The output directory. Each sample is written to
    /// `output_dir/<sample_id>/<sample_id>.<ext>`.
    output_dir: &'a Path,
    /// The output format of the transposed sequences.
    /// Use FASTA (`OutputFmt::Fasta` or `OutputFmt::FastaInt`).
    /// Transposed sequences are unaligned, so NEXUS and PHYLIP
    /// output would have an invalid character count.
    output_fmt: &'a OutputFmt,
}

impl<'a> SequenceTranspose<'a> {
    pub fn new(
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        output_dir: &'a Path,
        output_fmt: &'a OutputFmt,
    ) -> Self {
        Self {
            input_fmt,
            datatype,
            output_dir,
            output_fmt,
        }
    }

    /// Transpose sequences from locus alignments into per-sample files.
    /// Example:
    /// ```rust
    /// use std::path::PathBuf;
    /// use segul_core::core::sequence::transpose::SequenceTranspose;
    /// use segul_core::helper::types::{DataType, InputFmt, OutputFmt};
    /// use tempdir::TempDir;
    ///
    /// let files = vec![
    ///     PathBuf::from("tests/files/alignments/gene_1.nex"),
    ///     PathBuf::from("tests/files/alignments/gene_2.nex"),
    /// ];
    /// let output = TempDir::new("temp").unwrap();
    /// let handle = SequenceTranspose::new(
    ///     &InputFmt::Nexus,
    ///     &DataType::Dna,
    ///     output.path(),
    ///     &OutputFmt::Fasta,
    /// );
    /// handle.transpose(&files);
    /// assert!(output.path().join("ABCD").join("ABCD.fas").exists());
    /// ```
    pub fn transpose(&self, files: &[PathBuf]) {
        let spin = utils::set_spinner();
        spin.set_message("Transposing sequences...");
        let sequence_map = self.collect_sequences(files);
        self.write_transposed_sequences(&sequence_map);
        spin.finish_with_message("Finished transposing sequences!\n");
        self.print_output_info(sequence_map.len());
    }

    fn collect_sequences(&self, files: &[PathBuf]) -> HashMap<String, SeqMatrix> {
        let mut sequence_map: HashMap<String, SeqMatrix> = HashMap::new();
        files.iter().for_each(|file| {
            // We ignore the header because we only need the sequence
            let (matrix, _) = SeqParser::new(file, self.datatype).parse(self.input_fmt);
            let locus_id = self.get_locus_id(file);
            matrix.into_iter().for_each(|(sample_id, seq)| {
                // Add sample id as key if not exists
                // Or insert sequence (locus_id, sequence) to existing sample id
                sequence_map
                    .entry(sample_id)
                    .or_default()
                    .insert(locus_id.clone(), self.remove_gaps(&seq));
            });
        });
        sequence_map
    }

    fn remove_gaps(&self, seq: &str) -> String {
        seq.replace(['?', '-'], "")
    }

    fn write_transposed_sequences(&self, sequence_map: &HashMap<String, SeqMatrix>) {
        sequence_map.par_iter().for_each(|(sample_id, matrix)| {
            let output_path = self.get_output_path(sample_id);
            let mut header = Header::new();
            header.from_seq_matrix(matrix, false);
            let mut writer = SeqWriter::new(&output_path, matrix, &header);
            writer.write_sequence(self.output_fmt).unwrap_or_else(|e| {
                log::error!(
                    "Failed to write transposed sequence for sample {}: {}",
                    sample_id,
                    e
                );
            });
        });
    }

    fn get_output_path(&self, sample_id: &str) -> PathBuf {
        let output_path = self.output_dir.join(sample_id).join(sample_id);
        files::create_output_fname_from_path(&output_path, self.output_fmt)
    }

    fn get_locus_id(&self, file: &Path) -> String {
        file.file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string()
    }

    fn print_output_info(&self, sample_counts: usize) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Sample counts", utils::fmt_num(&sample_counts));
        log::info!("{:18}: {}", "Output dir", self.output_dir.display());
        self.print_output_fmt(self.output_fmt);
    }
}

#[cfg(test)]
mod tests {
    use std::fs;

    use tempdir::TempDir;

    use super::*;

    #[test]
    fn test_remove_gaps() {
        let handle = SequenceTranspose::new(
            &InputFmt::Fasta,
            &DataType::Dna,
            Path::new("."),
            &OutputFmt::Fasta,
        );
        assert_eq!(handle.remove_gaps("AT-G?C--"), "ATGC");
    }

    #[test]
    fn test_collect_sequences() {
        let files = vec![
            PathBuf::from("tests/files/alignments/gene_1.nex"),
            PathBuf::from("tests/files/alignments/gene_2.nex"),
        ];
        let handle = SequenceTranspose::new(
            &InputFmt::Nexus,
            &DataType::Dna,
            Path::new("."),
            &OutputFmt::Fasta,
        );
        let map = handle.collect_sequences(&files);
        let sample = map.get("ABCD").expect("Missing sample ABCD");
        assert!(sample.contains_key("gene_1"));
        assert!(sample.contains_key("gene_2"));
    }

    /// Write gappy FASTA locus alignments to `dir` and return their paths.
    /// `ABCF` is missing from `locus_2`.
    fn write_test_loci(dir: &Path) -> Vec<PathBuf> {
        let loci = [
            ("locus_1.fas", ">ABCD\nAT-GC?\n>ABCF\nATTGCA\n"),
            ("locus_2.fas", ">ABCD\n--GGCC\n"),
        ];
        loci.iter()
            .map(|(name, content)| {
                let path = dir.join(name);
                fs::write(&path, content).unwrap();
                path
            })
            .collect()
    }

    #[test]
    fn test_get_locus_id() {
        let handle = SequenceTranspose::new(
            &InputFmt::Nexus,
            &DataType::Dna,
            Path::new("."),
            &OutputFmt::Fasta,
        );
        let locus_id = handle.get_locus_id(Path::new("tests/files/alignments/gene_1.nex"));
        assert_eq!(locus_id, "gene_1");
    }

    #[test]
    fn test_get_output_path() {
        let output_dir = Path::new("output");
        let fasta = SequenceTranspose::new(
            &InputFmt::Nexus,
            &DataType::Dna,
            output_dir,
            &OutputFmt::Fasta,
        );
        assert_eq!(
            fasta.get_output_path("ABCD"),
            output_dir.join("ABCD").join("ABCD.fas")
        );
        let nexus = SequenceTranspose::new(
            &InputFmt::Nexus,
            &DataType::Dna,
            output_dir,
            &OutputFmt::Nexus,
        );
        assert_eq!(
            nexus.get_output_path("ABCD"),
            output_dir.join("ABCD").join("ABCD.nex")
        );
    }

    #[test]
    fn test_collect_sequences_all_samples() {
        let files = vec![
            PathBuf::from("tests/files/alignments/gene_1.nex"),
            PathBuf::from("tests/files/alignments/gene_2.nex"),
        ];
        let handle = SequenceTranspose::new(
            &InputFmt::Nexus,
            &DataType::Dna,
            Path::new("."),
            &OutputFmt::Fasta,
        );
        let map = handle.collect_sequences(&files);
        assert_eq!(map.len(), 2);
        // ABCF is only present in gene_2
        let sample = map.get("ABCF").expect("Missing sample ABCF");
        assert_eq!(sample.len(), 1);
        assert_eq!(sample.get("gene_2").unwrap(), "gatataat");
    }

    #[test]
    fn test_collect_sequences_removes_gaps() {
        let input = TempDir::new("temp").unwrap();
        let files = write_test_loci(input.path());
        let handle = SequenceTranspose::new(
            &InputFmt::Fasta,
            &DataType::Dna,
            Path::new("."),
            &OutputFmt::Fasta,
        );
        let map = handle.collect_sequences(&files);
        let abcd = map.get("ABCD").expect("Missing sample ABCD");
        assert_eq!(abcd.get("locus_1").unwrap(), "ATGC");
        assert_eq!(abcd.get("locus_2").unwrap(), "GGCC");
        let abcf = map.get("ABCF").expect("Missing sample ABCF");
        assert_eq!(abcf.len(), 1);
        assert_eq!(abcf.get("locus_1").unwrap(), "ATTGCA");
    }

    #[test]
    fn test_collect_sequences_keeps_locus_order() {
        let input = TempDir::new("temp").unwrap();
        let mut files = write_test_loci(input.path());
        files.reverse();
        let handle = SequenceTranspose::new(
            &InputFmt::Fasta,
            &DataType::Dna,
            Path::new("."),
            &OutputFmt::Fasta,
        );
        let map = handle.collect_sequences(&files);
        let loci: Vec<&String> = map.get("ABCD").unwrap().keys().collect();
        assert_eq!(loci, vec!["locus_2", "locus_1"]);
    }

    #[test]
    fn test_transpose_fasta_output() {
        let input = TempDir::new("temp").unwrap();
        let output = TempDir::new("temp").unwrap();
        let files = write_test_loci(input.path());
        let handle = SequenceTranspose::new(
            &InputFmt::Fasta,
            &DataType::Dna,
            output.path(),
            &OutputFmt::Fasta,
        );
        handle.transpose(&files);

        let abcd = output.path().join("ABCD").join("ABCD.fas");
        let abcf = output.path().join("ABCF").join("ABCF.fas");
        assert!(abcd.is_file());
        assert!(abcf.is_file());

        let (matrix, _) = SeqParser::new(&abcd, &DataType::Dna).parse(&InputFmt::Fasta);
        assert_eq!(matrix.len(), 2);
        assert_eq!(matrix.get("locus_1").unwrap(), "ATGC");
        assert_eq!(matrix.get("locus_2").unwrap(), "GGCC");

        let (matrix, _) = SeqParser::new(&abcf, &DataType::Dna).parse(&InputFmt::Fasta);
        assert_eq!(matrix.len(), 1);
        assert_eq!(matrix.get("locus_1").unwrap(), "ATTGCA");
    }

    #[test]
    fn test_transpose_interleaved_fasta_output() {
        let output = TempDir::new("temp").unwrap();
        let files = vec![
            PathBuf::from("tests/files/alignments/gene_1.nex"),
            PathBuf::from("tests/files/alignments/gene_2.nex"),
        ];
        let handle = SequenceTranspose::new(
            &InputFmt::Nexus,
            &DataType::Dna,
            output.path(),
            &OutputFmt::FastaInt,
        );
        handle.transpose(&files);

        let abcd = output.path().join("ABCD").join("ABCD.fas");
        assert!(abcd.is_file());
        let (matrix, _) = SeqParser::new(&abcd, &DataType::Dna).parse(&InputFmt::Fasta);
        assert_eq!(matrix.get("gene_1").unwrap(), "gatata");
        assert_eq!(matrix.get("gene_2").unwrap(), "gatataaa");
        assert!(output.path().join("ABCF").join("ABCF.fas").is_file());
    }
}
