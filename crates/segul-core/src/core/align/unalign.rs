//! Convert input alignments to unaligned sequence files.
//!
//!
//! # Example
//!
//! ```
//! use std::path::{Path, PathBuf};
//! use tempfile::TempDir;
//! use segul_core::helper::types::{DataType, InputFmt, OutputFmt, PartitionFmt};
//! use segul_core::core::align::unalign::UnalignAlignment;
//! use segul_core::helper::finder::SeqFileFinder;
//!
//! let input_fmt = InputFmt::Nexus;
//! let datatype = DataType::Dna;
//! let input_dir = Path::new("../../tests/files/alignments");
//! // Find matching alignment files in the input directory
//! let files = SeqFileFinder::new(Path::new(input_dir)).find(&input_fmt);
//! // Replace the temp directory with your own directory.
//! let output = TempDir::with_prefix("temp").unwrap();
//! let output_dir = output.path();
//! // Unalign feature only supports fasta or fasta-int output format
//! let output_fmt = OutputFmt::Fasta;
//! let handle = UnalignAlignment::new(
//!     &files,
//!     &input_fmt,
//!     &datatype,
//!     &output_dir,
//!     &output_fmt
//!     );
//! handle.unalign();

use colored::Colorize;
use rayon::prelude::*;
use std::path::{Path, PathBuf};

use crate::{
    core::OutputPrint,
    helper::{
        files,
        sequence::SeqParser,
        types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix},
        utils,
    },
    report::{
        DirSnapshot, ReportContext, RunReport, RunReportBuilder, report_path,
        spec::{AlignUnalignParams, TaskKind, canonical},
        stats::TaskCounts,
    },
    writer::sequences::SeqWriter,
};

/// Unalign alignment files
pub struct UnalignAlignment<'a> {
    /// Input alignment files
    input_files: &'a [PathBuf],
    /// Input format of alignment files
    input_fmt: &'a InputFmt,
    /// Data type of sequences
    datatype: &'a DataType,
    /// Output directory
    output_dir: &'a Path,
    /// Output format of sequence files. Always fasta
    output_fmt: &'a OutputFmt,
    report_context: ReportContext,
}

impl OutputPrint for UnalignAlignment<'_> {}

impl Default for UnalignAlignment<'_> {
    fn default() -> Self {
        Self {
            input_files: &[],
            input_fmt: &InputFmt::Fasta,
            datatype: &DataType::Dna,
            output_dir: Path::new(""),
            output_fmt: &OutputFmt::Fasta,
            report_context: ReportContext::default(),
        }
    }
}

impl<'a> UnalignAlignment<'a> {
    /// Create a new UnalignAlignment instance
    /// The output format is always fasta
    pub fn new(
        input_files: &'a [PathBuf],
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        output_dir: &'a Path,
        output_fmt: &'a OutputFmt,
    ) -> Self {
        Self {
            input_files,
            input_fmt,
            datatype,
            output_dir,
            output_fmt,
            report_context: ReportContext::default(),
        }
    }

    /// Attach frontend provenance and SHA-256 checksum settings for `report.json`.
    pub fn with_report_context(mut self, context: ReportContext) -> Self {
        self.report_context = context;
        self
    }

    /// Convert aligned sequences to unaligned sequences
    /// by removing gaps from each sequence
    pub fn unalign(&self) {
        self.unalign_with_report()
            .expect("Failed unaligning alignments");
    }

    /// Unalign sequences and write a unified `report.json`.
    pub fn unalign_with_report(&self) -> anyhow::Result<RunReport> {
        if self.output_fmt != &OutputFmt::Fasta && self.output_fmt != &OutputFmt::FastaInt {
            anyhow::bail!(
                "Unalign feature only supports fasta or fasta-int output format. Output format provided: {}",
                self.output_fmt
            );
        }
        let mut sorted = self.input_files.to_vec();
        alphanumeric_sort::sort_path_slice(&mut sorted);
        let params = AlignUnalignParams {
            input_format: canonical::input_fmt(self.input_fmt).to_string(),
            datatype: canonical::datatype(self.datatype).to_string(),
            output_format: canonical::output_fmt(self.output_fmt).to_string(),
        };
        let snapshot = DirSnapshot::capture(self.output_dir);
        let builder = RunReportBuilder::new(TaskKind::AlignUnalign, &self.report_context)
            .primary_inputs(&sorted)
            .parameters(&params)?;
        let spin = utils::set_spinner();
        spin.set_message("Converting un-aligned sequence files...");
        sorted.par_iter().for_each(|file| {
            let (matrix, header) = self.get_unalign(file);
            self.write_results(file, &matrix, &header);
        });
        spin.finish_with_message("Finished un-aligning alignments!\n");
        self.print_output_info();
        let written = snapshot.written_files(self.output_dir);
        let counts = TaskCounts {
            files_read: sorted.len(),
            files_written: written.len(),
            items_processed: Some(written.len()),
            items_skipped: None,
        };
        let report_file = report_path(self.output_dir, None);
        builder
            .outputs(
                self.output_dir,
                written.into_iter().map(|p| ("sequences", p)),
            )
            .statistics(&counts)?
            .write_to(&report_file)
    }

    fn get_unalign(&self, input: &Path) -> (SeqMatrix, Header) {
        let (mut matrix, header) =
            SeqParser::new(input, self.datatype).get_alignment(self.input_fmt);
        matrix.remove_bytes(|b| b == b'?' || b == b'-');

        (matrix, header)
    }

    #[cfg(test)]
    fn remove_gaps(&self, seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .copied()
            .filter(|&b| b != b'?' && b != b'-')
            .collect()
    }

    fn write_results(&self, input: &Path, matrix: &SeqMatrix, header: &Header) {
        let output_path = files::create_output_fname(self.output_dir, input, self.output_fmt);
        let mut writer = SeqWriter::new(&output_path, matrix, header);
        writer
            .write_sequence(self.output_fmt)
            .expect("Failed to write unaligned sequences");
    }

    fn print_output_info(&self) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Directory", self.output_dir.display());
        self.print_output_fmt(self.output_fmt);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unalign() {
        let seq = b"ATG-?C";
        let unalign = UnalignAlignment::default();
        let res = unalign.remove_gaps(seq);
        assert_eq!(res, b"ATGC");
    }
}
