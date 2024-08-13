//! Convert input alignments to unaligned sequence files.

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
        }
    }

    /// Convert aligned sequences to unaligned sequences
    /// by removing gaps from each sequence
    ///
    /// # Example
    ///
    /// ```
    /// use std::path::{Path, PathBuf};
    /// use tempdir::TempDir;
    /// use segul::helper::types::{DataType, InputFmt, OutputFmt, PartitionFmt};
    /// use segul::core::align::unalign::UnalignAlignment;
    /// use segul::helper::finder::SeqFileFinder;
    ///
    /// let input_fmt = InputFmt::Nexus;
    /// let datatype = DataType::Dna;
    /// let input_dir = Path::new("tests/files/alignments");
    /// // Find matching alignment files in the input directory
    /// let files = SeqFileFinder::new(Path::new(input_dir)).find(&input_fmt);
    /// // Replace the temp directory with your own directory.
    /// let output = TempDir::new("temp").unwrap();
    /// let output_dir = output.path();
    /// // Unalign feature only supports fasta or fasta-int output format
    /// let output_fmt = OutputFmt::Fasta;
    /// let handle = UnalignAlignment::new(
    ///     &files,
    ///     &input_fmt,
    ///     &datatype,
    ///     &output_dir,
    ///     &output_fmt
    ///     );
    /// handle.unalign();
    pub fn unalign(&self) {
        if self.output_fmt != &OutputFmt::Fasta && self.output_fmt != &OutputFmt::FastaInt {
            let msg = format!(
                "Unalign feature only supports fasta or fasta-int output format.\n\
                Output format provided: {}",
                self.output_fmt
            );
            log::warn!("{}", msg);
            return;
        }
        let spin = utils::set_spinner();
        spin.set_message("Converting un-aligned sequence files...");
        self.input_files.par_iter().for_each(|file| {
            let (matrix, header) = self.get_unalign(file);
            self.write_results(file, &matrix, &header);
        });
        spin.finish_with_message("Finished un-aligning alignments!\n");
        self.print_output_info();
    }

    fn get_unalign(&self, input: &Path) -> (SeqMatrix, Header) {
        let (mut matrix, header) =
            SeqParser::new(input, self.datatype).get_alignment(self.input_fmt);
        matrix.values_mut().for_each(|seq| {
            *seq = self.remove_gaps(seq);
        });

        (matrix, header)
    }

    // Iterate over map and replace '?' with '-' of each values
    fn remove_gaps(&self, seq: &str) -> String {
        seq.replace(['?', '-'], "")
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
        let seq = "ATG-?C";
        let unalign = UnalignAlignment::default();
        let res = unalign.remove_gaps(seq);
        assert!(res == String::from("ATGC"));
    }
}
