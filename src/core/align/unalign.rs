//! Convert an alignment to unaligned sequences.
//!
//! Support multiple alignments.

use std::path::{Path, PathBuf};

use colored::Colorize;
use rayon::prelude::*;

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

pub struct UnalignAlignment<'a> {
    input_files: &'a [PathBuf],
    input_fmt: &'a InputFmt,
    datatype: &'a DataType,
    output_dir: &'a Path,
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

    pub fn unalign(&self) {
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
