//! Convert sequence format from one to another.
use std::path::{Path, PathBuf};

use colored::Colorize;
use rayon::prelude::*;

use crate::core::OutputPrint;
use crate::helper::sequence::SeqParser;
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix};
use crate::helper::{files, utils};
use crate::writer::sequences::SeqWriter;

impl OutputPrint for AlignmentConversion<'_> {}

/// Convert sequence format
pub struct AlignmentConversion<'a> {
    input_fmt: &'a InputFmt,
    output_fmt: &'a OutputFmt,
    datatype: &'a DataType,
    sort: bool,
}

impl<'a> AlignmentConversion<'a> {
    pub fn new(
        input_fmt: &'a InputFmt,
        output_fmt: &'a OutputFmt,
        datatype: &'a DataType,
        sort: bool,
    ) -> Self {
        Self {
            input_fmt,
            output_fmt,
            datatype,
            sort,
        }
    }

    pub fn convert(&self, files: &[PathBuf], output: &Path) {
        let spin = utils::set_spinner();
        spin.set_message("Converting sequence format...");
        files.par_iter().for_each(|file| {
            let output_fname = files::create_output_fname(output, file, self.output_fmt);
            self.convert_any(file, &output_fname);
        });
        spin.finish_with_message("Finished converting sequence format!\n");
        self.print_output_info(output);
    }

    fn convert_any(&self, input: &Path, output: &Path) {
        if self.sort {
            self.convert_sorted(input, output);
        } else {
            self.convert_unsorted(input, output);
        }
    }

    fn convert_unsorted(&self, input: &Path, output: &Path) {
        let (matrix, header) = self.get_sequence(input);
        self.write_results(&matrix, header, output);
    }

    fn convert_sorted(&self, input: &Path, output: &Path) {
        let (mut matrix, header) = self.get_sequence(input);
        matrix.sort_keys();
        self.write_results(&matrix, header, output)
    }

    fn get_sequence(&self, input: &Path) -> (SeqMatrix, Header) {
        let seq = SeqParser::new(input, self.datatype);
        seq.parse(self.input_fmt)
    }

    fn write_results(&self, matrix: &SeqMatrix, header: Header, output: &Path) {
        let mut convert = SeqWriter::new(output, matrix, &header);
        convert
            .write_sequence(self.output_fmt)
            .expect("Failed writing output files");
    }

    fn print_output_info(&self, output: &Path) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Output dir", output.display());
        self.print_output_fmt(self.output_fmt);
    }
}
