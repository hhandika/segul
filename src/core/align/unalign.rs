//! Convert an alignment to unaligned sequences.
//!
//! Support multiple alignments.

use std::path::{Path, PathBuf};

use rayon::prelude::*;

use crate::{
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
    }

    fn get_unalign(&self, input: &Path) -> (SeqMatrix, Header) {
        let (mut matrix, header) =
            SeqParser::new(input, self.datatype).get_alignment(self.input_fmt);
        self.remove_gaps(&mut matrix);

        (matrix, header)
    }

    // Iterate over map and replace '?' with '-' of each values
    fn remove_gaps(&self, matrix: &mut SeqMatrix) {
        matrix.values_mut().for_each(|seq| {
            *seq = seq.replace(&['?', '-'], "-");
        });
    }

    fn write_results(&self, input: &Path, matrix: &SeqMatrix, header: &Header) {
        let output_path = files::create_output_fname(self.output_dir, input, self.output_fmt);
        let mut writer = SeqWriter::new(&output_path, matrix, header);
        writer
            .write_sequence(&self.output_fmt)
            .expect("Failed to write unaligned sequences");
    }
}
