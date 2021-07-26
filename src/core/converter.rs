use std::path::Path;

use indexmap::IndexMap;

use crate::writer::seqwriter::SeqWriter;

use crate::helper::sequence::Sequence;
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, PartitionFmt};

pub struct Converter<'a> {
    input: &'a Path,
    output: &'a Path,
    output_fmt: &'a OutputFmt,
    is_dir: bool,
    datatype: &'a DataType,
}

impl<'a> Converter<'a> {
    pub fn new(
        input: &'a Path,
        output: &'a Path,
        output_fmt: &'a OutputFmt,
        datatype: &'a DataType,
    ) -> Self {
        Self {
            input,
            output,
            output_fmt,
            datatype,
            is_dir: false,
        }
    }

    pub fn convert_unsorted(&mut self, input_fmt: &InputFmt) {
        let (matrix, header) = self.get_sequence(input_fmt);
        self.convert(&matrix, header);
    }

    pub fn convert_sorted(&mut self, input_fmt: &InputFmt) {
        let (mut matrix, header) = self.get_sequence(input_fmt);
        matrix.sort_keys();
        self.convert(&matrix, header)
    }

    pub fn set_isdir(&mut self, is_dir: bool) {
        self.is_dir = is_dir;
    }

    fn get_sequence(&mut self, input_fmt: &InputFmt) -> (IndexMap<String, String>, Header) {
        let seq = Sequence::new(self.input, self.datatype);
        seq.get(input_fmt)
    }

    fn convert(&self, matrix: &IndexMap<String, String>, header: Header) {
        let mut convert = SeqWriter::new(self.output, matrix, &header, None, &PartitionFmt::None);
        convert
            .write_sequence(self.output_fmt)
            .expect("Failed writing output files");

        if !self.is_dir {
            convert.print_save_path();
        }
    }
}
