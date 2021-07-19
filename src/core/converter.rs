use std::path::Path;

use indexmap::IndexMap;

use crate::writer::seqwriter::SeqWriter;

use crate::helper::common::{DataType, Header, InputFmt, OutputFmt, PartitionFmt};
use crate::helper::sequence::Sequence;

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
        let seq = self.get_sequence(input_fmt);
        self.convert(&seq.matrix, seq.header);
    }

    pub fn convert_sorted(&mut self, input_fmt: &InputFmt) {
        let mut seq = self.get_sequence(input_fmt);
        seq.matrix.sort_keys();
        self.convert(&seq.matrix, seq.header)
    }

    pub fn set_isdir(&mut self, is_dir: bool) {
        self.is_dir = is_dir;
    }

    fn get_sequence(&mut self, input_fmt: &InputFmt) -> Sequence {
        let mut sequence = Sequence::new();
        sequence.get(self.input, input_fmt, self.datatype);
        sequence
    }

    fn convert(&self, matrix: &IndexMap<String, String>, header: Header) {
        let mut convert = SeqWriter::new(self.output, matrix, header, None, &PartitionFmt::None);
        convert
            .write_sequence(self.output_fmt)
            .expect("Failed writing output files");

        if !self.is_dir {
            convert
                .print_save_path()
                .expect("Failed writing save path to stdout");
        }
    }
}
