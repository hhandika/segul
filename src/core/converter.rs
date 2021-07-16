use std::path::Path;

use indexmap::IndexMap;

use crate::writer::seqwriter::SeqWriter;

use crate::helper::common::{self, Header, InputFmt, OutputFmt, PartitionFormat};
use crate::parser::fasta::Fasta;
use crate::parser::nexus::Nexus;
use crate::parser::phylip::Phylip;

pub struct Converter<'a> {
    input: &'a Path,
    output: &'a Path,
    output_fmt: &'a OutputFmt,
    is_dir: bool,
}

impl<'a> Converter<'a> {
    pub fn new(input: &'a Path, output: &'a Path, output_fmt: &'a OutputFmt, is_dir: bool) -> Self {
        Self {
            input,
            output,
            output_fmt,
            is_dir,
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

    fn get_sequence(&mut self, input_fmt: &InputFmt) -> (IndexMap<String, String>, Header) {
        match input_fmt {
            InputFmt::Fasta => self.convert_fasta(),
            InputFmt::Nexus => self.convert_nexus(),
            InputFmt::Phylip => self.convert_phylip(false),
            InputFmt::PhylipInt => self.convert_phylip(true),
            InputFmt::Auto => {
                let input_fmt = common::infer_input_auto(self.input);
                self.get_sequence(&input_fmt)
            }
        }
    }

    fn convert_fasta(&mut self) -> (IndexMap<String, String>, Header) {
        let mut fas = Fasta::new(self.input);
        fas.read();
        (fas.matrix, fas.header)
    }

    fn convert_nexus(&mut self) -> (IndexMap<String, String>, Header) {
        let mut nex = Nexus::new(self.input);
        nex.read().expect("CANNOT READ NEXUS FILES");
        (nex.matrix, nex.header)
    }

    fn convert_phylip(&mut self, interleave: bool) -> (IndexMap<String, String>, Header) {
        let input = Path::new(self.input);
        let mut phy = Phylip::new(input, interleave);
        phy.read().expect("CANNOT READ PHYLIP FILES");
        (phy.matrix, phy.header)
    }

    fn convert(&self, matrix: &IndexMap<String, String>, header: Header) {
        let mut convert = SeqWriter::new(self.output, matrix, header, None, &PartitionFormat::None);
        convert
            .write_sequence(self.output_fmt)
            .expect("CANNOT WRITE OUTPUT FILES");

        if !self.is_dir {
            convert.print_save_path();
        }
    }
}
