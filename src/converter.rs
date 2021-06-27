use std::path::Path;

use indexmap::IndexMap;

use crate::writer::SeqWriter;

use crate::common::{Header, PartitionFormat, SeqFormat};
use crate::fasta::Fasta;
use crate::nexus::Nexus;
use crate::phylip::Phylip;

pub struct Converter<'a> {
    input: &'a Path,
    output: &'a Path,
    output_format: &'a SeqFormat,
    is_dir: bool,
}

impl<'a> Converter<'a> {
    pub fn new(
        input: &'a Path,
        output: &'a Path,
        output_format: &'a SeqFormat,
        is_dir: bool,
    ) -> Self {
        Self {
            input,
            output,
            output_format,
            is_dir,
        }
    }

    pub fn convert_unsorted(&mut self, input_format: &SeqFormat) {
        let (matrix, header) = self.get_sequence(input_format);
        self.convert(&matrix, header);
    }

    pub fn convert_sorted(&mut self, input_format: &SeqFormat) {
        let (mut matrix, header) = self.get_sequence(input_format);
        matrix.sort_keys();
        self.convert(&matrix, header)
    }

    fn get_sequence(&mut self, input_format: &SeqFormat) -> (IndexMap<String, String>, Header) {
        match input_format {
            SeqFormat::Fasta | SeqFormat::FastaInt => self.convert_fasta(),
            SeqFormat::Nexus | SeqFormat::NexusInt => self.convert_nexus(),
            SeqFormat::Phylip => self.convert_phylip(false),
            SeqFormat::PhylipInt => self.convert_phylip(true),
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
        let save_path = self.output.join(self.input.file_stem().unwrap());
        let mut convert = SeqWriter::new(&save_path, matrix, header, None, &PartitionFormat::None);
        convert
            .write_sequence(self.output_format)
            .expect("CANNOT WRITE OUTPUT FILES");

        if !self.is_dir {
            convert.display_save_path();
        }
    }
}
