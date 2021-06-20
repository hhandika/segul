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

    pub fn convert_fasta(&mut self) {
        let mut fas = Fasta::new(self.input);
        fas.read();
        self.convert(&fas.matrix, fas.header)
    }

    pub fn convert_nexus(&mut self) {
        let mut nex = Nexus::new(self.input);
        nex.read().expect("CANNOT READ NEXUS FILES");
        self.convert(&nex.matrix, nex.header);
    }

    pub fn convert_phylip(&mut self, interleave: bool) {
        let input = Path::new(self.input);
        let mut phy = Phylip::new(input, interleave);
        phy.read().expect("CANNOT READ PHYLIP FILES");
        self.convert(&phy.matrix, phy.header);
    }

    fn convert(&self, matrix: &IndexMap<String, String>, header: Header) {
        let save_path = self.output.join(self.input.file_stem().unwrap());
        let mut convert = SeqWriter::new(&save_path, matrix, header, None, &PartitionFormat::None);
        match self.output_format {
            SeqFormat::Nexus => convert.write_sequence(self.output_format),
            SeqFormat::Phylip => convert.write_sequence(self.output_format),
            SeqFormat::Fasta => convert.write_fasta(),
        }

        if !self.is_dir {
            convert.display_save_path();
        }
    }
}
