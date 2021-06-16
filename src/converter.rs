use std::path::Path;

use indexmap::IndexMap;

use crate::writer::SeqWriter;

use crate::common::{Header, OutputFormat, PartitionFormat};
use crate::fasta::Fasta;
use crate::nexus::Nexus;
use crate::phylip::Phylip;

pub struct Converter<'a> {
    input: &'a Path,
    output: &'a Path,
    output_format: &'a OutputFormat,
}

impl<'a> Converter<'a> {
    pub fn new(input: &'a Path, output: &'a Path, output_format: &'a OutputFormat) -> Self {
        Self {
            input,
            output,
            output_format,
        }
    }

    pub fn convert_fasta(&mut self) {
        let mut fasta = Fasta::new(self.input);
        fasta.read();
        let header = fasta.get_header();
        self.convert(&fasta.matrix, header)
    }

    pub fn convert_nexus(&mut self) {
        let mut nex = Nexus::new(self.input);
        nex.read().expect("CANNOT READ NEXUS FILES");
        let header = nex.get_header();
        self.convert(&nex.matrix, header);
    }

    pub fn convert_phylip(&mut self) {
        let input = Path::new(self.input);
        let mut phy = Phylip::new(input);
        phy.read().expect("CANNOT READ PHYLIP FILES");
        let header = phy.get_header();
        self.convert(&phy.matrix, header);
    }

    fn convert(&self, matrix: &IndexMap<String, String>, header: Header) {
        let save_path = self.output.join(self.input.file_stem().unwrap());
        let mut convert = SeqWriter::new(&save_path, matrix, header, None, &PartitionFormat::None);
        match self.output_format {
            OutputFormat::Nexus => convert.write_sequence(self.output_format),
            OutputFormat::Phylip => convert.write_sequence(self.output_format),
            OutputFormat::Fasta => convert.write_fasta(),
        }
    }
}
