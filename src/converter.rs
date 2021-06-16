use std::path::Path;

use indexmap::IndexMap;

use crate::writer::SeqWriter;

use crate::common::{Header, OutputFormat, PartitionFormat};
use crate::fasta::Fasta;

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
        let outpath = self.output.join(self.input.file_stem().unwrap());
        self.convert(&outpath, &fasta.matrix, header)
    }

    fn convert(&self, outpath: &Path, matrix: &IndexMap<String, String>, header: Header) {
        let mut convert = SeqWriter::new(outpath, matrix, header, None, &PartitionFormat::None);
        match self.output_format {
            OutputFormat::Nexus => convert.write_sequence(self.output_format),
            OutputFormat::Phylip => convert.write_sequence(self.output_format),
            _ => (),
        }
    }
}
