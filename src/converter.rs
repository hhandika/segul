use std::path::Path;

use indexmap::IndexMap;

use crate::writer::SeqWriter;

use crate::common::{Header, OutputFormat, PartitionFormat};
use crate::fasta::Fasta;

pub fn convert_fasta(input: &Path, output: &Path, output_format: OutputFormat) {
    let mut fasta = Fasta::new(input);
    fasta.read();
    let header = fasta.get_header();
    let outpath = output.join(input.file_stem().unwrap());
    convert(&outpath, &output_format, &fasta.matrix, header)
}

fn convert(
    outpath: &Path,
    output_format: &OutputFormat,
    matrix: &IndexMap<String, String>,
    header: Header,
) {
    let mut convert = SeqWriter::new(outpath, matrix, header, None, &PartitionFormat::None);
    match output_format {
        OutputFormat::Nexus => convert.write_sequence(output_format),
        OutputFormat::Phylip => convert.write_sequence(output_format),
        _ => (),
    }
}
