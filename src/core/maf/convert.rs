//! Convert MAF to a different alignment format.
//!
//! Include support to match name with a reference sequence.
use std::{
    collections::HashMap,
    fs::File,
    io::BufReader,
    path::{Path, PathBuf},
    sync::Mutex,
};

use colored::Colorize;
use indexmap::IndexMap;
use rayon::prelude::*;

use crate::{
    core::OutputPrint,
    helper::{
        files,
        types::{Header, OutputFmt, SeqMatrix},
        utils,
    },
    parser::{
        bed::BedParser,
        maf::{MafAlignment, MafParagraph, MafReader},
    },
    writer::sequences::SeqWriter,
};

pub struct MafConverter<'a> {
    input_files: &'a [PathBuf],
    name_source: &'a Path,
    name_from_bed: bool,
    output_dir: &'a Path,
    output_fmt: &'a OutputFmt,
}

impl OutputPrint for MafConverter<'_> {}

impl<'a> MafConverter<'a> {
    pub fn new(
        input_files: &'a [PathBuf],
        name_source: &'a Path,
        name_from_bed: bool,
        output_dir: &'a Path,
        output_fmt: &'a OutputFmt,
    ) -> Self {
        Self {
            input_files,
            name_source,
            name_from_bed,
            output_fmt,
            output_dir,
        }
    }

    pub fn convert(&self) {
        let spin = utils::set_spinner();
        spin.set_message("Converting MAF format...");
        if self.name_from_bed {
            self.parse_maf_from_bed();
        }
        spin.finish_with_message("Finished converting MAF format!\n");
        self.print_output_info();
    }

    fn parse_maf_from_bed(&self) {
        let names = self.get_name_from_bed(self.name_source);
        let parallel_names = Mutex::new(names);
        self.input_files.par_iter().for_each(|file| {
            let file = File::open(file).expect("Unable to open file");
            let buff = BufReader::new(file);
            let maf = MafReader::new(buff);
            maf.into_iter().for_each(|paragraph| match paragraph {
                MafParagraph::Alignment(aln) => {
                    let (matrix, header) = self.convert_to_seqmatrix(&aln);
                    // Target is usually the first sequence
                    let locked_names = parallel_names.lock().expect("Fail to access bed names");
                    let file_name = locked_names
                        .get(&aln.sequences[0].start)
                        .map(|name| Path::new(name))
                        .unwrap_or_else(|| Path::new(&aln.sequences[0].source));
                    let output = self.generate_output_path(file_name);
                    self.write_matrix(&matrix, &header, &output);
                }
                _ => (),
            });
        });
    }

    fn generate_output_path(&self, output_file: &Path) -> PathBuf {
        let output_fname =
            files::create_output_fname(self.output_dir, output_file, self.output_fmt);
        output_fname
    }

    fn convert_to_seqmatrix(&self, alignments: &MafAlignment) -> (SeqMatrix, Header) {
        let mut matrix: SeqMatrix = IndexMap::new();
        alignments.sequences.iter().for_each(|sample| {
            let seq = String::from_utf8_lossy(&sample.text).to_string();
            matrix.insert(sample.source.to_string(), seq);
        });
        let mut header = Header::new();
        header.from_seq_matrix(&matrix, true);
        (matrix, header)
    }

    fn write_matrix(&self, matrix: &SeqMatrix, header: &Header, output: &Path) {
        let mut writer = SeqWriter::new(output, matrix, header);
        writer
            .write_sequence(self.output_fmt)
            .expect("Failed writing output files");
    }

    fn get_name_from_bed(&self, bed: &Path) -> HashMap<usize, String> {
        let bed = BedParser::new(bed, false);
        let bed = bed.parser().expect("Unable to parse BED file");
        let mut names = HashMap::new();
        bed.iter().for_each(|record| match &record.name {
            Some(name) => {
                names.insert(record.chrom_start, name.to_string());
            }
            None => {
                names.insert(record.chrom_start, record.chrom.to_owned());
            }
        });
        names
    }

    fn print_output_info(&self) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Output dir", self.output_dir.display());
        self.print_output_fmt(self.output_fmt);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::maf::MafSequence;

    #[test]
    fn test_convert_to_seqmatrix() {
        let aln = MafAlignment {
            score: 0.0,
            quality: None,
            information: None,
            empty: None,
            sequences: vec![
                MafSequence {
                    source: "seq1".to_string(),
                    start: 1,
                    size: 10,
                    strand: '+',
                    src_size: 100,
                    text: b"ACGTACGTAC".to_vec(),
                },
                MafSequence {
                    source: "seq2".to_string(),
                    start: 1,
                    size: 10,
                    strand: '+',
                    src_size: 100,
                    text: b"ACGTACGTAC".to_vec(),
                },
            ],
        };
        let converter = MafConverter {
            input_files: &vec![],
            name_source: &Path::new(""),
            name_from_bed: false,
            output_fmt: &OutputFmt::Fasta,
            output_dir: &Path::new(""),
        };
        let (matrix, header) = converter.convert_to_seqmatrix(&aln);
        assert_eq!(matrix.len(), 2);
        assert_eq!(header.ntax, 2);
        assert_eq!(header.nchar, 10);
        assert_eq!(header.datatype, String::from("dna"));
    }

    #[test]
    fn test_generate_output_path() {
        let converter = MafConverter {
            input_files: &vec![],
            name_source: &Path::new(""),
            name_from_bed: false,
            output_fmt: &OutputFmt::Fasta,
            output_dir: &Path::new("output"),
        };
        let output = converter.generate_output_path(&Path::new("test.maf"));
        assert_eq!(output, Path::new("output/test.fas"));
    }
}
