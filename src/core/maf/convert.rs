//! Convert MAF to a different alignment format.
//!
//! Include support to match name with a reference sequence.
use std::{
    collections::HashMap,
    fs::{self, File},
    io::BufReader,
    path::{Path, PathBuf},
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
        bed::{BedParser, BedRecord},
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
        if self.name_from_bed {
            self.parse_maf_from_bed();
        } else {
            unreachable!(
                "Name source is not supported. \
            Use BED file instead and set the flag --from-bed"
            );
        }
        self.print_output_info();
    }

    fn parse_maf_from_bed(&self) {
        let spin = utils::set_spinner();
        spin.set_message("Converting MAF format...");
        let names = self.get_name_from_bed(self.name_source);
        self.input_files.par_iter().for_each(|file| {
            let file = File::open(file).expect("Unable to open file");
            let buff = BufReader::new(file);
            let maf = MafReader::new(buff);
            // Take the gene name from the BED file as a key
            // and store the alignment in a SeqMatrix
            let mut aln_collection: HashMap<String, SeqMatrix> = HashMap::new();
            let mut missing_refs = HashMap::new();
            maf.into_iter().for_each(|paragraph| match paragraph {
                MafParagraph::Alignment(aln) => {
                    let new_matrix = self.convert_to_seqmatrix(&aln);
                    // Target is usually the first sequence
                    let target = &aln.sequences[0];
                    let bed_name = names.get(&target.start);
                    match bed_name {
                        Some(name) => match aln_collection.get_mut(name) {
                            Some(matrix) => {
                                matrix.extend(new_matrix.to_owned());
                            }
                            None => {
                                aln_collection.insert(name.to_string(), new_matrix);
                            }
                        },
                        None => {
                            let name =
                                format!("{}-{}-{}", target.source, target.start, target.size,);
                            missing_refs.insert(name.replace(".", "_"), aln);
                        }
                    }
                }
                _ => (),
            });

            spin.set_message("Writing output files...");
            aln_collection.par_iter().for_each(|(target, matrix)| {
                let output = self.generate_output_path(self.output_dir, Path::new(target));
                let header = self.get_header(matrix);
                self.write_matrix(matrix, &header, &output);
            });
            spin.finish_with_message("Finished converting MAF format!\n");

            if !missing_refs.is_empty() {
                self.write_missing_refs(&missing_refs);
            }
        });
    }

    fn write_missing_refs(&self, missing_refs: &HashMap<String, MafAlignment>) {
        log::warn!(
            "{}: {}\n",
            "Missing references".yellow(),
            missing_refs.len()
        );
        let output_dir = self.output_dir.join("missing-refs");
        fs::create_dir_all(&output_dir).expect("Failed to create output directory");
        missing_refs.par_iter().for_each(|(name, aln)| {
            let output = self.generate_output_path(&output_dir, Path::new(name));
            let matrix = self.convert_to_seqmatrix(aln);
            let header = self.get_header(&matrix);
            self.write_matrix(&matrix, &header, &output);
        });
    }

    fn generate_output_path(&self, output_dir: &Path, output_file: &Path) -> PathBuf {
        let output_fname = files::create_output_fname(output_dir, output_file, self.output_fmt);
        output_fname
    }

    fn convert_to_seqmatrix(&self, alignments: &MafAlignment) -> SeqMatrix {
        let mut matrix: SeqMatrix = IndexMap::new();
        alignments.sequences.iter().for_each(|sample| {
            let seq = String::from_utf8_lossy(&sample.text).to_string();
            matrix.insert(sample.source.to_string(), seq);
        });
        matrix
    }

    fn get_header(&self, matrix: &SeqMatrix) -> Header {
        let mut header = Header::new();
        header.from_seq_matrix(matrix, true);
        header
    }

    fn write_matrix(&self, matrix: &SeqMatrix, header: &Header, output: &Path) {
        let mut writer = SeqWriter::new(output, matrix, header);
        writer
            .write_sequence(self.output_fmt)
            .expect("Failed writing output files");
    }

    fn get_name_from_bed(&self, bed: &Path) -> HashMap<usize, String> {
        let mut bed = BedParser::new(bed);
        let bed = bed.parse().expect("Unable to parse BED file");
        // Create a hashmap with the start position as key
        // and the gene name as value
        let mut names = HashMap::new();
        bed.iter().for_each(|record| {
            let gene_name = self.format_bed_name(record);
            names.insert(record.chrom_start, gene_name);
        });
        names
    }

    fn format_bed_name(&self, record: &BedRecord) -> String {
        let name = match &record.name {
            Some(name) => {
                format!(
                    "{}-{}-{}-{}",
                    name, record.chrom, record.chrom_start, record.chrom_end
                )
            }
            None => {
                format!(
                    "{}-{}-{}",
                    record.chrom, record.chrom_start, record.chrom_end
                )
            }
        };
        // Replace any dots in the name with underscores
        // to avoid issues with file extensions
        name.replace(".", "_")
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
    use crate::{helper::types::DnaStrand, parser::maf::MafSequence};

    #[test]
    fn test_convert_to_seqmatrix() {
        let aln = MafAlignment {
            score: None,
            quality: None,
            information: None,
            empty: None,
            sequences: vec![
                MafSequence {
                    source: "seq1".to_string(),
                    start: 1,
                    size: 10,
                    strand: DnaStrand::Forward,
                    src_size: 100,
                    text: b"ACGTACGTAC".to_vec(),
                },
                MafSequence {
                    source: "seq2".to_string(),
                    start: 1,
                    size: 10,
                    strand: DnaStrand::Forward,
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
        let matrix = converter.convert_to_seqmatrix(&aln);
        assert_eq!(matrix.len(), 2);
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
        let output_dir = Path::new("output");
        let output = converter.generate_output_path(output_dir, &Path::new("test.maf"));
        assert_eq!(output, Path::new("output/test.fas"));
        let gene_name = String::from("hoxa1-chrom1-1-10");
        let output = converter.generate_output_path(output_dir, Path::new(&gene_name));
        assert_eq!(output, Path::new("output/hoxa1-chrom1-1-10.fas"));
    }

    #[test]
    fn test_format_bed_name() {
        let record = BedRecord::new("chr1".to_string(), 1, 10, Some("gene".to_string()));
        let converter = MafConverter {
            input_files: &vec![],
            name_source: &Path::new(""),
            name_from_bed: false,
            output_fmt: &OutputFmt::Fasta,
            output_dir: &Path::new("output"),
        };
        let name = converter.format_bed_name(&record);
        assert_eq!(name, "gene-chr1-1-10");
    }
}
