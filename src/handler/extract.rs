use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use colored::Colorize;
use indexmap::IndexMap;
use rayon::prelude::*;
use regex::Regex;

use crate::handler::OutputPrint;
use crate::helper::sequence::{SeqCheck, SeqParser};
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix};
use crate::helper::{filenames, utils};
use crate::writer::sequences::SeqWriter;

impl OutputPrint for Extract<'_> {}

pub enum ExtractOpts {
    Regex(String),
    Id(Vec<String>),
    None,
}

pub struct Extract<'a> {
    input_fmt: &'a InputFmt,
    opts: &'a ExtractOpts,
    datatype: &'a DataType,
}

impl<'a> Extract<'a> {
    pub fn new(opts: &'a ExtractOpts, input_fmt: &'a InputFmt, datatype: &'a DataType) -> Self {
        Self {
            opts,
            input_fmt,
            datatype,
        }
    }

    pub fn extract_sequences(&self, files: &[PathBuf], output: &Path, output_fmt: &OutputFmt) {
        let file_counts = AtomicUsize::new(0);
        let spin = utils::set_spinner();
        spin.set_message("Extracting sequences with matching IDs...");
        files.par_iter().for_each(|file| {
            let (seq, _) = SeqParser::new(file, self.datatype).get(self.input_fmt);
            let matrix = self.get_matrix(seq);
            if !matrix.is_empty() {
                let header = self.get_header(&matrix);
                let outname = filenames::create_output_fname(output, file, output_fmt);
                let mut writer = SeqWriter::new(&outname, &matrix, &header);
                writer
                    .write_sequence(output_fmt)
                    .expect("Failed writing the output files");
                file_counts.fetch_add(1, Ordering::Relaxed);
            }
        });
        spin.finish_with_message("Finished extracting sequences!\n");
        let counts = file_counts.load(Ordering::Relaxed);
        assert!(counts > 0, "No matching IDs found!");
        self.print_output_info(&counts, output, output_fmt);
    }

    fn get_matrix(&self, seqmat: SeqMatrix) -> SeqMatrix {
        let mut matrix: SeqMatrix = IndexMap::new();
        match self.opts {
            ExtractOpts::Regex(re) => seqmat.iter().for_each(|(id, seq)| {
                let matched_id = self.match_id(id, re);
                if matched_id {
                    matrix.insert(id.to_string(), seq.to_string());
                }
            }),
            ExtractOpts::Id(ids) => seqmat.iter().for_each(|(id, seq)| {
                ids.iter().for_each(|match_id| {
                    if match_id == id {
                        matrix.insert(id.to_string(), seq.to_string());
                    }
                })
            }),
            ExtractOpts::None => panic!("Please, specify a matching parameter!"),
        };
        matrix
    }

    fn match_id(&self, id: &str, re: &str) -> bool {
        let re = Regex::new(re).expect("Failed capturing nexus commands");
        re.is_match(id)
    }

    fn get_header(&self, matrix: &SeqMatrix) -> Header {
        let mut seq_info = SeqCheck::new();
        seq_info.check(matrix);
        let mut header = Header::new();
        header.aligned = seq_info.is_alignment;
        header.nchar = seq_info.longest;
        header.ntax = matrix.len();
        header
    }

    fn print_output_info(&self, file_counts: &usize, output: &Path, output_fmt: &OutputFmt) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "File counts", utils::fmt_num(file_counts));
        log::info!("{:18}: {}", "Output dir", output.display());
        self.print_output_fmt(output_fmt);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_match_id() {
        let id = "Bunomys_penitus";
        let re = "(?i)(Penitus)$";
        let extract = Extract::new(&ExtractOpts::None, &InputFmt::Fasta, &DataType::Dna);
        assert_eq!(true, extract.match_id(id, re));
    }

    #[test]
    fn test_get_matrix_re() {
        let re = ExtractOpts::Regex(String::from("(?i)(celebensis)"));
        let file = Path::new("tests/files/complete.nex");
        let extract = Extract::new(&re, &InputFmt::Nexus, &DataType::Dna);
        let (seq, _) = SeqParser::new(file, extract.datatype).get(extract.input_fmt);
        let matrix = extract.get_matrix(seq);
        assert_eq!(2, matrix.len());
    }

    #[test]
    fn test_get_matrix_id() {
        let re = ExtractOpts::Id(vec![String::from("Taeromys_calitrichus_NMVZ27408")]);
        let file = Path::new("tests/files/complete.nex");
        let extract = Extract::new(&re, &InputFmt::Nexus, &DataType::Dna);
        let (seq, _) = SeqParser::new(file, extract.datatype).get(extract.input_fmt);
        let matrix = extract.get_matrix(seq);
        assert_eq!(1, matrix.len());
    }
}
