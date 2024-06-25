//! Extract sequences from a collection of sequence files.
//!
//! Support extraction of sequences based on matching IDs or regular expressions.
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use colored::Colorize;
use indexmap::IndexMap;
use rayon::prelude::*;
use regex::Regex;

use crate::handler::OutputPrint;
use crate::helper::sequence::{SeqCheck, SeqParser};
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix};
use crate::helper::{files, utils};
use crate::writer::sequences::SeqWriter;

impl OutputPrint for Extract<'_> {}

/// Extraction options
pub enum ExtractOpts {
    /// Use regular expression to match sequence IDs
    Regex(String),
    /// Match sequence IDs from a list of IDs
    Id(Vec<String>),
    None,
}

/// Extract sequences from a collection of sequence files.
/// See each type for available options.
pub struct Extract<'a> {
    /// Input format of the sequence files
    input_fmt: &'a InputFmt,
    /// Extraction options
    opts: &'a ExtractOpts,
    /// Data type of the sequences
    datatype: &'a DataType,
    /// Output directory
    output_dir: &'a Path,
    /// Output format of the sequence files
    output_fmt: &'a OutputFmt,
}

impl<'a> Extract<'a> {
    /// Create a new `Extract` instance
    pub fn new(
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        opts: &'a ExtractOpts,
        output_dir: &'a Path,
        output_fmt: &'a OutputFmt,
    ) -> Self {
        Self {
            input_fmt,
            datatype,
            opts,
            output_dir,
            output_fmt,
        }
    }

    /// Extract sequences with matching IDs or regular expressions
    pub fn extract_sequences(&self, files: &[PathBuf]) {
        let file_counts = AtomicUsize::new(0);
        let spin = utils::set_spinner();
        spin.set_message("Extracting sequences with matching IDs...");
        files.par_iter().for_each(|file| {
            let (seq, _) = SeqParser::new(file, self.datatype).parse(self.input_fmt);
            let matrix = self.get_matrix(seq);
            if !matrix.is_empty() {
                let header = self.get_header(&matrix);
                let output_name =
                    files::create_output_fname(self.output_dir, file, self.output_fmt);
                let mut writer = SeqWriter::new(&output_name, &matrix, &header);
                writer
                    .write_sequence(self.output_fmt)
                    .expect("Failed writing the output files");
                file_counts.fetch_add(1, Ordering::Relaxed);
            }
        });
        spin.finish_with_message("Finished extracting sequences!\n");
        let counts = file_counts.load(Ordering::Relaxed);
        assert!(counts > 0, "No matching IDs found!");
        self.print_output_info(&counts, self.output_dir, self.output_fmt);
    }

    fn get_matrix(&self, matrix: SeqMatrix) -> SeqMatrix {
        let mut new_matrix: SeqMatrix = IndexMap::new();
        match self.opts {
            ExtractOpts::Regex(re) => matrix.iter().for_each(|(id, seq)| {
                let matched_id = self.match_id(id, re);
                if matched_id {
                    new_matrix.insert(id.to_string(), seq.to_string());
                }
            }),
            ExtractOpts::Id(ids) => matrix.iter().for_each(|(id, seq)| {
                ids.iter().for_each(|match_id| {
                    if match_id == id {
                        new_matrix.insert(id.to_string(), seq.to_string());
                    }
                })
            }),
            ExtractOpts::None => panic!("Please, specify a matching parameter!"),
        };
        new_matrix
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
        let output_dir = Path::new("tests/files");
        let output_fmt = OutputFmt::Fasta;
        let extract = Extract::new(
            &InputFmt::Fasta,
            &DataType::Dna,
            &ExtractOpts::None,
            &output_dir,
            &output_fmt,
        );
        assert!(extract.match_id(id, re));
    }

    #[test]
    fn test_get_matrix_re() {
        let re = ExtractOpts::Regex(String::from("(?i)(celebensis)"));
        let file = Path::new("tests/files/complete.nex");
        let output_dir = Path::new("tests/files");
        let output_fmt = OutputFmt::Fasta;
        let extract = Extract::new(
            &InputFmt::Nexus,
            &DataType::Dna,
            &re,
            output_dir,
            &output_fmt,
        );
        let (seq, _) = SeqParser::new(file, extract.datatype).parse(extract.input_fmt);
        let matrix = extract.get_matrix(seq);
        assert_eq!(2, matrix.len());
    }

    #[test]
    fn test_get_matrix_id() {
        let re = ExtractOpts::Id(vec![String::from("Taeromys_calitrichus_NMVZ27408")]);
        let file = Path::new("tests/files/complete.nex");
        let output_dir = Path::new("tests/files");
        let output_fmt = OutputFmt::Fasta;
        let extract = Extract::new(
            &InputFmt::Nexus,
            &DataType::Dna,
            &re,
            output_dir,
            &output_fmt,
        );
        let (seq, _) = SeqParser::new(file, extract.datatype).parse(extract.input_fmt);
        let matrix = extract.get_matrix(seq);
        assert_eq!(1, matrix.len());
    }
}
