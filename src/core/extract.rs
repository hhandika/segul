use std::fs;
use std::path::{Path, PathBuf};

use ansi_term::Colour::Yellow;
use indexmap::IndexMap;
use regex::Regex;

use crate::helper::sequence::{SeqCheck, Sequence};
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, PartitionFmt, SeqMatrix};
use crate::helper::utils;
use crate::writer::sequences::SeqWriter;

pub enum Params {
    Regex(String),
    File(PathBuf),
    Id(Vec<String>),
    None,
}

pub struct Extract<'a> {
    input_fmt: &'a InputFmt,
    params: &'a Params,
    datatype: &'a DataType,
}

impl<'a> Extract<'a> {
    pub fn new(params: &'a Params, input_fmt: &'a InputFmt, datatype: &'a DataType) -> Self {
        Self {
            params,
            input_fmt,
            datatype,
        }
    }

    pub fn extract_sequences(&self, files: &[PathBuf], output: &Path, output_fmt: &OutputFmt) {
        fs::create_dir_all(output).expect("Failed creating output directory");
        let mut file_counts = 0;
        let spin = utils::set_spinner();
        spin.set_message("Extracting sequences with matching IDs...");
        files.iter().for_each(|file| {
            let (seq, _) = Sequence::new(file, self.datatype).get(self.input_fmt);
            let matrix = 
            if !matrix.is_empty() {
                let header = self.get_header(&matrix);
                let outname = output.join(
                    file.file_name()
                        .expect("Failed parsing filename for output file"),
                );
                let mut writer =
                    SeqWriter::new(&outname, &matrix, &header, None, &PartitionFmt::None);
                writer
                    .write_sequence(output_fmt)
                    .expect("Failed writing the output files");
                file_counts += 1;
            }
        });
        spin.finish_with_message("Finished finding matching IDs!\n");
        assert!(file_counts > 0, "No matching IDs found!");
        self.print_output_info(file_counts);
    }

    fn match_id(&self, id: &str, re: &str) -> bool {
        let re = Regex::new(re).expect("Failed capturing nexus commands");
        re.is_match(id)
    }

    fn get_matrix(&self, seqmat: SeqMatrix) -> SeqMatrix {
        let mut matrix: SeqMatrix = IndexMap::new();
        match self.params {
            Params::Regex(re) => seqmat.iter().for_each(|(id, seq)| {
                let matched_id = self.match_id(id, re);
                if matched_id {
                    matrix.insert(id.to_string(), seq.to_string());
                }
            }),
            Params::File(path) => println!("Path: {}\n", path.display()),
            Params::Id(ids) => mat.iter().for_each(|(id, seq)| {
                ids.iter().for_each(|match_id| {
                    if match_id == id {
                        matrix.insert(id.to_string(), seq.to_string());
                    }
                })
            }),
            _ => unreachable!("Please, specify a matching parameter!"),
        };
        matrix
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

    fn print_output_info(&self, file_counts: usize) {
        log::info!("{}", Yellow.paint("Output"));
        log::info!("{:18}: {}", "File counts", file_counts);
    }
}
