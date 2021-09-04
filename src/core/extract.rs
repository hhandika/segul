use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use ansi_term::Colour::Yellow;
use indexmap::IndexMap;
use rayon::prelude::*;
use regex::Regex;

use crate::helper::sequence::{SeqCheck, Sequence};
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, PartitionFmt, SeqMatrix};
use crate::helper::utils;
use crate::writer::sequences::SeqWriter;

pub enum Params {
    Regex(String),
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
        let file_counts = AtomicUsize::new(0);
        let spin = utils::set_spinner();
        spin.set_message("Extracting sequences with matching IDs...");
        files.par_iter().for_each(|file| {
            let (seq, _) = Sequence::new(file, self.datatype).get(self.input_fmt);
            let matrix = self.get_matrix(seq);
            if !matrix.is_empty() {
                let header = self.get_header(&matrix);
                let outname = self.get_output_names(output, file, output_fmt);
                let mut writer =
                    SeqWriter::new(&outname, &matrix, &header, None, &PartitionFmt::None);
                writer
                    .write_sequence(output_fmt)
                    .expect("Failed writing the output files");
                file_counts.fetch_add(1, Ordering::Relaxed);
            }
        });
        spin.finish_with_message("Finished extracting sequences!\n");
        let counts = file_counts.load(Ordering::Relaxed);
        assert!(counts > 0, "No matching IDs found!");
        self.print_output_info(counts, output);
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
            Params::Id(ids) => seqmat.iter().for_each(|(id, seq)| {
                ids.iter().for_each(|match_id| {
                    if match_id == id {
                        matrix.insert(id.to_string(), seq.to_string());
                    }
                })
            }),
            Params::None => panic!("Please, specify a matching parameter!"),
        };
        matrix
    }

    fn get_output_names(&self, dir: &Path, file: &Path, output_fmt: &OutputFmt) -> PathBuf {
        let path = dir.join(
            file.file_name()
                .expect("Failed parsing filename for output file"),
        );
        match output_fmt {
            OutputFmt::Fasta | OutputFmt::FastaInt => path.with_extension("fas"),
            OutputFmt::Nexus | OutputFmt::NexusInt => path.with_extension("nex"),
            OutputFmt::Phylip | OutputFmt::PhylipInt => path.with_extension("phy"),
        }
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

    fn print_output_info(&self, file_counts: usize, output: &Path) {
        log::info!("{}", Yellow.paint("Output"));
        log::info!("{:18}: {}", "File counts", file_counts);
        log::info!("{:18}: {}", "Output dir", output.display())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_match_id() {
        let id = "Bunomys_penitus";
        let re = "(?i)(Penitus)$";
        let extract = Extract::new(&Params::None, &InputFmt::Fasta, &DataType::Dna);
        assert_eq!(true, extract.match_id(id, re));
    }
}
