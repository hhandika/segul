use std::path::{Path, PathBuf};

use colored::Colorize;
use rayon::prelude::*;
use regex::Regex;

use crate::handler::OutputPrint;
use crate::helper::filenames;
use crate::helper::finder::IDs;
use crate::helper::sequence::{SeqCheck, SeqParser};
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix};
use crate::helper::utils;
use crate::writer::sequences::SeqWriter;

impl OutputPrint for Remove<'_> {}

pub enum RemoveOpts {
    Id(Vec<String>),
    Regex(String),
}

pub struct Remove<'a> {
    input_fmt: &'a InputFmt,
    datatype: &'a DataType,
    outdir: &'a Path,
    output_fmt: &'a OutputFmt,
    opts: &'a RemoveOpts,
}

impl<'a> Remove<'a> {
    pub fn new(
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        outdir: &'a Path,
        output_fmt: &'a OutputFmt,
        opts: &'a RemoveOpts,
    ) -> Self {
        Self {
            input_fmt,
            datatype,
            outdir,
            output_fmt,
            opts,
        }
    }

    pub fn remove(&self, files: &[PathBuf]) {
        let spin = utils::set_spinner();
        spin.set_message("Removing sequences...");
        match self.opts {
            RemoveOpts::Id(ids) => self.par_remove(files, ids),
            RemoveOpts::Regex(re) => {
                let ids = self.find_matching_ids(files, re);
                self.par_remove(files, &ids);
            }
        }
        spin.finish_with_message("Finished removing sequences!\n");
        self.print_output_info();
    }

    fn find_matching_ids(&self, files: &[PathBuf], re: &str) -> Vec<String> {
        let ids = IDs::new(files, self.input_fmt, self.datatype).id_unique();
        let re = Regex::new(re).expect("Failed parsing regex");
        let mut matching_ids = Vec::with_capacity(ids.len());
        ids.iter().for_each(|id| {
            if re.is_match(id) {
                matching_ids.push(id.to_string());
            }
        });
        matching_ids.shrink_to_fit();
        matching_ids
    }

    fn par_remove(&self, files: &[PathBuf], ids: &[String]) {
        files.par_iter().for_each(|file| {
            let (matrix, header) = self.remove_sequence(file, ids);
            if !matrix.is_empty() {
                self.write_output(&matrix, &header, file);
            }
        })
    }

    fn write_output(&self, matrix: &SeqMatrix, header: &Header, file: &Path) {
        let outpath = filenames::create_output_fname(self.outdir, file, self.output_fmt);
        let mut writer = SeqWriter::new(&outpath, matrix, header);
        writer
            .write_sequence(self.output_fmt)
            .expect("Failed writing output sequence");
    }

    fn remove_sequence(&self, fpath: &Path, ids: &[String]) -> (SeqMatrix, Header) {
        let (mut matrix, header) = SeqParser::new(fpath, self.datatype).parse(self.input_fmt);
        ids.iter()
            .for_each(|id| if matrix.remove(id).is_some() {});

        let fnl_header = if !matrix.is_empty() && header.ntax != matrix.len() {
            self.get_header(&matrix)
        } else {
            header
        };

        (matrix, fnl_header)
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

    fn print_output_info(&self) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Output dir", self.outdir.display());
        self.print_output_fmt(self.output_fmt);
    }
}

#[cfg(test)]
mod test {
    use super::*;

    macro_rules! input {
        ($remove: ident) => {
            let input_fmt = InputFmt::Fasta;
            let datatype = DataType::Dna;

            let opts = RemoveOpts::Regex(String::from("^abc"));
            let outdir = Path::new(".");
            let output_fmt = OutputFmt::Fasta;
            let $remove = Remove::new(&input_fmt, &datatype, outdir, &output_fmt, &opts);
        };
    }

    #[test]
    fn test_remove_seq() {
        let ids = vec![String::from("ABCD")];
        input!(remove);
        let file = Path::new("tests/files/simple.fas");
        let (_, header) = remove.remove_sequence(file, &ids);
        assert_eq!(header.ntax, 1);
    }

    #[test]
    fn test_remove_regex() {
        let re = String::from("(?i)^abc");
        input!(remove);
        let file = PathBuf::from("tests/files/simple.fas");
        let files = [file.clone()];
        let ids = remove.find_matching_ids(&files, &re);
        let (_, header) = remove.remove_sequence(&file, &ids);
        assert_eq!(header.ntax, 1);
    }
}
