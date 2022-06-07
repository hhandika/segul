use std::path::{Path, PathBuf};

use ansi_term::Colour::Yellow;
use rayon::prelude::*;
use regex::Regex;

use crate::handler::OutputPrint;
use crate::helper::filenames;
use crate::helper::sequence::SeqParser;
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix};
use crate::helper::utils;
use crate::writer::sequences::SeqWriter;
use indexmap::IndexMap;

impl OutputPrint for Rename<'_> {}

macro_rules! process_files {
    ($self: ident, $files: ident, $func: ident, $outdir: ident, $output_fmt: ident, $($input: tt)*) => {
        $files.par_iter().for_each(|file| {
            let (matrix, header) = $self.$func(file, $($input)*);
            let outpath = filenames::create_output_fname($outdir, file, $output_fmt);
            $self.write_output(&matrix, &header, &outpath, $output_fmt);
        });
    };
}

pub enum RenameOpts {
    RnId(Vec<(String, String)>),   // Rename ID using tabulated file
    RmStr(String),                 // Remove characters in seq id using string input
    RmRegex(String, bool),         // Similar to RmStr but using regex as input
    RpStr(String, String),         // Replace characters in seq id using string input
    RpRegex(String, String, bool), // Similar to ReplaceStr but using regex as input
}

pub struct Rename<'a> {
    input_fmt: &'a InputFmt,
    datatype: &'a DataType,
    opts: &'a RenameOpts,
}

impl<'a> Rename<'a> {
    pub fn new(input_fmt: &'a InputFmt, datatype: &'a DataType, opts: &'a RenameOpts) -> Self {
        Self {
            input_fmt,
            datatype,
            opts,
        }
    }

    pub fn dry_run(&self) {
        let names = match self.opts {
            RenameOpts::RnId(names) => names,
            _ => unimplemented!(),
        };
        log::info!("{:18}: {}", "New ID count", names.len());
        log::info!("{:18}: Dry run\n", "Status");
        log::info!("{}", Yellow.paint("Results"));
        names
            .iter()
            .for_each(|(origin, destination)| log::info!("{} -> {}", origin, destination));
        println!();
    }

    pub fn rename(&self, files: &[PathBuf], outdir: &Path, output_fmt: &OutputFmt) {
        let spin = utils::set_spinner();
        spin.set_message("Batch renaming dna sequence IDs...");
        match self.opts {
            RenameOpts::RnId(names) => {
                process_files!(self, files, replace_id, outdir, output_fmt, names);
            }
            RenameOpts::RmStr(input_str) => {
                process_files!(self, files, replace_str, outdir, output_fmt, input_str, "");
            }
            RenameOpts::RmRegex(input_re, is_all) => {
                process_files!(self, files, replace_re, outdir, output_fmt, input_re, "", is_all);
            }
            RenameOpts::RpStr(input_str, output_str) => {
                process_files!(
                    self,
                    files,
                    replace_str,
                    outdir,
                    output_fmt,
                    input_str,
                    output_str
                );
            }
            RenameOpts::RpRegex(input_re, output_re, is_all) => {
                process_files!(
                    self, files, replace_re, outdir, output_fmt, input_re, output_re, is_all
                );
            }
        }
        spin.finish_with_message("Finished batch renaming dna sequence IDs!\n");
        self.print_output_info(outdir, output_fmt);
    }

    fn write_output(
        &self,
        matrix: &SeqMatrix,
        header: &Header,
        outpath: &Path,
        output_fmt: &OutputFmt,
    ) {
        let mut writer = SeqWriter::new(&outpath, &matrix, &header);
        writer
            .write_sequence(output_fmt)
            .expect("Failed writing output sequence");
    }

    fn replace_id(&self, file: &Path, names: &[(String, String)]) -> (SeqMatrix, Header) {
        let (mut matrix, header) = SeqParser::new(file, self.datatype).get(self.input_fmt);
        let original_size = matrix.len();
        names.iter().for_each(|(origin, destination)| {
            let values = matrix.remove(origin);
            if let Some(value) = values {
                matrix.insert(destination.to_string(), value);
            }
        });

        assert_eq!(
            original_size,
            matrix.len(),
            "Failed renaming files. New ID counts does not match original ID counts. \
         Original ID counts: {}. New ID counts: {}",
            original_size,
            matrix.len()
        );
        (matrix, header)
    }

    fn replace_str(&self, file: &Path, str_from: &str, str_to: &str) -> (SeqMatrix, Header) {
        let (matrix, header) = SeqParser::new(file, self.datatype).get(self.input_fmt);
        let mut new_matrix = IndexMap::with_capacity(matrix.len());
        matrix.iter().for_each(|(id, seq)| {
            if id.contains(str_from) {
                let new_id = id.replace(str_from, str_to);
                new_matrix.insert(new_id, seq.to_string());
            } else {
                new_matrix.insert(id.to_string(), seq.to_string());
            }
        });

        (new_matrix, header)
    }

    fn replace_re(
        &self,
        file: &Path,
        re_from: &str,
        str_to: &str,
        all: &bool,
    ) -> (SeqMatrix, Header) {
        let (matrix, header) = SeqParser::new(file, self.datatype).get(self.input_fmt);
        let mut new_matrix = IndexMap::with_capacity(matrix.len());
        matrix.iter().for_each(|(id, seq)| {
            let re = Regex::new(re_from).expect("Failed parsing regex");
            let new_id = if *all {
                re.replace_all(id, str_to)
            } else {
                re.replace(id, str_to)
            };
            new_matrix.insert(new_id.to_string(), seq.to_string());
        });

        (new_matrix, header)
    }

    fn print_output_info(&self, output: &Path, output_fmt: &OutputFmt) {
        log::info!("{}", Yellow.paint("Output"));
        log::info!("{:18}: {}", "Output dir", output.display());
        self.print_output_fmt(output_fmt);
    }
}

#[cfg(test)]
mod test {
    use super::*;

    macro_rules! input {
        ($rename: ident, $file: ident) => {
            let input_fmt = InputFmt::Fasta;
            let datatype = DataType::Dna;
            let $file = Path::new("tests/files/simple.fas");
            let opts = RenameOpts::RmStr(String::from("test"));
            let $rename = Rename::new(&input_fmt, &datatype, &opts);
        };
    }

    #[test]
    fn test_rename_id() {
        input!(rename, file);
        let names = [(String::from("ABCD"), String::from("WXYZ"))];
        let (seq, _) = rename.replace_id(&file, &names);
        assert_eq!(seq.len(), 2);
        assert_eq!(seq.get("WXYZ"), Some(&String::from("AGTATG")));
        assert_eq!(seq.get("ABCD"), None);
    }

    #[test]
    fn test_rename_rm_str() {
        input!(rename, file);
        let (seq, _) = rename.replace_str(&file, "BC", "");
        assert_eq!(seq.get("AD"), Some(&String::from("AGTATG")));
    }

    #[test]
    fn test_rename_rm_re() {
        input!(rename, file);
        let (seq, _) = rename.replace_re(&file, "^A", "", &false);
        let (seq2, _) = rename.replace_re(&file, "[^ABC]", "", &false);
        assert_eq!(seq.get("BCD"), Some(&String::from("AGTATG")));
        assert_eq!(seq2.get("ABC"), Some(&String::from("AGTATG")));
    }
}
