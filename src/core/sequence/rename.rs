use std::path::{Path, PathBuf};

use colored::Colorize;
use indexmap::IndexMap;
use indexmap::IndexSet;
use rayon::prelude::*;
use regex::Regex;

use crate::core::OutputPrint;
use crate::helper::files;
use crate::helper::finder::IDs;
use crate::helper::sequence::SeqParser;
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix};
use crate::helper::utils;
use crate::writer::sequences::SeqWriter;

impl OutputPrint for SequenceRenaming<'_> {}

macro_rules! process_files {
    ($self: ident, $files: ident, $func: ident, $($input: tt)*) => {
        $files.par_iter().for_each(|file| {
            let (matrix, header) = $self.$func(file, $($input)*);
            $self.write_output(&matrix, &header, file);
        });
    };
}

macro_rules! rm_id {
    ($new_ids: ident, $ids: ident) => {
        $new_ids.iter().for_each(|(old, _)| {
            $ids.shift_remove(old);
        });
    };
}

pub enum SeqRenamingParameters {
    /// Rename ID using tabulated file
    RnId(Vec<(String, String)>),
    /// Remove characters in seq id using string input
    RmStr(String),
    /// Similar to RmStr but using regex as input
    RmRegex(String, bool),
    /// Replace characters in seq id using string input
    RpStr(String, String),
    /// Similar to ReplaceStr but using regex as input
    RpRegex(String, String, bool),

    None,
}

pub struct SequenceRenamingDry<'a> {
    input_fmt: &'a InputFmt,
    datatype: &'a DataType,
    opts: &'a SeqRenamingParameters,
}

impl<'a> SequenceRenamingDry<'a> {
    pub fn new(
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        opts: &'a SeqRenamingParameters,
    ) -> Self {
        Self {
            input_fmt,
            datatype,
            opts,
        }
    }

    pub fn dry_run(&self, files: &[PathBuf]) {
        let spin = utils::set_spinner();
        spin.set_message("Processing dna sequence IDs (DRY-RUN)...");
        let mut ids = IDs::new(files, self.input_fmt, self.datatype).id_unique();
        let new_ids = match self.opts {
            SeqRenamingParameters::RnId(names) => self.replace_id(&mut ids, names),
            SeqRenamingParameters::RmStr(input_str) => self.replace_str(&mut ids, input_str, ""),
            SeqRenamingParameters::RmRegex(input_re, is_all) => {
                self.replace_re(&mut ids, input_re, "", is_all)
            }
            SeqRenamingParameters::RpStr(from, to) => self.replace_str(&mut ids, from, to),
            SeqRenamingParameters::RpRegex(from, to, is_all) => {
                self.replace_re(&mut ids, from, to, is_all)
            }
            SeqRenamingParameters::None => unreachable!("Missing rename parameters"),
        };
        spin.finish_with_message("Finished processing (DRY-RUN)!\n");

        // Print results
        log::info!("{}", "Results".yellow());
        log::info!("{:18}: {}", "Renamed ID counts", new_ids.len());
        new_ids.iter().for_each(|(old, new)| {
            log::info!("{:18}: {} {} {}", "[Rename]", old, "->".green(), new);
        });

        // Print remaining unchanged ids
        log::info!("");
        log::info!("{:18}: {}", "Unchanged counts", ids.len());
        if !ids.is_empty() {
            ids.iter().for_each(|id| {
                log::info!("{:18}: {}", "[Unchanged]", id);
            });
        }
        println!();
    }

    fn replace_id(
        &self,
        ids: &mut IndexSet<String>,
        names: &[(String, String)],
    ) -> Vec<(String, String)> {
        let mut new_ids: Vec<(String, String)> = Vec::new();
        names.iter().for_each(|(old, new)| {
            let is_id = ids.shift_remove(old);
            if is_id {
                new_ids.push((old.to_string(), new.to_string()));
            }
        });
        new_ids
    }

    fn replace_str(
        &self,
        ids: &mut IndexSet<String>,
        from: &str,
        to: &str,
    ) -> Vec<(String, String)> {
        let mut new_ids: Vec<(String, String)> = Vec::new();
        ids.iter().for_each(|id| {
            if id.contains(from) {
                let new_id = id.replace(from, to);
                new_ids.push((id.to_string(), new_id));
            }
        });

        rm_id!(new_ids, ids);

        new_ids
    }

    fn replace_re(
        &self,
        ids: &mut IndexSet<String>,
        from: &str,
        to: &str,
        all: &bool,
    ) -> Vec<(String, String)> {
        let mut new_ids: Vec<(String, String)> = Vec::new();
        ids.iter().for_each(|id| {
            let re = Regex::new(from).expect("Failed parsing regex");
            let new_id = if *all {
                re.replace_all(id, to)
            } else {
                re.replace(id, to)
            };
            let changed_id = id.to_string();
            if new_id != changed_id {
                new_ids.push((changed_id, new_id.to_string()));
            }
        });

        rm_id!(new_ids, ids);

        new_ids
    }
}

pub struct SequenceRenaming<'a> {
    input_fmt: &'a InputFmt,
    datatype: &'a DataType,
    output_dir: &'a Path,
    output_fmt: &'a OutputFmt,
    opts: &'a SeqRenamingParameters,
}

impl<'a> SequenceRenaming<'a> {
    pub fn new(
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        output_dir: &'a Path,
        output_fmt: &'a OutputFmt,
        opts: &'a SeqRenamingParameters,
    ) -> Self {
        Self {
            input_fmt,
            datatype,
            output_dir,
            output_fmt,
            opts,
        }
    }

    pub fn rename(&self, files: &[PathBuf]) {
        let spin = utils::set_spinner();
        spin.set_message("Batch renaming dna sequence IDs...");
        match self.opts {
            SeqRenamingParameters::RnId(names) => {
                process_files!(self, files, replace_id, names);
            }
            SeqRenamingParameters::RmStr(input_str) => {
                process_files!(self, files, replace_str, input_str, "");
            }
            SeqRenamingParameters::RmRegex(input_re, is_all) => {
                process_files!(self, files, replace_re, input_re, "", is_all);
            }
            SeqRenamingParameters::RpStr(from, to) => {
                process_files!(self, files, replace_str, from, to);
            }
            SeqRenamingParameters::RpRegex(from, to, is_all) => {
                process_files!(self, files, replace_re, from, to, is_all);
            }
            SeqRenamingParameters::None => unreachable!("Missing rename parameters"),
        }
        spin.finish_with_message("Finished batch renaming dna sequence IDs!\n");
        self.print_output_info();
    }

    fn write_output(&self, matrix: &SeqMatrix, header: &Header, file: &Path) {
        let output_path = files::create_output_fname(self.output_dir, file, self.output_fmt);
        let mut writer = SeqWriter::new(&output_path, matrix, header);
        writer
            .write_sequence(self.output_fmt)
            .expect("Failed writing output sequence");
    }

    fn replace_id(&self, file: &Path, names: &[(String, String)]) -> (SeqMatrix, Header) {
        let (mut matrix, mut header) = SeqParser::new(file, self.datatype).parse(self.input_fmt);
        let original_size = matrix.len();
        names.iter().for_each(|(origin, destination)| {
            if let Some(value) = matrix.shift_remove(origin) {
                // We check the sequence length if the ID is not unique
                // to avoid overwriting the sequence with a shorter one.
                if let Some(existing_value) = matrix.get(destination) {
                    log::warn!(
                        "\n{}! ID {} already exists. \
                    Keeping the longest sequences.",
                        "Warning".yellow(),
                        destination
                    );
                    if existing_value.len() > value.len() {
                        matrix.insert(origin.to_string(), value);
                    }
                } else {
                    matrix.insert(destination.to_string(), value);
                }
            }
        });

        if original_size != matrix.len() {
            log::warn!(
                "Duplicate IDs found! Original ID count: {}. \
            New ID count: {}. Total duplicates: {}\n",
                original_size,
                matrix.len(),
                original_size - matrix.len()
            );
            header.update(&matrix);
        }
        (matrix, header)
    }

    fn replace_str(&self, file: &Path, from: &str, to: &str) -> (SeqMatrix, Header) {
        let (matrix, header) = SeqParser::new(file, self.datatype).parse(self.input_fmt);
        let mut new_matrix = IndexMap::with_capacity(matrix.len());
        matrix.iter().for_each(|(id, seq)| {
            if id.contains(from) {
                let new_id = id.replace(from, to);
                new_matrix.insert(new_id, seq.to_string());
            } else {
                new_matrix.insert(id.to_string(), seq.to_string());
            }
        });

        (new_matrix, header)
    }

    fn replace_re(&self, file: &Path, from: &str, to: &str, all: &bool) -> (SeqMatrix, Header) {
        let (matrix, header) = SeqParser::new(file, self.datatype).parse(self.input_fmt);
        let mut new_matrix = IndexMap::with_capacity(matrix.len());
        matrix.iter().for_each(|(id, seq)| {
            let re = Regex::new(from).expect("Failed parsing regex");
            let new_id = if *all {
                re.replace_all(id, to)
            } else {
                re.replace(id, to)
            };
            new_matrix.insert(new_id.to_string(), seq.to_string());
        });

        (new_matrix, header)
    }

    fn print_output_info(&self) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Output dir", self.output_dir.display());
        self.print_output_fmt(self.output_fmt);
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
            let opts = SeqRenamingParameters::RmStr(String::from("test"));
            let output_dir = Path::new(".");
            let output_fmt = OutputFmt::Fasta;
            let $rename =
                SequenceRenaming::new(&input_fmt, &datatype, output_dir, &output_fmt, &opts);
        };
    }

    #[test]
    fn test_rename_id() {
        input!(rename, file);
        let names = [(String::from("ABCD"), String::from("WXYZ"))];
        let (seq, _) = rename.replace_id(file, &names);
        assert_eq!(seq.len(), 2);
        assert_eq!(seq.get("WXYZ"), Some(&String::from("AGTATG")));
        assert_eq!(seq.get("ABCD"), None);
    }

    #[test]
    fn test_rename_rm_str() {
        input!(rename, file);
        let (seq, _) = rename.replace_str(file, "BC", "");
        assert_eq!(seq.get("AD"), Some(&String::from("AGTATG")));
    }

    #[test]
    fn test_rename_rm_re() {
        input!(rename, file);
        let (seq, _) = rename.replace_re(file, "^A", "", &false);
        let (seq2, _) = rename.replace_re(file, "[^ABC]", "", &false);
        assert_eq!(seq.get("BCD"), Some(&String::from("AGTATG")));
        assert_eq!(seq2.get("ABC"), Some(&String::from("AGTATG")));
    }
}
