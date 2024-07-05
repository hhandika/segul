//! Create unique IDs from sequence alignment files and map them to the alignment files.
use std::ffi::OsStr;
use std::fs::{self};
use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;

use colored::Colorize;
use indexmap::IndexSet;
use rayon::prelude::*;

use crate::helper::finder::IDs;
use crate::helper::sequence::SeqParser;
use crate::helper::types::{DataType, InputFmt};
use crate::helper::utils;
use crate::writer::text::IdWriter;

/// The `Id` struct is used to generate unique IDs
/// from sequence alignment files and map them to the alignment files.
pub struct SequenceID<'a> {
    /// The sequence alignment files.
    files: &'a [PathBuf],
    /// The input format of the sequence alignment files.
    pub input_fmt: &'a InputFmt,
    /// The data type of the sequence alignment files.
    pub datatype: &'a DataType,
    /// The output path for the unique IDs.
    pub output: &'a Path,
    /// The prefix for the output file (optional)
    /// If provided, the output file will be named as `prefix_id.txt`.
    /// If not provided, the output file will use the default name `id.txt`.
    /// or `prefix_map.csv` if mapping is enabled.
    pub prefix: Option<&'a str>,
}

impl<'a> SequenceID<'a> {
    pub fn new(
        files: &'a [PathBuf],
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        output: &'a Path,
        prefix: Option<&'a str>,
    ) -> Self {
        Self {
            files,
            input_fmt,
            datatype,
            output,
            prefix: prefix,
        }
    }

    /// Generate unique IDs from sequence alignment files.
    /// The unique IDs are written to a file as a text file.
    /// The file will be named as `id.txt` or `prefix_id.txt` if prefix is provided.
    /// Example:
    /// ```rust
    /// use std::path::{Path, PathBuf};
    /// use segul::core::sequence::id::SequenceID;
    /// use segul::helper::types::{DataType, InputFmt};
    /// use tempdir::TempDir;
    ///
    /// let alignment_2 = PathBuf::from("tests/files/concat/gene_2.nex");
    /// let alignment_1 = PathBuf::from("tests/files/concat/gene_1.nex");
    /// let files = vec![alignment_1, alignment_2];
    /// let output = TempDir::new("tempt").unwrap();
    /// let handle = SequenceID::new(&files, &InputFmt::Auto, &DataType::Dna, Path::new(output.path()), None);
    /// handle.get_unique();
    /// assert!(output.path().join("id.txt").exists());
    /// ```

    pub fn get_unique(&self) {
        fs::create_dir_all(self.output.parent().expect("Failed getting parent path"))
            .expect("Failed creating output dir");
        let spin = utils::set_spinner();
        spin.set_message("Indexing IDs..");
        let ids = self.get_unique_id(self.files);
        spin.finish_with_message("DONE!\n");
        let writer = IdWriter::new(self.output, &ids, self.prefix);
        writer.write_unique_id().expect("Failed writing results");
        self.print_output(ids.len());
    }

    pub fn map_id(&self) {
        let spin = utils::set_spinner();
        spin.set_message("Mapping IDs..");
        let ids = self.get_unique_id(self.files);
        let mapped_ids = self.par_map_id(self.files, &ids);
        let writer = IdWriter::new(self.output, &ids, self.prefix);
        writer
            .write_unique_id()
            .expect("Failed writing unique IDs to file");
        writer
            .write_mapped_id(&mapped_ids)
            .expect("Failed writing mapped ID to file");
        spin.finish_with_message("DONE!\n");
        self.print_output(ids.len());
    }

    fn get_unique_id(&self, files: &[PathBuf]) -> IndexSet<String> {
        let mut id = IDs::new(files, self.input_fmt, self.datatype).id_unique();
        id.sort();
        id
    }

    fn par_map_id(&self, files: &[PathBuf], ids: &IndexSet<String>) -> Vec<IdRecords> {
        let (sender, receiver) = channel();
        files.par_iter().for_each_with(sender, |s, file| {
            s.send(self.map_id_to_aln(file, ids))
                .expect("Error in mapping IDs");
        });
        let mut records: Vec<IdRecords> = receiver.iter().collect();
        records.par_sort_by(|a, b| alphanumeric_sort::compare_str(&a.name, &b.name));
        records
    }

    fn map_id_to_aln(&self, file: &Path, ids: &IndexSet<String>) -> IdRecords {
        let fstem = self.get_aln_name(file);
        let mut rec = IdRecords::new(fstem, ids.len());
        let (seq, _) = SeqParser::new(file, self.datatype).parse(self.input_fmt);
        ids.iter().for_each(|id| {
            let is_id_present = seq.contains_key(id);
            rec.records.push(is_id_present);
        });
        rec
    }

    fn get_aln_name(&self, file: &Path) -> String {
        file.file_stem()
            .and_then(OsStr::to_str)
            .expect("Failed getting file stem for mapping IDs")
            .to_string()
    }

    fn print_output(&self, ids: usize) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Total unique IDs", ids);
        log::info!("{:18}: {}", "Output dir", self.output.display());
    }
}

pub struct IdRecords {
    pub name: String,
    pub records: Vec<bool>,
}

impl IdRecords {
    fn new(name: String, size: usize) -> Self {
        Self {
            name,
            records: Vec::with_capacity(size),
        }
    }
}
