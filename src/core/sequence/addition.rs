//! Add sequence to a file.
//!
//! Support adding many sequences at once.
//! The destination file and the sequences can
//! be in different formats, but the file name,
//! excluding the extension, (file stem) must be the same.
//! Input sequences can also be an alignment files.
//! However, SEGUL will output unaligned sequences
//! by excluding the gaps and missing data
//! from the the resulting file.

use std::{
    collections::HashMap,
    path::{Path, PathBuf},
    sync::Mutex,
};

use colored::Colorize;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::{
    core::OutputPrint,
    helper::{
        files,
        sequence::SeqParser,
        types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix},
        utils,
    },
    writer::sequences::SeqWriter,
};

pub struct SequenceAddition<'a> {
    /// Input sequence files to be added.
    input_files: &'a [PathBuf],
    /// Input sequence file format.
    input_fmt: &'a InputFmt,
    /// Data type of the sequences.
    datatype: &'a DataType,
    /// Destination file to add the sequences.
    output: &'a Path,
    /// Output file format.
    output_fmt: &'a OutputFmt,
}

impl OutputPrint for SequenceAddition<'_> {}

impl<'a> SequenceAddition<'a> {
    pub fn new(
        input_files: &'a [PathBuf],
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        output: &'a Path,
        output_fmt: &'a OutputFmt,
    ) -> Self {
        Self {
            input_files,
            input_fmt,
            datatype,
            output,
            output_fmt,
        }
    }

    pub fn add(&self, dest_file: &[PathBuf], dest_fmt: &InputFmt) {
        let spinner = utils::set_spinner();
        spinner.set_message("Adding sequences...");
        let counter = self.add_sequences(dest_file, dest_fmt);
        spinner.finish_with_message("Finished adding sequences.");
        self.print_output_info(&counter);
    }

    fn add_sequences(&self, dest_file: &[PathBuf], dest_fmt: &InputFmt) -> SequenceCounter {
        let dest_collection = self.create_dest_library(dest_file);
        let counter = Mutex::new(SequenceCounter::new(self.input_files.len()));
        self.input_files.par_iter().for_each(|input| {
            let input_stem = self.get_file_stem(input);
            let input_matrix = self.get_matrix(input, self.input_fmt);
            let dest_file = dest_collection
                .get(&input_stem)
                .expect("Failed to get destination file.");
            let mut dest_matrix = self.get_matrix(dest_file, dest_fmt);
            input_matrix.iter().for_each(|(name, sequence)| {
                if dest_matrix.contains_key(name) {
                    log::warn!(
                        "Sequence {} already exists in the {} file. Skipping...",
                        name,
                        dest_file.to_string_lossy()
                    );
                    counter
                        .lock()
                        .expect("Failed to lock counter.")
                        .skip_sequence(sequence);
                } else {
                    dest_matrix.insert(name.to_string(), sequence.to_string());
                    counter
                        .lock()
                        .expect("Failed to lock counter.")
                        .add_sequence(sequence);
                }
            });
            self.clean_missing_data(&mut dest_matrix);
            self.write_output(&dest_matrix, dest_file);
        });
        let mut counter = counter.into_inner().expect("Failed to get counter.");
        counter.calculate_mean();
        counter
    }

    fn create_dest_library(&self, dest_file: &[PathBuf]) -> HashMap<String, PathBuf> {
        let dest_collection: Mutex<HashMap<String, PathBuf>> = Mutex::new(HashMap::new());
        dest_file.par_iter().for_each(|file| {
            let file_stem = self.get_file_stem(file);
            let mut dest_collection = dest_collection
                .lock()
                .expect("Failed to lock dest_collection.");
            dest_collection.insert(file_stem, file.clone());
        });
        dest_collection
            .into_inner()
            .expect("Failed to get dest_collection.")
    }

    fn get_file_stem(&self, file: &Path) -> String {
        file.file_stem()
            .expect("Failed to get file stem.")
            .to_string_lossy()
            .to_string()
    }

    fn get_matrix(&self, input: &Path, input_fmt: &InputFmt) -> SeqMatrix {
        let (seq, _) = SeqParser::new(input, self.datatype).parse(input_fmt);
        seq
    }

    fn clean_missing_data(&self, matrix: &mut SeqMatrix) {
        matrix.values_mut().for_each(|seq| {
            *seq = seq.replace(['?', '-'], "");
        });
    }

    fn write_output(&self, final_matrix: &SeqMatrix, file: &Path) {
        let output_path = files::create_output_fname(self.output, file, self.output_fmt);
        let mut header: Header = Header::new();
        header.from_seq_matrix(final_matrix, false);
        let mut writer = SeqWriter::new(&output_path, final_matrix, &header);
        writer
            .write_sequence(self.output_fmt)
            .expect("Failed to write sequences.");
    }

    fn print_output_info(&self, counter: &SequenceCounter) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Directory", self.output.display());
        log::info!("{:18}: {}", "Total files", counter.total_input_files);
        log::info!("{:18}: {}", "Total sequences", counter.total_sequence_added);
        log::info!("{:18}: {}", "Total added", counter.total_sequence_added);
        log::info!("{:18}: {}", "Total skipped", counter.skip_counter);
        log::info!("{:18}: {:.2}", "Mean sequences", counter.mean_sequence);
        log::info!("{:18}: {:.2}", "Mean length", counter.mean_length);
        log::info!(
            "{:18}: {:.2}",
            "Mean added length",
            counter.mean_added_length
        );
    }
}

#[derive(Debug, Serialize, Deserialize)]
struct SequenceCounter {
    total_input_files: usize,
    total_sequence_added: usize,
    skip_counter: usize,
    total_added: usize,
    total_sequence_count: usize,
    mean_sequence: f64,
    mean_length: f64,
    total_length: usize,
    total_added_length: usize,
    mean_added_length: f64,
}

impl SequenceCounter {
    fn new(total_input_files: usize) -> Self {
        Self {
            total_input_files,
            total_sequence_added: 0,
            skip_counter: 0,
            total_added: 0,
            total_sequence_count: 0,
            mean_sequence: 0.0,
            mean_length: 0.0,
            total_length: 0,
            total_added_length: 0,
            mean_added_length: 0.0,
        }
    }

    fn add_sequence(&mut self, sequence: &str) {
        self.total_added += 1;
        self.total_sequence_count += 1;
        self.total_sequence_added += 1;
        self.total_length += sequence.len();
        self.total_added_length += sequence.len();
    }

    fn skip_sequence(&mut self, sequence: &str) {
        self.total_sequence_count += 1;
        self.skip_counter += 1;
        self.total_length += sequence.len();
    }

    fn calculate_mean(&mut self) {
        self.mean_sequence = self.total_sequence_added as f64 / self.total_input_files as f64;
        self.mean_length = self.total_length as f64 / self.total_sequence_count as f64;
        self.mean_added_length = self.total_added_length as f64 / self.total_sequence_added as f64;
    }
}

#[cfg(test)]
mod tests {
    use tempdir::TempDir;

    use super::*;

    #[test]
    fn test_stat_counter() {
        let mut counter = SequenceCounter::new(2);
        counter.add_sequence("ATCG");
        counter.add_sequence("ATCG");
        counter.skip_sequence("ATCG");
        counter.calculate_mean();
        assert_eq!(counter.total_sequence_count, 3);
        assert_eq!(counter.total_sequence_added, 2);
        assert_eq!(counter.skip_counter, 1);
        assert_eq!(counter.mean_sequence, 1.0);
        assert_eq!(counter.mean_length, 4.0);
        assert_eq!(counter.total_length, 12);
        assert_eq!(counter.total_added_length, 8);
        assert_eq!(counter.mean_added_length, 4.0);
    }

    #[test]
    fn test_sequence_addition() {
        let input_files = vec![
            PathBuf::from("tests/files/gappy/gene_1.nex"),
            PathBuf::from("tests/files/gappy/gene_2.nex"),
        ];
        let dest_files = vec![
            PathBuf::from("tests/files/alignments/gene_1.nex"),
            PathBuf::from("tests/files/alignments/gene_2.nex"),
        ];
        let output = TempDir::new("temp").unwrap();
        let addition = SequenceAddition::new(
            &input_files,
            &InputFmt::Auto,
            &DataType::Dna,
            output.path(),
            &OutputFmt::Nexus,
        );
        addition.add(&dest_files, &InputFmt::Auto);
        let counter = addition.add_sequences(&dest_files, &InputFmt::Auto);
        assert_eq!(counter.total_input_files, 2);
        assert_eq!(counter.total_sequence_added, 2);
        assert_eq!(counter.skip_counter, 3);
        let output_files = output.path().read_dir().unwrap();
        assert_eq!(output_files.count(), 2);
        output.close().unwrap();
    }
}
