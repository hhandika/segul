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
    collections::{HashMap, HashSet},
    path::{Path, PathBuf},
    sync::{Mutex, RwLock},
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
    /// Include or not include skipped files.
    /// If true, skipped files will be written to the output.
    /// in the same format as the output.
    added_only: bool,
}

impl Default for SequenceAddition<'_> {
    fn default() -> Self {
        Self {
            input_files: &[],
            input_fmt: &InputFmt::Auto,
            datatype: &DataType::Dna,
            output: Path::new("output"),
            output_fmt: &OutputFmt::Fasta,
            added_only: false,
        }
    }
}

impl OutputPrint for SequenceAddition<'_> {}

impl<'a> SequenceAddition<'a> {
    pub fn new(
        input_files: &'a [PathBuf],
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        output: &'a Path,
        output_fmt: &'a OutputFmt,
        added_only: bool,
    ) -> Self {
        Self {
            input_files,
            input_fmt,
            datatype,
            output,
            output_fmt,
            added_only,
        }
    }

    pub fn add(&self, dest_file: &[PathBuf], dest_fmt: &InputFmt) {
        let spinner = utils::set_spinner();
        spinner.set_message("Adding sequences...");
        let counter = self.add_sequences(dest_file, dest_fmt);
        spinner.finish_with_message("Finished adding sequences.\n");
        let mut total_written = 0;
        if !self.added_only {
            let skipped_files = self.get_skipped_files(&counter, dest_file);
            total_written = self.write_skip_files(&skipped_files);
        }
        self.print_output_info(&counter, total_written);
    }

    fn add_sequences(&self, dest_file: &[PathBuf], dest_fmt: &InputFmt) -> SequenceCounter {
        let dest_collection = self.create_dest_library(dest_file);
        let counter = Mutex::new(SequenceCounter::new(
            self.input_files.len(),
            dest_file.len(),
        ));
        self.input_files.par_iter().for_each(|input| {
            let input_stem = self.get_file_stem(input);
            match dest_collection.get(&input_stem) {
                Some(dest_file) => {
                    let input_matrix = self.get_matrix(input, self.input_fmt);
                    let dest_matrix =
                        self.create_final_matrix(input_matrix, dest_file, dest_fmt, &counter);
                    match dest_matrix {
                        Some(matrix) => {
                            self.write_output(&matrix, dest_file);
                            counter
                                .lock()
                                .expect("Failed to lock counter.")
                                .add_file(&input_stem);
                        }
                        None => {
                            counter
                                .lock()
                                .expect("Failed to lock counter.")
                                .skip_file(&input_stem);
                        }
                    }
                }
                None => {
                    log::warn!("No destination file found for {}. Skipping...", &input_stem);
                    counter
                        .lock()
                        .expect("Failed to lock counter.")
                        .skip_file(&input_stem);
                }
            };
        });
        let mut counter = counter.into_inner().expect("Failed to get counter.");
        counter.calculate_mean();
        counter
    }

    fn get_skipped_files(&self, counter: &SequenceCounter, dest_files: &[PathBuf]) -> Vec<PathBuf> {
        dest_files
            .iter()
            .filter(|file| !counter.added_files.contains(&self.get_file_stem(file)))
            .cloned()
            .collect()
    }

    fn write_skip_files(&self, skipped_files: &[PathBuf]) -> usize {
        let counter = RwLock::new(0);
        skipped_files.par_iter().for_each(|file| {
            let final_matrix = self.get_matrix(file, self.input_fmt);
            self.write_output(&final_matrix, file);
            *counter.write().expect("Failed to write counter.") += 1;
        });
        let counter = *counter.read().expect("Failed to read counter.");
        counter
    }

    fn create_final_matrix(
        &self,
        input_matrix: SeqMatrix,
        dest_file: &Path,
        dest_fmt: &InputFmt,
        counter: &Mutex<SequenceCounter>,
    ) -> Option<SeqMatrix> {
        let mut dest_matrix = self.get_matrix(dest_file, dest_fmt);
        let mut added_count = 0;
        input_matrix.iter().for_each(|(name, sequence)| {
            if dest_matrix.contains_key(name) {
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
                added_count += 1;
            }
        });
        if added_count > 0 {
            self.clean_missing_data(&mut dest_matrix);
            Some(dest_matrix)
        } else {
            None
        }
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

    fn print_output_info(&self, counter: &SequenceCounter, total_written: usize) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Directory", self.output.display());
        self.print_output_fmt(self.output_fmt);

        log::info!("\n{}", "File Summary".yellow());
        log::info!("{:18}: {}", "Total input files", counter.total_input_files);
        log::info!(
            "{:18}: {}",
            "Total dest files",
            counter.total_destination_files
        );
        log::info!("{:18}: {}", "Files skipped", counter.skipped_file_counts);
        log::info!("{:18}: {}", "Files added", counter.total_file_added);
        log::info!(
            "{:18}: {}\n",
            "Files written",
            counter.total_file_added + total_written
        );
        log::info!("{}", "Sequences Summary".yellow());
        log::info!(
            "{:18}: {}",
            "Total sequences",
            counter.total_sequence_counts
        );
        log::info!("{:18}: {}", "Sequence skipped", counter.skipped_sequences);
        log::info!("{:18}: {}", "Sequence added", counter.total_sequence_added);
        log::info!("{:18}: {:.2}", "Mean length", counter.mean_length);
        if counter.total_sequence_added > 0 {
            log::info!("{:18}: {:.2}", "Mean added", counter.mean_added_sequences);
            log::info!(
                "{:18}: {:.2}",
                "Mean added length",
                counter.mean_added_length
            );
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
struct SequenceCounter {
    total_input_files: usize,
    total_destination_files: usize,
    total_file_added: usize,
    /// Store file stem
    added_files: HashSet<String>,
    /// Store file stem
    skipped_files: Vec<String>,
    skipped_file_counts: usize,
    skipped_sequences: usize,
    total_sequence_counts: usize,
    total_sequence_added: usize,
    mean_added_sequences: f64,
    mean_length: f64,
    total_length: usize,
    total_added_length: usize,
    mean_added_length: f64,
}

impl SequenceCounter {
    fn new(total_input_files: usize, total_destination_files: usize) -> Self {
        Self {
            total_input_files,
            total_destination_files,
            total_file_added: 0,
            added_files: HashSet::new(),
            skipped_files: Vec::new(),
            skipped_file_counts: 0,
            total_sequence_added: 0,
            skipped_sequences: 0,
            total_sequence_counts: 0,
            mean_added_sequences: 0.0,
            mean_length: 0.0,
            total_length: 0,
            total_added_length: 0,
            mean_added_length: 0.0,
        }
    }

    fn calculate_mean(&mut self) {
        if self.total_file_added > 0 {
            self.mean_added_sequences =
                self.total_sequence_added as f64 / self.total_sequence_counts as f64;
            self.mean_added_length =
                self.total_added_length as f64 / self.total_sequence_added as f64;
        }
        self.mean_length = self.total_length as f64 / self.total_sequence_counts as f64;
    }

    fn add_sequence(&mut self, sequence: &str) {
        self.total_sequence_counts += 1;
        self.total_sequence_added += 1;
        self.total_length += sequence.len();
        self.total_added_length += sequence.len();
    }

    fn skip_sequence(&mut self, sequence: &str) {
        self.total_sequence_counts += 1;
        self.skipped_sequences += 1;
        self.total_length += sequence.len();
    }

    fn skip_file(&mut self, file_stem: &str) {
        self.skipped_file_counts += 1;
        self.skipped_files.push(file_stem.to_string());
    }

    fn add_file(&mut self, file_stem: &str) {
        self.total_file_added += 1;
        self.added_files.insert(file_stem.to_string());
    }
}

#[cfg(test)]
mod tests {
    use tempdir::TempDir;

    use super::*;

    #[test]
    fn test_stat_counter() {
        let mut counter = SequenceCounter::new(2, 2);
        counter.add_sequence("ATCG");
        counter.add_sequence("ATCG");
        counter.skip_sequence("ATCG");
        counter.skip_sequence("ATCG");
        counter.add_file("file");
        counter.calculate_mean();
        assert_eq!(counter.total_input_files, 2);
        assert_eq!(counter.total_sequence_counts, 4);
        assert_eq!(counter.total_sequence_added, 2);
        assert_eq!(counter.skipped_sequences, 2);
        assert_eq!(counter.total_file_added, 1);
        assert_eq!(counter.skipped_file_counts, 0);
        assert_eq!(counter.mean_added_sequences, 0.5);
        assert_eq!(counter.mean_length, 4.0);
        assert_eq!(counter.total_length, 16);
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
            &OutputFmt::Fasta,
            false,
        );
        let counter = addition.add_sequences(&dest_files, &InputFmt::Auto);
        assert_eq!(counter.total_input_files, 2);
        assert_eq!(counter.total_sequence_added, 2);
        assert_eq!(counter.skipped_sequences, 3);
        let output_files = output.path().read_dir().unwrap();
        assert_eq!(output_files.count(), 2);
        output.close().unwrap();
    }

    #[test]
    fn test_finding_skipped_files() {
        let mut counter = SequenceCounter::new(2, 3);
        counter.add_file("gene_1");
        let dest_files = vec![
            PathBuf::from("tests/files/alignments/gene_1.nex"),
            PathBuf::from("tests/files/alignments/gene_2.nex"),
            PathBuf::from("tests/files/alignments/gene_3.nex"),
        ];
        let addition = SequenceAddition::default();
        let skipped_files = addition.get_skipped_files(&counter, &dest_files);
        assert_eq!(skipped_files.len(), 2);
    }
}
