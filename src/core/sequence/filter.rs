//! Filter sequence based on selected criteria.
//!
//! Unlike the `alignment filter` that works on the alignment level,
//! the `sequence filter` will filter sequences within the alignment.
//! The filtering criteria are:
//! - Total percentage of gaps in a sequence. The sequence will be removed
//! if the percentage of gaps is higher than the threshold.
//! - Minimum sequence length. Filter sequence that has an equal or more sequence length
//! than the specified number.

use std::{
    path::{Path, PathBuf},
    sync::atomic::{AtomicUsize, Ordering},
};

use colored::Colorize;
use rayon::prelude::*;

use crate::{
    helper::{
        concat::ConcatParams,
        files,
        sequence::SeqParser,
        types::{DataType, Header, InputFmt, OutputFmt, PartitionFmt, SeqMatrix},
        utils,
    },
    writer::sequences::SeqWriter,
};

/// Available sequence filtering options.
pub enum SequenceFilteringOptions {
    /// Filtered sequences based the percentage of gaps in a sequence.
    /// Remove sequences that have a percentage
    /// of gaps higher than the threshold.
    /// The threshold is calculated as the number
    /// of gaps times the sequence length.
    /// For example, if the threshold is 0.5
    /// and the sequence length is 10,
    /// the sequence will be removed
    /// if it has more than five gaps.
    /// If the threshold calculation results in a float,
    /// the float will be floored.
    /// For example, 5.1, 5.5, or 5.9, will be 5.
    GapThreshold(f64),
    /// Filter sequences based on the minimum sequence length.
    /// Remove sequences that have a length less than the specified number.
    MinSequenceLength(usize),
}

pub struct SequenceFiltering<'a> {
    /// List of input files.
    files: &'a [PathBuf],
    /// Input format.
    input_fmt: &'a InputFmt,
    /// Data type.
    datatype: &'a DataType,
    /// Output directory.
    output: &'a Path,
    /// Output format.
    output_fmt: &'a OutputFmt,
    /// Choice of filtering options.
    params: &'a SequenceFilteringOptions,
    /// Concatenation parameters.
    concat: Option<ConcatParams>,
}

impl<'a> SequenceFiltering<'a> {
    pub fn new(
        files: &'a [PathBuf],
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        output: &'a Path,
        output_fmt: &'a OutputFmt,
        params: &'a SequenceFilteringOptions,
    ) -> Self {
        Self {
            files,
            input_fmt,
            datatype,
            output,
            output_fmt,
            params,
            concat: None,
        }
    }

    /// Filter sequences based on the selected criteria.
    /// The resulting files are alignments with filtered sequences.
    /// Example:
    /// ```
    /// use std::path::{Path, PathBuf};
    /// use tempdir::TempDir;
    /// use segul::helper::types::{DataType, InputFmt, OutputFmt, PartitionFmt};
    /// use segul::core::sequence::filter::{SequenceFiltering, SequenceFilteringOptions};
    /// use segul::helper::finder::SeqFileFinder;
    ///
    ///
    /// let input_fmt = InputFmt::Nexus;
    /// let datatype = DataType::Dna;
    /// let input_dir = Path::new("tests/files/concat");
    /// let files = SeqFileFinder::new(Path::new(input_dir)).find(&input_fmt);
    /// // Replace the temp directory with your own directory.
    /// let output = TempDir::new("tempt").unwrap();
    /// let output_fmt = OutputFmt::Nexus;
    /// let params = SequenceFilteringOptions::MinSequenceLength(7);
    /// let handle = SequenceFiltering::new(&files, &input_fmt, &datatype, &output.path(), &output_fmt,  &params);
    /// handle.filter();
    pub fn filter(&self) {
        let spinner = utils::set_spinner();
        spinner.set_message("Filtering sequences...");
        let filtered_aln = match self.params {
            SequenceFilteringOptions::GapThreshold(threshold) => {
                self.filter_gappy_sequences(threshold)
            }
            SequenceFilteringOptions::MinSequenceLength(min_length) => {
                self.filter_sequences_by_length(min_length)
            }
        };
        spinner.finish_with_message("Finished filtering sequences!\n");
        if filtered_aln == 0 {
            log::warn!("No matching sequences were found!!");
        } else {
            self.print_output_info(filtered_aln);
        }
    }

    /// Setter for the concatenation parameters.
    pub fn set_concat(&mut self, partition_fmt: &PartitionFmt, prefix: &Path) {
        self.concat = Some(ConcatParams::new(self.output_fmt, partition_fmt, prefix));
    }

    fn filter_gappy_sequences(&self, threshold: &f64) -> usize {
        let counter = AtomicUsize::new(0);
        self.files.par_iter().for_each(|file| {
            let (mut matrix, mut header) = self.get_alignment(file);
            self.remove_gappy_sequences(&mut matrix, &mut header, threshold);
            if header.ntax > 0 {
                self.write_sequence(file, &matrix, &header);
                counter.fetch_add(1, Ordering::Relaxed);
            }
        });

        counter.into_inner()
    }

    fn write_sequence(&self, file: &Path, matrix: &SeqMatrix, header: &Header) {
        let output_path = files::create_output_fname(&self.output, file, &self.output_fmt);
        let mut seq_writer = SeqWriter::new(&output_path, matrix, header);
        seq_writer
            .write_sequence(self.output_fmt)
            .expect("Failed writing the output file");
    }

    fn print_output_info(&self, counter: usize) {
        log::info!("{}", "Output".yellow());
        log::info!("Directory: {}", self.output.display());
        log::info!("Format: {:?}", self.output_fmt);
        log::info!("Total files: {}", counter);
    }

    fn filter_sequences_by_length(&self, length: &usize) -> usize {
        let counter = AtomicUsize::new(0);
        self.files.par_iter().for_each(|file| {
            let (mut matrix, mut header) = self.get_alignment(file);
            matrix.retain(|_, seq| seq.len() >= *length);
            header.ntax = matrix.len();
            if header.ntax > 0 {
                self.write_sequence(file, &matrix, &header);
                counter.fetch_add(1, Ordering::Relaxed);
            }
        });

        counter.into_inner()
    }

    fn get_alignment(&self, file: &Path) -> (SeqMatrix, Header) {
        let alignment = SeqParser::new(file, &self.datatype);
        alignment.get_alignment(self.input_fmt)
    }

    fn remove_gappy_sequences(&self, matrix: &mut SeqMatrix, header: &mut Header, threshold: &f64) {
        let alignment_len = header.nchar;
        let max_gaps = self.count_max_gaps(alignment_len, *threshold);
        let mut to_remove = Vec::new();
        // Find sequences with gaps higher than the threshold.
        matrix.iter().for_each(|(name, seq)| {
            let gaps = seq.bytes().filter(|&c| c == b'-' || c == b'?').count();
            if gaps <= max_gaps {
                to_remove.push(String::from(name));
            }
        });

        to_remove.iter().for_each(|name| {
            matrix.retain(|n, _| n == name);
        });

        header.ntax = matrix.len();
    }

    fn count_max_gaps(&self, length: usize, threshold: f64) -> usize {
        (length as f64 * threshold).floor() as usize
    }
}

#[cfg(test)]
mod tests {
    use crate::helper::finder::SeqFileFinder;

    use super::*;
    use tempdir::TempDir;

    macro_rules! setup {
        ($input_dir: ident, $handle: ident, $params: ident, $output: ident) => {
            let files = SeqFileFinder::new($input_dir).find(&InputFmt::Nexus);
            let input_fmt = InputFmt::Nexus;
            let datatype = DataType::Dna;
            let $output = TempDir::new("output").unwrap();
            let $handle = SequenceFiltering::new(
                &files,
                &input_fmt,
                &datatype,
                $output.path(),
                &OutputFmt::Nexus,
                &$params,
            );
        };
    }

    #[test]
    fn test_filter_sequences_by_length() {
        let dir = Path::new("tests/files/concat");
        let params = SequenceFilteringOptions::MinSequenceLength(7);
        setup!(dir, handle, params, output);
        handle.filter();
        let output_files = SeqFileFinder::new(output.path()).find(&InputFmt::Nexus);
        assert_eq!(output_files.len(), 1);
    }

    #[test]
    fn test_gap_calculation() {
        let input_dir = Path::new("tests/files/gappy");
        let params = SequenceFilteringOptions::MinSequenceLength(7);
        setup!(input_dir, handle, params, output);
        let threshold = 0.5;
        let test_file = input_dir.join("gene_1.nex");
        let (mut matrix, mut header) = handle.get_alignment(&test_file);
        handle.remove_gappy_sequences(&mut matrix, &mut header, &threshold);
        assert_eq!(matrix.len(), 1);
    }

    #[test]
    fn test_threshold_calculation() {
        let input_dir = Path::new("tests/files/gappy");
        let params = SequenceFilteringOptions::GapThreshold(0.5);
        setup!(input_dir, handle, params, output);
        let threshold = 0.5;
        let alignment_len = 11;
        let max_gaps = handle.count_max_gaps(alignment_len, threshold);
        assert_eq!(max_gaps, 5);
        let threshold = 0.5;
        let alignment_len = 10;
        let max_gaps = handle.count_max_gaps(alignment_len, threshold);
        assert_eq!(max_gaps, 5);
    }

    #[test]
    fn test_floor_calculation() {
        let float_1: f64 = 5.1;
        let float_2: f64 = 5.5;
        let float_3: f64 = 5.9;
        let res_1 = float_1.floor() as usize;
        let res_2 = float_2.floor() as usize;
        let res_3 = float_3.floor() as usize;
        assert_eq!(res_1, 5);
        assert_eq!(res_2, 5);
        assert_eq!(res_3, 5);
    }
}
