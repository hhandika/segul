//! Filter sequence based on user-defined criteria.
//!
//! Unlike the `alignment filter` that works on the alignment level,
//! the `sequence filter` will filter sequences within the alignment.
//! This is useful when you want to remove sequences that are too short
//! or have too many gaps.

use std::{
    path::{Path, PathBuf},
    sync::atomic::{AtomicUsize, Ordering},
};

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
    /// Total percentage of gaps in a sequence.
    /// The sequence will be removed
    /// if the percentage of gaps is higher than the threshold.
    /// Calculated as the number of gaps divided by the sequence length.
    GapThreshold(f64),
    /// Minimum sequence length.
    /// Filter sequence that has an equal or more sequence length
    /// than the specified number.
    MinSequenceLength(usize),
}

pub struct SequenceFiltering<'a> {
    files: &'a [PathBuf],
    input_fmt: &'a InputFmt,
    datatype: &'a DataType,
    output: &'a Path,
    output_fmt: &'a OutputFmt,
    params: &'a SequenceFilteringOptions,
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
        assert!(filtered_aln > 0, "No alignments left after filtering!");
        spinner.finish_with_message("Finished filtering sequences!\n");
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

    // Remove sequences that have a percentage of gaps higher than the threshold.
    // The threshold is calculated as the number of gaps times the sequence length.
    // For example, if the threshold is 0.5, and the sequence length is 10,
    // the sequence will be removed if it has more than 5 gaps.
    // When the result is in decimal, it will be floor to the nearest integer.
    // For example, if the threshold is 0.5 and the sequence length is 11,
    // the threshold will be 5.5, and the sequence will be removed if it has more than 5 gaps.
    fn remove_gappy_sequences(&self, matrix: &mut SeqMatrix, header: &mut Header, threshold: &f64) {
        let alignment_len = header.nchar;
        let max_gaps = self.count_max_gaps(alignment_len, *threshold);
        let mut to_remove = Vec::new();
        // Find sequences with gaps higher than the threshold.
        matrix.iter().for_each(|(name, seq)| {
            let gaps = seq.bytes().filter(|&c| c == b'-' || c == b'?').count();
            if gaps > max_gaps {
                to_remove.push(String::from(name));
            }
        });

        to_remove.iter().for_each(|name| {
            matrix.remove(name);
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

    #[test]
    fn test_filter_sequences_by_length() {
        let dir = "tests/files/concat";
        let files = SeqFileFinder::new(Path::new(dir)).find(&InputFmt::Nexus);
        let input_fmt = InputFmt::Nexus;
        let datatype = DataType::Dna;
        let output = TempDir::new("output").unwrap();
        let params = SequenceFilteringOptions::MinSequenceLength(7);
        let handle = SequenceFiltering::new(
            &files,
            &input_fmt,
            &datatype,
            output.path(),
            &OutputFmt::Nexus,
            &params,
        );
        handle.filter();
        let files = SeqFileFinder::new(output.path()).find(&InputFmt::Nexus);
        assert_eq!(files.len(), 1);
    }
}
