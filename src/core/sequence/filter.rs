//! Filter sequence based on selected criteria.
//!
//! Unlike the `alignment filter` that works on the alignment level,
//! the `sequence filter` will filter sequences within the alignment.
//! The filtering criteria are:
//! - Total percentage of gaps in a sequence. The sequence will be removed
//!     if the percentage of gaps is higher than the threshold.
//! - Minimum sequence length. Filter sequence that has an equal or more sequence length
//!     than the specified number.

use std::{
    path::{Path, PathBuf},
    sync::atomic::{AtomicUsize, Ordering},
};

use colored::Colorize;
use rayon::prelude::*;

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

macro_rules! filter_by_length {
    ($self: ident, $length: ident, $filter: ident) => {{
        let counter = AtomicUsize::new(0);
        $self.files.par_iter().for_each(|file| {
            let (mut matrix, mut header) = $self.get_sequence_matrix(file);
            matrix.retain(|_, seq| $self.$filter(seq, $length));
            header.ntax = matrix.len();
            if header.ntax > 0 {
                $self.write_sequence(file, &matrix, &header);
                counter.fetch_add(1, Ordering::Relaxed);
            }
        });

        counter.into_inner()
    }};
}

/// Available sequence filtering options.
pub enum SeqFilteringParameters {
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
    PercentMaxGap(f64),
    /// Filter sequences based on the minimum sequence length without gaps.
    /// Remove sequences that have a sequence length less than the specified number.
    MinSequenceLength(usize),
    MaxSeequenceLength(usize),
    None,
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
    params: &'a SeqFilteringParameters,
}

impl OutputPrint for SequenceFiltering<'_> {}

impl<'a> SequenceFiltering<'a> {
    pub fn new(
        files: &'a [PathBuf],
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        output: &'a Path,
        output_fmt: &'a OutputFmt,
        params: &'a SeqFilteringParameters,
    ) -> Self {
        Self {
            files,
            input_fmt,
            datatype,
            output,
            output_fmt,
            params,
        }
    }

    /// Filter sequences based on the selected criteria.
    /// The resulting files are alignments with filtered sequences.
    /// Example:
    /// ```
    /// use std::path::{Path, PathBuf};
    /// use tempdir::TempDir;
    /// use segul::helper::types::{DataType, InputFmt, OutputFmt, PartitionFmt};
    /// use segul::core::sequence::filter::{SequenceFiltering, SeqFilteringParameters};
    /// use segul::helper::finder::SeqFileFinder;
    ///
    ///
    /// let input_fmt = InputFmt::Nexus;
    /// let datatype = DataType::Dna;
    /// let input_dir = Path::new("tests/files/alignments");
    /// let files = SeqFileFinder::new(Path::new(input_dir)).find(&input_fmt);
    /// // Replace the temp directory with your own directory.
    /// let output = TempDir::new("tempt").unwrap();
    /// let output_fmt = OutputFmt::Nexus;
    /// let params = SeqFilteringParameters::MinSequenceLength(7);
    /// let handle = SequenceFiltering::new(&files, &input_fmt, &datatype, &output.path(), &output_fmt,  &params);
    /// handle.filter();
    pub fn filter(&self) {
        let spinner = utils::set_spinner();
        spinner.set_message("Filtering sequences...");
        let filtered_aln = match self.params {
            SeqFilteringParameters::PercentMaxGap(threshold) => {
                self.filter_gappy_sequences(threshold)
            }
            SeqFilteringParameters::MinSequenceLength(min_length) => {
                self.filter_sequences_by_min_length(min_length)
            }
            SeqFilteringParameters::MaxSeequenceLength(max_length) => {
                self.filter_sequences_by_max_length(max_length)
            }
            SeqFilteringParameters::None => {
                log::warn!("No filtering parameters were provided!");
                0
            }
        };
        spinner.finish_with_message("Finished filtering sequences!\n");
        if filtered_aln == 0 {
            log::warn!("No matching sequences were found!");
        } else {
            self.print_output_info(filtered_aln);
        }
    }

    fn filter_gappy_sequences(&self, threshold: &f64) -> usize {
        let counter = AtomicUsize::new(0);
        self.files.par_iter().for_each(|file| {
            let (mut matrix, mut header) = self.get_sequence_matrix(file);
            self.remove_gappy_sequences(&mut matrix, &mut header, threshold);
            if header.ntax > 0 {
                self.write_sequence(file, &matrix, &header);
                counter.fetch_add(1, Ordering::Relaxed);
            }
        });

        counter.into_inner()
    }

    fn filter_sequences_by_min_length(&self, length: &usize) -> usize {
        filter_by_length!(self, length, is_longer_sequence)
    }

    fn filter_sequences_by_max_length(&self, length: &usize) -> usize {
        filter_by_length!(self, length, is_shorter_sequence)
    }

    fn is_longer_sequence(&self, sequence: &str, length: &usize) -> bool {
        self.count_non_gaps(sequence) >= *length
    }

    fn is_shorter_sequence(&self, sequence: &str, length: &usize) -> bool {
        self.count_non_gaps(sequence) <= *length
    }

    fn write_sequence(&self, file: &Path, matrix: &SeqMatrix, header: &Header) {
        let output_path = files::create_output_fname(self.output, file, self.output_fmt);
        let mut seq_writer = SeqWriter::new(&output_path, matrix, header);
        seq_writer
            .write_sequence(self.output_fmt)
            .expect("Failed writing the output file");
    }

    fn print_output_info(&self, counter: usize) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Directory", self.output.display());
        self.print_output_fmt(self.output_fmt);
        log::info!("{:18}: {}", "Total files", counter);
    }

    fn get_sequence_matrix(&self, file: &Path) -> (SeqMatrix, Header) {
        let sequence = SeqParser::new(file, self.datatype);
        sequence.parse(self.input_fmt)
    }

    fn remove_gappy_sequences(&self, matrix: &mut SeqMatrix, header: &mut Header, threshold: &f64) {
        let alignment_len = header.nchar;
        let max_gaps = self.count_max_gaps(alignment_len, *threshold);
        matrix.retain(|_, seq| self.count_gaps(seq) <= max_gaps);

        header.ntax = matrix.len();
    }

    fn count_non_gaps(&self, seq: &str) -> usize {
        seq.bytes().filter(|&c| c != b'-' && c != b'?').count()
    }

    fn count_gaps(&self, seq: &str) -> usize {
        seq.bytes().filter(|&c| c == b'-' || c == b'?').count()
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
        let dir = Path::new("tests/files/alignments");
        let params = SeqFilteringParameters::MinSequenceLength(7);
        setup!(dir, handle, params, output);
        handle.filter();
        let output_files = SeqFileFinder::new(output.path()).find(&InputFmt::Nexus);
        assert_eq!(output_files.len(), 1);
    }

    #[test]
    fn test_filter_sequences_by_gaps() {
        let dir = Path::new("tests/files/gappy");
        let params = SeqFilteringParameters::PercentMaxGap(0.5);
        setup!(dir, handle, params, output);
        handle.filter();
        let output_files = SeqFileFinder::new(output.path()).find(&InputFmt::Nexus);
        assert_eq!(output_files.len(), 4);
    }

    #[test]
    fn test_gap_calculation() {
        let input_dir = Path::new("tests/files/gappy");
        let params = SeqFilteringParameters::MinSequenceLength(7);
        setup!(input_dir, handle, params, output);
        let threshold = 0.5;
        let test_file = input_dir.join("gene_1.nex");
        let (mut matrix, mut header) = handle.get_sequence_matrix(&test_file);
        handle.remove_gappy_sequences(&mut matrix, &mut header, &threshold);
        assert_eq!(matrix.len(), 1);
    }

    #[test]
    fn test_threshold_calculation() {
        let input_dir = Path::new("tests/files/gappy");
        let params = SeqFilteringParameters::PercentMaxGap(0.5);
        setup!(input_dir, handle, params, output);
        let threshold = 0.5;
        let alignment_len = 11;
        let max_gaps = handle.count_max_gaps(alignment_len, threshold);
        assert_eq!(max_gaps, 5);
        let threshold = 0.5;
        let alignment_len = 10;
        let max_gaps = handle.count_max_gaps(alignment_len, threshold);
        assert_eq!(max_gaps, 5);
        let seq = "ataggata--??nn";
        let gaps = handle.count_gaps(seq);
        assert_eq!(gaps, 4);
        let non_gaps = handle.count_non_gaps(seq);
        assert_eq!(non_gaps, 10);
        let seq_2 = "ataggata--";
        let non_gap_2 = handle.count_non_gaps(seq_2);
        assert_eq!(non_gap_2, 8);
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
