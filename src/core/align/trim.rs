//! Remove sites given a threshold of missing data.

use std::path::{Path, PathBuf};

use crate::helper::types::{DataType, InputFmt, OutputFmt};

pub enum TrimmingParameters {
    /// Trim based on a threshold of missing data
    MissingData(f64),
    /// Trim based on a threshold of Parsimony Informative sites (PIS)
    /// PIS is the number of sites that have at least two different states
    /// and at least two of the states have a minimum frequency of 2
    ParsInf(usize),
}

/// Trim alignment data structure
pub struct AlignmentTrimming<'a> {
    /// List of input files
    pub input_files: &'a [PathBuf],
    /// Input format of alignment files
    pub input_fmt: &'a InputFmt,
    /// Data type of sequences
    pub datatype: &'a DataType,
    /// Output directory
    pub output_dir: &'a Path,
    /// Output format of sequence files
    pub output_fmt: &'a OutputFmt,
    /// Trimming parameters
    pub params: &'a TrimmingParameters,
}

impl Default for AlignmentTrimming<'_> {
    fn default() -> Self {
        Self {
            input_files: &[],
            input_fmt: &InputFmt::Fasta,
            datatype: &DataType::Dna,
            output_dir: Path::new(""),
            output_fmt: &OutputFmt::Fasta,
            params: &TrimmingParameters::MissingData(0.1),
        }
    }
}

impl<'a> AlignmentTrimming<'a> {
    /// Create a new AlignmentTrimming instance
    pub fn new(
        input_files: &'a [PathBuf],
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        output_dir: &'a Path,
        output_fmt: &'a OutputFmt,
        params: &'a TrimmingParameters,
    ) -> Self {
        Self {
            input_files,
            input_fmt,
            datatype,
            output_dir,
            output_fmt,
            params,
        }
    }

    /// Trim alignment files based on the given parameters
    pub fn trim(&self) {
        match self.params {
            TrimmingParameters::MissingData(threshold) => {
                println!(
                    "Trimming based on missing data with threshold: {}",
                    threshold
                );
            }
            TrimmingParameters::ParsInf(threshold) => {
                println!(
                    "Trimming based on Parsimony Informative Sites with threshold: {}",
                    threshold
                );
            }
        }
    }
}
