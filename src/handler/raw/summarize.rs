//! A handler for summarizing raw sequence data.

use std::path::PathBuf;

// use noodles::fastq;

use crate::helper::types::{RawReadFmt, SummaryMode};

/// Include support for any compressed or uncompressed fastq files.
pub struct RawSummaryHandler<'a> {
    pub inputs: &'a [PathBuf],
    pub input_fmt: &'a RawReadFmt,
    pub mode: &'a SummaryMode,
}

impl<'a> RawSummaryHandler<'a> {
    pub fn new(inputs: &'a [PathBuf], input_fmt: &'a RawReadFmt, mode: &'a SummaryMode) -> Self {
        Self {
            inputs,
            input_fmt,
            mode,
        }
    }

    pub fn summarize(&self) {
        self.inputs.iter().for_each(|p| println!("{}", p.display()))
    }
}
