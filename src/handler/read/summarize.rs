//! A handler for summarizing raw sequence data.
//!
//! The handler accept input in FASTQ or compressed FASTQ (.gzip)
//! It provide three mode to generate summary statistics
//! 1. Minimal: Read count only
//! 2. Default: Essential statistics,
//! such as read counts, base counts, gc, at, and n content, and qscore statistics

use std::{
    path::{Path, PathBuf},
    sync::mpsc::channel,
};

use colored::Colorize;

use indicatif::ProgressBar;
use rayon::prelude::*;

use crate::{
    helper::{
        types::{SeqReadFmt, SummaryMode},
        utils::set_spinner,
    },
    stats::fastq::{FastqSummary, FastqSummaryMin},
    writer::read::ReadSummaryWriter,
};

/// Include support for any compressed or uncompressed fastq files.
pub struct ReadSummaryHandler<'a> {
    pub inputs: &'a mut [PathBuf],
    pub input_fmt: &'a SeqReadFmt,
    /// Summary mode
    /// * `Minimal` - Only write the number of reads in each file
    /// * `Default` - Write all the essential summary statistics
    /// * `Complete` - Write all the summary statistics
    pub mode: &'a SummaryMode,
    pub output: &'a Path,
}

impl<'a> ReadSummaryHandler<'a> {
    pub fn new(
        inputs: &'a mut [PathBuf],
        input_fmt: &'a SeqReadFmt,
        mode: &'a SummaryMode,
        output: &'a Path,
    ) -> Self {
        Self {
            inputs,
            input_fmt,
            mode,
            output,
        }
    }

    pub fn summarize(&mut self) {
        let spin = set_spinner();
        spin.set_message("Calculating summary of fastq files");
        match self.mode {
            SummaryMode::Minimal => {
                let mut records = self.par_summarize_minimal();
                self.write_record_min(&spin, &mut records);
            }
            SummaryMode::Default => {
                let mut records = self.par_summarize_default();
                self.write_records(&spin, &mut records);
            }
            SummaryMode::Complete => {
                let mut records = self.par_summarize_complete();
                self.write_records(&spin, &mut records);
            }
        }
        spin.finish_with_message("Finished processing fastq files\n");
        self.print_output_info();
    }

    fn par_summarize_default(&self) -> Vec<FastqSummary> {
        let (sender, receiver) = channel();

        self.inputs.par_iter().for_each_with(sender, |s, p| {
            let record = self.summarize_default(p);
            s.send(record)
                .expect("Failed parallel processing fastq files");
        });

        receiver.iter().collect()
    }

    fn par_summarize_complete(&self) -> Vec<FastqSummary> {
        let (sender, receiver) = channel();

        self.inputs.par_iter().for_each_with(sender, |s, p| {
            let record = self.summarize_complete(p);
            s.send(record)
                .expect("Failed parallel processing fastq files");
        });

        receiver.iter().collect()
    }

    fn par_summarize_minimal(&self) -> Vec<FastqSummaryMin> {
        let (sender, receiver) = channel();

        self.inputs.par_iter().for_each_with(sender, |s, p| {
            let summary = self.summarize_minimal(p, self.input_fmt);
            s.send(summary)
                .expect("Failed parallel processing fastq files");
        });

        receiver.iter().collect()
    }

    fn write_record_min(&mut self, spin: &ProgressBar, records: &mut [FastqSummaryMin]) {
        let writer = ReadSummaryWriter::new(self.output);
        spin.set_message("Writing records\n");
        writer
            .write_read_count_only(&records)
            .expect("Failed writing to file");
    }

    fn write_records(&mut self, spin: &ProgressBar, records: &mut [FastqSummary]) {
        // Sort records by file name
        records.sort_by(|a, b| a.path.cmp(&b.path));
        spin.set_message("Writing records\n");
        let writer = ReadSummaryWriter::new(self.output);
        writer.write(records).expect("Failed writing to file");
    }

    fn summarize_minimal(&self, p: &Path, input_fmt: &SeqReadFmt) -> FastqSummaryMin {
        let mut summary = FastqSummaryMin::new(p);
        summary.summarize(input_fmt);
        summary
    }

    fn summarize_default(&self, path: &Path) -> FastqSummary {
        let mut summary = FastqSummary::new(path);
        summary.summarize(self.input_fmt);
        summary
    }

    fn summarize_complete(&self, path: &Path) -> FastqSummary {
        let mut summary = FastqSummary::new(path);
        let mapped_records = summary.summarize_map(self.input_fmt);
        let writer = ReadSummaryWriter::new(self.output);
        writer.write_per_read_records(path, &mapped_records.reads, &mapped_records.qscores);
        summary
    }

    fn print_output_info(&self) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Dir", self.output.display());
    }
}
