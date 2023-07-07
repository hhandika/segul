//! Generate summary statistics for raw-read sequences
//!
//! Support FASTQ and compressed FASTQ in gunzip format.

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

/// Generate read summary statistics.
///
/// It accept input in FASTQ or compressed FASTQ (.gzip).
/// The resulting file is in Comma-separated (.csv) format.
/// It accepts three mode to generate summary statistics:
/// 1. Minimal: read count only
/// 2. Default: essential statistics,
/// such as read counts, base counts, gc, at, and n content, and qscore statistics
/// 3. Complete: all the essential plus summary
/// statistics per position in read for each file.
pub struct ReadSummaryHandler<'a> {
    /// Input path.
    pub inputs: &'a mut [PathBuf],
    /// Read sequence format.
    /// In Auto, it will try to infer the format
    /// based on the file extension.
    /// For example, sequence_1.fastq.gz
    /// will result in SeqReadFmt::Gzip
    pub input_fmt: &'a SeqReadFmt,
    /// Summary statistic mode
    pub mode: &'a SummaryMode,
    /// Output path. No extension required
    pub output: &'a Path,
}

impl<'a> ReadSummaryHandler<'a> {
    /// Create a new ReadSummaryHandler instance.
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

    /// Generate summary statistics for fastq files.
    /// # Arguments
    /// * `path` - A mutable slice of PathBuf that holds the fastq files.
    /// * `input_fmt` - The fastq input format.
    /// * `mode` - The summary mode.
    /// * `output` - The output path.
    /// # Example
    /// ```
    /// use std::path::{Path, PathBuf};
    /// use segul::handler::read::summarize::ReadSummaryHandler;
    /// use segul::helper::types::{SeqReadFmt, SummaryMode};
    /// use tempdir::TempDir;
    ///
    /// let mut files = vec![
    ///    PathBuf::from("tests/files/raw/read_1.fastq"),
    ///    PathBuf::from("tests/files/raw/read_2.fastq"),
    /// ];
    /// let output = TempDir::new("tempt").unwrap();
    /// let spinner = segul::helper::utils::set_spinner();
    /// let mut handler = ReadSummaryHandler::new(
    ///     &mut files,
    ///     &SeqReadFmt::Auto,
    ///     &SummaryMode::Default,
    ///     Path::new(output.path()),
    /// );
    /// handler.summarize();
    /// ```
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
