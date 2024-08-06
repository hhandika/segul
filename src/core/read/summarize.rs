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
    stats::fastq::{FastqMappedRead, FastqSummary, FastqSummaryMin},
    writer::read::{ReadPosSummaryWriter, ReadSummaryWriter},
};

/// Generate read summary statistics.
///
/// It accept input in FASTQ or compressed FASTQ (.gzip).
/// The resulting file is in Comma-separated (.csv) format.
/// It accepts three mode to generate summary statistics:
/// 1. Minimal: read count only
/// 2. Default: essential statistics, such as read counts, base counts, gc, at, and n content, and qscore statistics
/// 3. Complete: all the essential plus summary
/// statistics per position in read for each file.
pub struct GenomicReadSummary<'a> {
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
    /// Output file prefix
    pub prefix: Option<&'a str>,
}

impl<'a> GenomicReadSummary<'a> {
    /// Create a new ReadSummaryHandler instance.
    pub fn new(
        inputs: &'a mut [PathBuf],
        input_fmt: &'a SeqReadFmt,
        mode: &'a SummaryMode,
        output: &'a Path,
        prefix: Option<&'a str>,
    ) -> Self {
        Self {
            inputs,
            input_fmt,
            mode,
            output,
            prefix,
        }
    }

    /// Generate summary statistics for fastq files.
    /// # Arguments
    /// * `path` - A mutable slice of PathBuf that holds the fastq files.
    /// * `input_fmt` - The fastq input format.
    /// * `mode` - The summary mode.
    /// * `output` - The output path.
    /// # Example
    /// ```rust
    /// use std::path::{Path, PathBuf};
    /// use segul::core::read::summarize::GenomicReadSummary;
    /// use segul::helper::types::{SeqReadFmt, SummaryMode};
    /// use tempdir::TempDir;
    ///
    /// let mut files = vec![
    ///    PathBuf::from("tests/files/raw/read_1.fastq"),
    ///    PathBuf::from("tests/files/raw/read_2.fastq"),
    /// ];
    /// let output = TempDir::new("tempt").unwrap();
    /// let spinner = segul::helper::utils::set_spinner();
    /// let mut handle = GenomicReadSummary::new(
    ///     &mut files,
    ///     &SeqReadFmt::Auto,
    ///     &SummaryMode::Default,
    ///     Path::new(output.path()),
    ///     None,
    /// );
    /// handle.summarize();
    /// ```
    pub fn summarize(&self) {
        let spin = set_spinner();
        spin.set_message("Calculating summary of fastq files");
        match self.mode {
            SummaryMode::Minimal => {
                let mut records = self.par_summarize_minimal();
                self.write_record_min(&spin, &mut records);
            }
            SummaryMode::Default => {
                let mut records = self.par_summarize_default();
                self.write_record_default(&spin, &mut records);
            }
            SummaryMode::Complete => {
                let all_records = self.par_summarize_complete();
                let (mut records, read_records): (Vec<FastqSummary>, Vec<FastqMappedRead>) =
                    all_records.into_iter().unzip();
                self.write_record_complete(&spin, &mut records, &read_records);
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

    fn summarize_default(&self, path: &Path) -> FastqSummary {
        let mut summary = FastqSummary::new(path);
        summary.summarize(self.input_fmt);
        summary
    }

    fn par_summarize_complete(&self) -> Vec<(FastqSummary, FastqMappedRead)> {
        let (sender, receiver) = channel();

        self.inputs.par_iter().for_each_with(sender, |s, p| {
            let record = self.summarize_complete(p);
            s.send(record)
                .expect("Failed parallel processing fastq files");
        });

        receiver.iter().collect()
    }

    fn summarize_complete(&self, path: &Path) -> (FastqSummary, FastqMappedRead) {
        let mut summary = FastqSummary::new(path);
        let mapped_records = summary.summarize_map(self.input_fmt);
        (summary, mapped_records)
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

    fn summarize_minimal(&self, p: &Path, input_fmt: &SeqReadFmt) -> FastqSummaryMin {
        let mut summary = FastqSummaryMin::new(p);
        summary.summarize(input_fmt);
        summary
    }

    fn write_record_min(&self, spin: &ProgressBar, records: &mut [FastqSummaryMin]) {
        let writer = ReadSummaryWriter::new(self.output, self.prefix);
        spin.set_message("Writing records\n");
        writer
            .write_read_count_only(records)
            .expect("Failed writing to file");
    }

    fn write_record_default(&self, spin: &ProgressBar, records: &mut [FastqSummary]) {
        // Sort records by file name
        records.sort_by(|a, b| a.path.cmp(&b.path));
        spin.set_message("Writing records\n");
        let writer = ReadSummaryWriter::new(self.output, self.prefix);
        writer.write(records).expect("Failed writing to file");
    }

    fn write_record_complete(
        &self,
        spin: &ProgressBar,
        records: &mut [FastqSummary],
        read_records: &[FastqMappedRead],
    ) {
        // Sort records by file name
        records.sort_by(|a, b| a.path.cmp(&b.path));
        let writer = ReadSummaryWriter::new(self.output, self.prefix);
        writer.write(records).expect("Failed writing to file");

        spin.set_message("Writing records\n");
        let pos_writer = ReadPosSummaryWriter::new(self.output, self.prefix);
        pos_writer
            .write(read_records)
            .expect("Failed writing to file");
    }

    fn print_output_info(&self) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Dir", self.output.display());
    }
}

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use tempdir::TempDir;

    use crate::core::read::summarize::GenomicReadSummary;
    use crate::helper::types::{SeqReadFmt, SummaryMode};
    use crate::stats::fastq::{FastqMappedRead, FastqSummary};

    #[test]
    fn test_summarize() {
        let mut files = vec![
            PathBuf::from("tests/files/raw/read_1.fastq"),
            PathBuf::from("tests/files/raw/read_2.fastq"),
        ];
        let output = TempDir::new("tempt").unwrap();
        let handler = GenomicReadSummary::new(
            &mut files,
            &SeqReadFmt::Auto,
            &SummaryMode::Default,
            output.path(),
            None,
        );
        handler.summarize();
        assert!(output.path().exists());
    }

    #[test]
    fn test_read_count_only() {
        let mut files = vec![
            PathBuf::from("tests/files/raw/read_1.fastq"),
            PathBuf::from("tests/files/raw/read_2.fastq"),
        ];
        let output = TempDir::new("tempt").unwrap();
        let handler = GenomicReadSummary::new(
            &mut files,
            &SeqReadFmt::Auto,
            &SummaryMode::Minimal,
            output.path(),
            None,
        );
        let records = handler.par_summarize_complete();
        let (_, pos): (Vec<FastqSummary>, Vec<FastqMappedRead>) = records.into_iter().unzip();
        pos.iter().for_each(|p| {
            assert_eq!(p.reads.len(), 36);
            assert_eq!(p.qscores.len(), 36);
        });

        assert_eq!(pos.len(), 2);
    }
}
