//! Generate summary statistics for raw-read sequences
//!
//! Support FASTQ and compressed FASTQ in gunzip format.

use std::path::{Path, PathBuf};

use anyhow::Result;
use colored::Colorize;
use indicatif::ProgressBar;
use rayon::prelude::*;

use crate::{
    helper::{
        types::{PhredOffset, SeqReadFmt, SummaryMode},
        utils::set_spinner,
    },
    stats::fastq::{
        DEFAULT_DUP_SKETCH_SIZE, FastqAnalyzer, FastqSummary, FastqSummaryMin, count_reads,
        report::{ReadReport, SampleReadReport},
    },
    writer::read::{
        ReadDistSummaryWriter, ReadPosSummaryWriter, ReadReportWriter, ReadSummaryWriter,
    },
};

// More threads than this per file rarely speed up BGZF decoding,
// and each one adds blocks held in memory.
const MAX_FILE_THREADS: usize = 8;

/// Files written by [`GenomicReadSummary::summarize`].
#[derive(Debug, Clone)]
pub struct ReadSummaryFiles {
    /// Per-file summary. Its columns depend on `mode`.
    pub summary: PathBuf,
    pub mode: SummaryMode,
    /// Per-position archive (`read-pos-summary.zip`), written in Complete mode.
    pub positions: Option<PathBuf>,
    /// Distribution archive (`read-dist-summary.zip`), written in Complete mode.
    pub distributions: Option<PathBuf>,
    /// All statistics as JSON (`read-report.json`), written in Complete mode.
    pub report: Option<PathBuf>,
}

/// Generate read summary statistics.
///
/// It accept input in FASTQ or compressed FASTQ (.gzip).
/// The resulting file is in Comma-separated (.csv) format.
/// It accepts three mode to generate summary statistics:
/// 1. Minimal: read count only
/// 2. Default: essential statistics, such as read counts, base counts, gc, at, and n content, and qscore statistics
/// 3. Complete: all the essential statistics plus statistics per read
///    position, read length, quality, and GC distributions, duplication,
///    overrepresented sequences, and adapter content for each file.
///
/// Each file is read once. Large files are split into blocks that are
/// summarized in parallel.
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
    /// Quality score encoding. Detected from the reads by default.
    pub phred: PhredOffset,
    /// Distinct sequences kept to estimate duplication in Complete mode.
    /// Zero counts every distinct sequence.
    pub dup_sketch_size: usize,
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
            phred: PhredOffset::Auto,
            dup_sketch_size: DEFAULT_DUP_SKETCH_SIZE,
        }
    }

    /// Keep this many distinct sequences to estimate duplication, or zero
    /// to count every one. See [`FastqAnalyzer::with_dup_sketch_size`].
    pub fn with_dup_sketch_size(mut self, size: usize) -> Self {
        self.dup_sketch_size = size;
        self
    }

    /// Use a fixed quality score offset instead of detecting it.
    pub fn with_phred(mut self, phred: PhredOffset) -> Self {
        self.phred = phred;
        self
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
    /// use segul_core::core::read::summarize::GenomicReadSummary;
    /// use segul_core::helper::types::{SeqReadFmt, SummaryMode};
    /// use tempfile::TempDir;
    ///
    /// let mut files = vec![
    ///    PathBuf::from("tests/files/raw/read_1.fastq"),
    ///    PathBuf::from("tests/files/raw/read_2.fastq"),
    /// ];
    /// let output = TempDir::with_prefix("tempt").unwrap();
    /// let spinner = segul_core::helper::utils::set_spinner();
    /// let mut handle = GenomicReadSummary::new(
    ///     &mut files,
    ///     &SeqReadFmt::Auto,
    ///     &SummaryMode::Default,
    ///     Path::new(output.path()),
    ///     None,
    /// );
    /// handle.summarize();
    /// ```
    pub fn summarize(&self) -> ReadSummaryFiles {
        self.try_summarize()
            .unwrap_or_else(|error| panic!("Failed summarizing reads: {:#}", error))
    }

    /// Same as [`Self::summarize`], but returns errors from invalid input files.
    pub fn try_summarize(&self) -> Result<ReadSummaryFiles> {
        let spin = set_spinner();
        spin.set_message("Calculating summary of fastq files");
        let files = match self.mode {
            SummaryMode::Minimal => {
                let records = self.par_summarize_minimal()?;
                ReadSummaryFiles {
                    summary: self.write_record_min(&spin, &records),
                    mode: *self.mode,
                    positions: None,
                    distributions: None,
                    report: None,
                }
            }
            SummaryMode::Default => {
                let records = self.par_summarize(false)?;
                ReadSummaryFiles {
                    summary: self.write_record_default(&spin, &records, false),
                    mode: *self.mode,
                    positions: None,
                    distributions: None,
                    report: None,
                }
            }
            SummaryMode::Complete => {
                let records = self.par_summarize(true)?;
                self.write_record_complete(&spin, records)?
            }
        };
        spin.finish_with_message("Finished processing fastq files\n");
        self.print_output_info();
        Ok(files)
    }

    /// Threads for each file when the pool is shared between files.
    fn file_threads(&self) -> usize {
        (rayon::current_num_threads() / self.inputs.len().max(1)).clamp(1, MAX_FILE_THREADS)
    }

    fn par_summarize(&self, complete: bool) -> Result<Vec<FastqSummary>> {
        let threads = self.file_threads();
        let mut records = self
            .inputs
            .par_iter()
            .map(|path| {
                FastqAnalyzer::new(path, self.phred, complete, threads)
                    .with_dup_sketch_size(self.dup_sketch_size)
                    .analyze()
                    .map(|report| FastqSummary::from_report(path, report))
            })
            .collect::<Result<Vec<FastqSummary>>>()?;
        records.sort_by(|a, b| a.path.cmp(&b.path));
        self.print_totals(
            records
                .iter()
                .map(|r| (r.metrics().read_count, r.metrics().base_count)),
        );
        Ok(records)
    }

    fn par_summarize_minimal(&self) -> Result<Vec<FastqSummaryMin>> {
        let threads = self.file_threads();
        let mut records = self
            .inputs
            .par_iter()
            .map(|path| {
                let mut summary = FastqSummaryMin::new(path);
                summary.read_count = count_reads(path, threads)?;
                Ok(summary)
            })
            .collect::<Result<Vec<FastqSummaryMin>>>()?;
        records.sort_by(|a, b| a.path.cmp(&b.path));
        Ok(records)
    }

    fn write_record_min(&self, spin: &ProgressBar, records: &[FastqSummaryMin]) -> PathBuf {
        let writer = ReadSummaryWriter::new(self.output, self.prefix);
        spin.set_message("Writing records\n");
        writer
            .write_read_count_only(records)
            .expect("Failed writing to file")
    }

    fn write_record_default(
        &self,
        spin: &ProgressBar,
        records: &[FastqSummary],
        complete: bool,
    ) -> PathBuf {
        spin.set_message("Writing records\n");
        let writer = ReadSummaryWriter::new(self.output, self.prefix);
        writer
            .write_summary(records, complete)
            .expect("Failed writing to file")
    }

    fn write_record_complete(
        &self,
        spin: &ProgressBar,
        records: Vec<FastqSummary>,
    ) -> Result<ReadSummaryFiles> {
        let summary = self.write_record_default(spin, &records, true);
        let report = ReadReport {
            samples: records.into_iter().map(|r| r.report).collect(),
        };
        let samples: &[SampleReadReport] = &report.samples;
        let positions = ReadPosSummaryWriter::new(self.output, self.prefix).write(samples)?;
        let distributions = ReadDistSummaryWriter::new(self.output, self.prefix).write(samples)?;
        let report = ReadReportWriter::new(self.output, self.prefix).write(&report)?;
        Ok(ReadSummaryFiles {
            summary,
            mode: *self.mode,
            positions: Some(positions),
            distributions: Some(distributions),
            report: Some(report),
        })
    }

    fn print_totals(&self, counts: impl Iterator<Item = (u64, u64)>) {
        let (reads, bases) = counts.fold((0, 0), |(r, b), (reads, bases)| (r + reads, b + bases));
        log::info!("{:18}: {}", "Total reads", reads);
        log::info!("{:18}: {}", "Total bases", bases);
    }

    fn print_output_info(&self) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Dir", self.output.display());
    }
}

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use tempfile::TempDir;

    use crate::core::read::summarize::GenomicReadSummary;
    use crate::helper::types::{SeqReadFmt, SummaryMode};
    use crate::stats::fastq::report::ReadReport;

    #[test]
    fn test_summarize() {
        let mut files = vec![
            PathBuf::from("tests/files/raw/read_1.fastq"),
            PathBuf::from("tests/files/raw/read_2.fastq"),
        ];
        let output = TempDir::with_prefix("tempt").unwrap();
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
        let output = TempDir::with_prefix("tempt").unwrap();
        let handler = GenomicReadSummary::new(
            &mut files,
            &SeqReadFmt::Auto,
            &SummaryMode::Minimal,
            output.path(),
            None,
        );
        let records = handler.par_summarize_minimal().unwrap();
        assert_eq!(records.len(), 2);
        assert!(records.iter().all(|r| r.read_count == 2));
    }

    #[test]
    fn test_complete_outputs() {
        let mut files = vec![
            PathBuf::from("tests/files/raw/read_2.fastq"),
            PathBuf::from("tests/files/raw/read_1.fastq"),
        ];
        let output = TempDir::with_prefix("tempt").unwrap();
        let handler = GenomicReadSummary::new(
            &mut files,
            &SeqReadFmt::Auto,
            &SummaryMode::Complete,
            output.path(),
            Some("run"),
        );
        let files = handler.try_summarize().unwrap();
        assert!(files.summary.ends_with("run_default-read-summary.csv"));
        assert!(files.positions.unwrap().exists());
        assert!(files.distributions.unwrap().exists());
        let report = ReadReport::from_path(&files.report.unwrap()).unwrap();
        // Samples are sorted by path.
        assert_eq!(report.samples[0].sample, "read_1");
        assert_eq!(report.samples[0].positions.len(), 36);
    }
}
