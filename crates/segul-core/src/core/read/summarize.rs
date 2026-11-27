//! Generate summary statistics for raw-read sequences
//!
//! Support FASTQ and compressed FASTQ in gunzip format.

use std::path::{Path, PathBuf};

use anyhow::Result;
use colored::Colorize;
use indicatif::ProgressBar;
use rayon::prelude::*;

use crate::report::spec::{ReadSummaryParams, TaskKind, canonical};
use crate::report::{ReportContext, RunReportBuilder, report_path};
use crate::writer::read::{ReadDistSummaryWriter, ReadPosSummaryWriter, ReadSummaryWriter};
use crate::{
    helper::{
        progress::{self, ProgressReporter},
        types::{PhredOffset, SeqReadFmt, SummaryMode},
        utils::set_spinner,
    },
    stats::fastq::{
        DEFAULT_DUP_SKETCH_SIZE, FastqAnalyzer, FastqSummary, FastqSummaryMin,
        count_reads_with_progress,
        qc::QcProfile,
        report::{ReadMetrics, ReadReport, SampleReadReport, sample_name},
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
/// It accepts three modes to generate summary statistics:
/// 1. Minimal: read count only
/// 2. Basic: per-file and per-read summary statistics in a single CSV
/// 3. Complete (default): all Basic statistics plus per-position and
///    relative-position profiles, length x quality, distributions,
///    duplication, overrepresented sequences, adapter content, and JSON report.
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
    /// Quality-check profile stored in `report.json`.
    pub qc_profile: QcProfile,
    /// Receives per-file progress and can cancel the run.
    progress: Option<&'a dyn ProgressReporter>,
    report_context: ReportContext,
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
            qc_profile: QcProfile::Auto,
            progress: None,
            report_context: ReportContext::default(),
        }
    }

    /// Attach frontend provenance and SHA-256 checksum settings for `report.json`.
    pub fn with_report_context(mut self, context: ReportContext) -> Self {
        self.report_context = context;
        self
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

    /// Set the QC profile recorded in `read-report.json`.
    pub fn with_qc_profile(mut self, profile: QcProfile) -> Self {
        self.qc_profile = profile;
        self
    }

    /// Report per-file progress to `reporter`, which can also cancel the run.
    pub fn with_progress(mut self, reporter: &'a dyn ProgressReporter) -> Self {
        self.progress = Some(reporter);
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
    ///     &SummaryMode::Basic,
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
        let params = ReadSummaryParams {
            input_format: canonical::read_fmt(self.input_fmt).to_string(),
            mode: canonical::summary_mode(self.mode).to_string(),
            phred: canonical::phred(self.phred).to_string(),
            dup_sketch_size: self.dup_sketch_size,
            qc_profile: self.qc_profile.to_string(),
            prefix: self.prefix.map(str::to_string),
        };
        let mut sorted_inputs = self.inputs.to_vec();
        alphanumeric_sort::sort_path_slice(&mut sorted_inputs);
        let builder = RunReportBuilder::new(TaskKind::ReadSummary, &self.report_context)
            .primary_inputs(&sorted_inputs)
            .parameters(&params)?;
        let spin = set_spinner();
        spin.set_message("Calculating summary of fastq files");
        progress::report_start(self.progress, sorted_inputs.len());
        let files = match self.mode {
            SummaryMode::Minimal => {
                let records = self.par_summarize_minimal(&sorted_inputs)?;
                let summary = self.write_record_min(&spin, &records);
                let report = self.build_minimal_report(&records);
                let report_out = report_path(self.output, self.prefix);
                builder
                    .outputs(self.output, [("summary", summary.clone())])
                    .statistics(&report)?
                    .write_to(&report_out)?;
                ReadSummaryFiles {
                    summary,
                    mode: *self.mode,
                    positions: None,
                    distributions: None,
                    report: Some(report_out),
                }
            }
            SummaryMode::Basic => {
                let records = self.par_summarize(&sorted_inputs, false)?;
                let summary = self.write_record_basic(&spin, &records, false);
                let report = ReadReport::from_samples(
                    records.into_iter().map(|r| r.report).collect(),
                    self.qc_profile,
                );
                let report_out = report_path(self.output, self.prefix);
                builder
                    .outputs(self.output, [("summary", summary.clone())])
                    .statistics(&report)?
                    .write_to(&report_out)?;
                ReadSummaryFiles {
                    summary,
                    mode: *self.mode,
                    positions: None,
                    distributions: None,
                    report: Some(report_out),
                }
            }
            SummaryMode::Complete => {
                let records = self.par_summarize(&sorted_inputs, true)?;
                self.write_record_complete(&spin, records, builder)?
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

    fn par_summarize(&self, inputs: &[PathBuf], complete: bool) -> Result<Vec<FastqSummary>> {
        let threads = self.file_threads();
        let mut records = inputs
            .par_iter()
            .map(|path| {
                progress::track(self.progress, path, || {
                    FastqAnalyzer::new(path, self.phred, complete, threads)
                        .with_dup_sketch_size(self.dup_sketch_size)
                        .with_progress(self.progress)
                        .analyze()
                        .map(|report| FastqSummary::from_report(path, report))
                })
            })
            .collect::<Result<Vec<FastqSummary>>>()?;
        records.sort_by(|a, b| alphanumeric_sort::compare_path(&a.path, &b.path));
        self.print_totals(
            records
                .iter()
                .map(|r| (r.metrics().read_count, r.metrics().base_count)),
        );
        Ok(records)
    }

    fn par_summarize_minimal(&self, inputs: &[PathBuf]) -> Result<Vec<FastqSummaryMin>> {
        let threads = self.file_threads();
        let mut records = inputs
            .par_iter()
            .map(|path| {
                progress::track(self.progress, path, || {
                    let mut summary = FastqSummaryMin::new(path);
                    summary.read_count = count_reads_with_progress(path, threads, self.progress)?;
                    Ok(summary)
                })
            })
            .collect::<Result<Vec<FastqSummaryMin>>>()?;
        records.sort_by(|a, b| alphanumeric_sort::compare_path(&a.path, &b.path));
        Ok(records)
    }

    fn build_minimal_report(&self, records: &[FastqSummaryMin]) -> ReadReport {
        let samples = records
            .iter()
            .map(|rec| SampleReadReport {
                file_path: rec.path.display().to_string(),
                file_name: rec.file_name.clone(),
                sample: sample_name(&rec.path),
                metrics: ReadMetrics {
                    read_count: rec.read_count as u64,
                    ..Default::default()
                },
                ..Default::default()
            })
            .collect();
        ReadReport::from_samples(samples, self.qc_profile)
    }

    fn write_record_min(&self, spin: &ProgressBar, records: &[FastqSummaryMin]) -> PathBuf {
        let writer = ReadSummaryWriter::new(self.output, self.prefix);
        spin.set_message("Writing records\n");
        progress::report_stage(self.progress, "Writing records");
        writer
            .write_read_count_only(records)
            .expect("Failed writing to file")
    }

    fn write_record_basic(
        &self,
        spin: &ProgressBar,
        records: &[FastqSummary],
        complete: bool,
    ) -> PathBuf {
        spin.set_message("Writing records\n");
        progress::report_stage(self.progress, "Writing records");
        let writer = ReadSummaryWriter::new(self.output, self.prefix);
        writer
            .write_summary(records, complete)
            .expect("Failed writing to file")
    }

    fn write_record_complete(
        &self,
        spin: &ProgressBar,
        records: Vec<FastqSummary>,
        builder: RunReportBuilder,
    ) -> Result<ReadSummaryFiles> {
        let summary = self.write_record_basic(spin, &records, true);
        let report = ReadReport::from_samples(
            records.into_iter().map(|r| r.report).collect(),
            self.qc_profile,
        );
        self.print_cohort_info(&report);
        let samples: &[SampleReadReport] = &report.samples;
        let positions = ReadPosSummaryWriter::new(self.output, self.prefix).write(samples)?;
        let distributions =
            ReadDistSummaryWriter::new(self.output, self.prefix).write_report(&report)?;
        let report_out = report_path(self.output, self.prefix);
        builder
            .outputs(
                self.output,
                [
                    ("summary", summary.clone()),
                    ("positions", positions.clone()),
                    ("distributions", distributions.clone()),
                ],
            )
            .statistics(&report)?
            .write_to(&report_out)?;
        Ok(ReadSummaryFiles {
            summary,
            mode: *self.mode,
            positions: Some(positions),
            distributions: Some(distributions),
            report: Some(report_out),
        })
    }

    fn print_cohort_info(&self, report: &ReadReport) {
        if !report.pairs.is_empty() {
            let out_of_sync = report.pairs.iter().filter(|p| !p.in_sync).count();
            log::info!("{:18}: {}", "Matched pairs", report.pairs.len());
            if out_of_sync > 0 {
                log::warn!(
                    "{:18}: {} pair(s) have mismatched R1/R2 read counts",
                    "Pair sync",
                    out_of_sync
                );
            }
        }
        if let Some(cohort) = &report.cohort
            && cohort.sample_count > 1
        {
            log::info!("{:18}: {:.2}", "Depth CV", cohort.depth_cv);
            log::info!("{:18}: {:.2}", "Depth Gini", cohort.depth_gini);
            if !cohort.low_yield_samples.is_empty() {
                log::warn!(
                    "{:18}: {}",
                    "Low-yield samples",
                    cohort.low_yield_samples.join(", ")
                );
            }
        }
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
    use std::path::{Path, PathBuf};
    use std::sync::atomic::{AtomicUsize, Ordering};

    use tempfile::TempDir;

    use crate::core::read::summarize::GenomicReadSummary;
    use crate::helper::progress::{CancelAfter, Cancelled, ProgressReporter};
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
            &SummaryMode::Basic,
            output.path(),
            None,
        );
        handler.summarize();
        assert!(output.path().exists());
    }

    #[derive(Default)]
    struct CountingReporter {
        total: AtomicUsize,
        finished: AtomicUsize,
        cancelled: bool,
    }

    impl ProgressReporter for CountingReporter {
        fn on_start(&self, total_files: usize) {
            self.total.store(total_files, Ordering::Relaxed);
        }

        fn on_file_finished(&self, _path: &Path) {
            self.finished.fetch_add(1, Ordering::Relaxed);
        }

        fn is_cancelled(&self) -> bool {
            self.cancelled
        }
    }

    #[test]
    fn test_progress() {
        let mut files = vec![
            PathBuf::from("tests/files/raw/read_1.fastq"),
            PathBuf::from("tests/files/raw/read_2.fastq"),
        ];
        let output = TempDir::with_prefix("tempt").unwrap();
        let reporter = CountingReporter::default();
        GenomicReadSummary::new(
            &mut files,
            &SeqReadFmt::Auto,
            &SummaryMode::Basic,
            output.path(),
            None,
        )
        .with_progress(&reporter)
        .try_summarize()
        .unwrap();
        assert_eq!(reporter.total.load(Ordering::Relaxed), 2);
        assert_eq!(reporter.finished.load(Ordering::Relaxed), 2);
    }

    #[test]
    fn test_progress_cancelled() {
        let mut files = vec![PathBuf::from("tests/files/raw/read_1.fastq")];
        let output = TempDir::with_prefix("tempt").unwrap();
        let reporter = CountingReporter {
            cancelled: true,
            ..Default::default()
        };
        let error = GenomicReadSummary::new(
            &mut files,
            &SeqReadFmt::Auto,
            &SummaryMode::Basic,
            output.path(),
            None,
        )
        .with_progress(&reporter)
        .try_summarize()
        .unwrap_err();
        assert!(error.downcast_ref::<Cancelled>().is_some());
        assert_eq!(reporter.finished.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn test_progress_cancelled_mid_file() {
        for mode in [SummaryMode::Minimal, SummaryMode::Complete] {
            let mut files = vec![PathBuf::from("tests/files/raw/read_1.fastq")];
            let output = TempDir::with_prefix("tempt").unwrap();
            // The first check lets the file start; the next one is inside the file.
            let reporter = CancelAfter::new(1);
            let error =
                GenomicReadSummary::new(&mut files, &SeqReadFmt::Auto, &mode, output.path(), None)
                    .with_progress(&reporter)
                    .try_summarize()
                    .unwrap_err();
            assert!(error.downcast_ref::<Cancelled>().is_some());
            assert_eq!(reporter.finished.load(Ordering::Relaxed), 0);
        }
    }

    #[test]
    fn test_read_count_only() {
        let files = vec![
            PathBuf::from("tests/files/raw/read_1.fastq"),
            PathBuf::from("tests/files/raw/read_2.fastq"),
        ];
        let output = TempDir::with_prefix("tempt").unwrap();
        let mut files_copy = files.clone();
        let handler = GenomicReadSummary::new(
            &mut files_copy,
            &SeqReadFmt::Auto,
            &SummaryMode::Minimal,
            output.path(),
            None,
        );
        let records = handler.par_summarize_minimal(&files).unwrap();
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
        assert!(files.summary.ends_with("run_read-summary.csv"));
        assert!(files.positions.unwrap().exists());
        assert!(files.distributions.unwrap().exists());
        let report = ReadReport::from_path(&files.report.unwrap()).unwrap();
        // Samples are sorted by path.
        assert_eq!(report.samples[0].sample, "read_1");
        assert_eq!(report.samples[0].positions.len(), 36);
    }
}
