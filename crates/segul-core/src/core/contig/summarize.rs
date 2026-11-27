//! Handler for summarizing contig files
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};
use colored::Colorize;
use rayon::prelude::*;

use crate::report::spec::{ContigSummaryParams, TaskKind, canonical};
use crate::report::{ReportContext, RunReportBuilder, report_path};
use crate::writer::contigs::ContigSummaryWriter;
use crate::{
    helper::{
        progress::{self, ProgressReporter},
        types::ContigFmt,
        utils::set_spinner,
    },
    stats::{contig_report::ContigReport, contigs::ContigSummary},
};

/// Files written by [`ContigSummaryHandler::summarize`].
#[derive(Debug, Clone)]
pub struct ContigSummaryFiles {
    /// Per-file statistics (`contig-summary.csv`).
    pub summary: PathBuf,
    /// Unified provenance and statistics JSON (`report.json`).
    pub report: PathBuf,
}

pub struct ContigSummaryHandler<'a> {
    /// Contig files
    files: &'a [PathBuf],
    /// Contig sequence format
    input_fmt: &'a ContigFmt,
    /// Output path
    output: &'a Path,
    /// Output file prefix
    prefix: Option<&'a str>,
    /// Receives per-file progress and can cancel the run.
    progress: Option<&'a dyn ProgressReporter>,
    report_context: ReportContext,
}

impl<'a> ContigSummaryHandler<'a> {
    pub fn new(
        files: &'a [PathBuf],
        input_fmt: &'a ContigFmt,
        output: &'a Path,
        prefix: Option<&'a str>,
    ) -> Self {
        Self {
            files,
            input_fmt,
            output,
            prefix,
            progress: None,
            report_context: ReportContext::default(),
        }
    }

    /// Attach frontend provenance and SHA-256 checksum settings for `report.json`.
    pub fn with_report_context(mut self, context: ReportContext) -> Self {
        self.report_context = context;
        self
    }

    /// Report per-file progress to `reporter`, which can also cancel the run.
    pub fn with_progress(mut self, reporter: &'a dyn ProgressReporter) -> Self {
        self.progress = Some(reporter);
        self
    }

    /// Summarize the contig files. Returns the written files.
    pub fn summarize(&self) -> ContigSummaryFiles {
        self.try_summarize()
            .unwrap_or_else(|error| panic!("Failed summarizing contigs: {:#}", error))
    }

    /// Same as [`Self::summarize`], but returns errors, including
    /// [`progress::Cancelled`] when the progress reporter stops the run.
    pub fn try_summarize(&self) -> Result<ContigSummaryFiles> {
        let params = ContigSummaryParams {
            input_format: canonical::contig_fmt(self.input_fmt).to_string(),
            prefix: self.prefix.map(str::to_string),
        };
        let mut sorted_files = self.files.to_vec();
        alphanumeric_sort::sort_path_slice(&mut sorted_files);
        let builder = RunReportBuilder::new(TaskKind::ContigSummary, &self.report_context)
            .primary_inputs(&sorted_files)
            .parameters(&params)?;
        let spin = set_spinner();
        spin.set_message("Calculating summary of contig files");
        let contig_summary = self.summarize_contigs(&sorted_files)?;
        let writer = ContigSummaryWriter::new(&contig_summary, self.output, self.prefix);
        spin.set_message("Writing records\n");
        progress::report_stage(self.progress, "Writing records");
        let summary = writer.write()?;
        let report = report_path(self.output, self.prefix);
        let contig_report = ContigReport::new(&contig_summary);
        builder
            .outputs(self.output, [("summary", summary.clone())])
            .statistics(&contig_report)?
            .write_to(&report)?;
        spin.finish_with_message("Finished processing contig files\n");
        self.print_output_info(&report);
        Ok(ContigSummaryFiles { summary, report })
    }

    fn summarize_contigs(&self, files: &[PathBuf]) -> Result<Vec<ContigSummary>> {
        if files.is_empty() {
            bail!("No contig files found")
        }
        progress::report_start(self.progress, files.len());
        files
            .par_iter()
            .map(|path| progress::track(self.progress, path, || Ok(self.process_contigs(path))))
            .collect()
    }

    fn process_contigs(&self, input: &Path) -> ContigSummary {
        let mut summary = ContigSummary::new();
        summary.summarize(input, self.input_fmt);
        summary
    }

    fn print_output_info(&self, report: &Path) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Dir", self.output.display());
        log::info!("{:18}: {}", "Report", report.display());
    }
}

#[cfg(test)]
mod test {
    use tempfile::TempDir;

    use super::*;

    #[test]
    fn test_summary_files() {
        let output = TempDir::with_prefix("contig").unwrap();
        let files: Vec<PathBuf> = ["contigs1.fa", "contigs2.fa"]
            .iter()
            .map(|f| PathBuf::from("tests/files/contigs").join(f))
            .collect();
        let written =
            ContigSummaryHandler::new(&files, &ContigFmt::Fasta, output.path(), Some("run"))
                .summarize();
        assert!(written.summary.is_file());
        assert_eq!(written.report, output.path().join("run_report.json"));
        let report = ContigReport::from_path(&written.report).unwrap();
        assert_eq!(report.summary.samples, 2);
        assert_eq!(report.samples[1].name, "contigs2");
    }
}
