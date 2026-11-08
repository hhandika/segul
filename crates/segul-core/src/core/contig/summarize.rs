//! Handler for summarizing contig files
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};
use colored::Colorize;
use rayon::prelude::*;

use crate::writer::contigs::ContigSummaryWriter;
use crate::{
    helper::{
        progress::{self, ProgressReporter},
        types::ContigFmt,
        utils::set_spinner,
    },
    stats::contigs::ContigSummary,
};

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
        }
    }

    /// Report per-file progress to `reporter`, which can also cancel the run.
    pub fn with_progress(mut self, reporter: &'a dyn ProgressReporter) -> Self {
        self.progress = Some(reporter);
        self
    }

    /// Summarize the contig files. Returns the path of the summary file.
    pub fn summarize(&self) -> PathBuf {
        self.try_summarize()
            .unwrap_or_else(|error| panic!("Failed summarizing contigs: {:#}", error))
    }

    /// Same as [`Self::summarize`], but returns errors, including
    /// [`progress::Cancelled`] when the progress reporter stops the run.
    pub fn try_summarize(&self) -> Result<PathBuf> {
        let spin = set_spinner();
        spin.set_message("Calculating summary of contig files");
        let contig_summary = self.summarize_contigs()?;
        let writer = ContigSummaryWriter::new(&contig_summary, self.output, self.prefix);
        spin.set_message("Writing records\n");
        progress::report_stage(self.progress, "Writing records");
        let output = writer.write()?;
        spin.finish_with_message("Finished processing contig files\n");
        self.print_input_info();
        Ok(output)
    }

    fn summarize_contigs(&self) -> Result<Vec<ContigSummary>> {
        if self.files.is_empty() {
            bail!("No contig files found")
        }
        progress::report_start(self.progress, self.files.len());
        self.files
            .par_iter()
            .map(|path| progress::track(self.progress, path, || Ok(self.process_contigs(path))))
            .collect()
    }

    fn process_contigs(&self, input: &Path) -> ContigSummary {
        let mut summary = ContigSummary::new();
        summary.summarize(input, self.input_fmt);
        summary
    }

    fn print_input_info(&self) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Dir", self.output.display());
    }
}
