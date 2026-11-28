//! Convert phylogenetic tree formats between Newick, Extended Newick, and NEXUS.

use std::path::{Path, PathBuf};

use colored::Colorize;
use rayon::prelude::*;

use crate::helper::progress::{self, ProgressReporter};
use crate::helper::utils;
use crate::parser::tree::{ParseOptions, TreeFmt, TreeOutputFmt, TreeReader};
use crate::report::spec::{TaskKind, TreeConvertParams, canonical};
use crate::report::stats::TaskCounts;
use crate::report::{DirSnapshot, ReportContext, RunReport, RunReportBuilder, report_path};
use crate::writer::tree::{TreeWriteOptions, TreeWriter};

/// Convert phylogenetic tree files between supported tree formats.
pub struct TreeConversion<'a> {
    input_fmt: &'a TreeFmt,
    output_fmt: &'a TreeOutputFmt,
    parse_opts: ParseOptions,
    write_opts: TreeWriteOptions,
    report_context: ReportContext,
    progress: Option<&'a dyn ProgressReporter>,
}

impl<'a> TreeConversion<'a> {
    /// Create a new conversion task from `input_fmt` to `output_fmt`.
    pub fn new(input_fmt: &'a TreeFmt, output_fmt: &'a TreeOutputFmt) -> Self {
        Self {
            input_fmt,
            output_fmt,
            parse_opts: ParseOptions::default(),
            write_opts: TreeWriteOptions::default(),
            report_context: ReportContext::default(),
            progress: None,
        }
    }

    /// Customize parser options (e.g., skipping annotations for maximum throughput).
    pub fn with_parse_options(mut self, opts: ParseOptions) -> Self {
        self.parse_opts = opts;
        self
    }

    /// Customize writer options (e.g., `use_translate`, `write_annotations`, `gzip`).
    pub fn with_write_options(mut self, opts: TreeWriteOptions) -> Self {
        self.write_opts = opts;
        self
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

    /// Convert `files` into `output` directory without emitting `report.json`.
    pub fn convert(&self, files: &[PathBuf], output: &Path) -> anyhow::Result<usize> {
        progress::report_start(self.progress, files.len());
        let mut total_trees = 0usize;
        let counts: Vec<usize> = files
            .par_iter()
            .map(|file| self.convert_single(file, output))
            .collect::<anyhow::Result<Vec<_>>>()?;
        for n in counts {
            total_trees += n;
        }
        Ok(total_trees)
    }

    /// Convert `files` into `output` directory and write a unified `report.json`.
    pub fn convert_with_report(
        &self,
        files: &[PathBuf],
        output: &Path,
    ) -> anyhow::Result<RunReport> {
        let mut sorted = files.to_vec();
        alphanumeric_sort::sort_path_slice(&mut sorted);
        let params = TreeConvertParams {
            input_format: canonical::tree_fmt(self.input_fmt).to_string(),
            output_format: canonical::tree_output_fmt(self.output_fmt).to_string(),
            translate: self.write_opts.use_translate,
            strip_annotations: !self.write_opts.write_annotations,
            gzip: self.write_opts.gzip,
        };
        let snapshot = DirSnapshot::capture(output);
        let builder = RunReportBuilder::new(TaskKind::TreeConvert, &self.report_context)
            .primary_inputs(&sorted)
            .parameters(&params)?;

        progress::report_start(self.progress, sorted.len());
        let spin = utils::set_spinner();
        spin.set_message("Converting phylogenetic trees...");
        let tree_counts: Vec<usize> = sorted
            .par_iter()
            .map(|file| self.convert_single(file, output))
            .collect::<anyhow::Result<Vec<_>>>()?;
        let total_trees: usize = tree_counts.into_iter().sum();
        spin.finish_with_message("Finished converting phylogenetic trees!\n");

        self.print_output_info(output, total_trees);

        let written = snapshot.written_files(output);
        let counts = TaskCounts {
            files_read: sorted.len(),
            files_written: written.len(),
            items_processed: Some(total_trees),
            items_skipped: None,
        };
        let report_file = report_path(output, None);
        builder
            .outputs(output, written.into_iter().map(|p| ("tree", p)))
            .statistics(&counts)?
            .write_to(&report_file)
    }

    fn convert_single(&self, input: &Path, output_dir: &Path) -> anyhow::Result<usize> {
        let output_path =
            output_tree_path(output_dir, input, *self.output_fmt, self.write_opts.gzip);
        progress::track(self.progress, input, || {
            let reader = TreeReader::open(input, *self.input_fmt, self.parse_opts.clone())?;
            let writer = TreeWriter::new(&output_path, *self.output_fmt, self.write_opts.clone());
            // Checked per tree, so a cancelled run stops mid-file.
            let trees = reader.map(|tree| {
                progress::check_cancelled(self.progress)?;
                tree
            });
            writer.write_stream(trees)
        })
    }

    fn print_output_info(&self, output: &Path, total_trees: usize) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Output dir", output.display());
        log::info!("{:18}: {}", "Output format", self.output_fmt);
        log::info!("{:18}: {}", "Trees written", utils::fmt_num(&total_trees));
    }
}

/// Build an output path in `output_dir`, stripping any trailing `.gz`/`.gzip`
/// from `input_file` before applying `output_fmt`'s extension (and re-adding
/// `.gz` when `gzip` is enabled).
pub fn output_tree_path(
    output_dir: &Path,
    input_file: &Path,
    output_fmt: TreeOutputFmt,
    gzip: bool,
) -> PathBuf {
    let name = input_file
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("tree");
    let stripped = strip_gzip_suffix(name);
    let stem_path = Path::new(stripped).with_extension(output_fmt.extension());
    let mut final_name = stem_path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("tree")
        .to_string();
    if gzip {
        final_name.push_str(".gz");
    }
    output_dir.join(final_name)
}

fn strip_gzip_suffix(name: &str) -> &str {
    if name.len() > 3 && name[name.len() - 3..].eq_ignore_ascii_case(".gz") {
        &name[..name.len() - 3]
    } else if name.len() > 5 && name[name.len() - 5..].eq_ignore_ascii_case(".gzip") {
        &name[..name.len() - 5]
    } else {
        name
    }
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::Ordering;

    use tempfile::TempDir;

    use super::*;
    use crate::helper::progress::CancelAfter;

    const MULTI_NWK: &str = "../../tests/files/trees/multi.nwk";

    #[test]
    fn test_progress_reports_each_file() {
        let output = TempDir::with_prefix("tree-convert").unwrap();
        let files = vec![PathBuf::from(MULTI_NWK)];
        let reporter = CancelAfter::new(usize::MAX);
        let trees = TreeConversion::new(&TreeFmt::Auto, &TreeOutputFmt::Nexus)
            .with_progress(&reporter)
            .convert(&files, output.path())
            .unwrap();
        assert!(trees > 1);
        assert_eq!(reporter.finished.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn test_cancelled_mid_file() {
        let output = TempDir::with_prefix("tree-convert").unwrap();
        let files = vec![PathBuf::from(MULTI_NWK)];
        // The first check passes before the file starts, then the first tree cancels.
        let reporter = CancelAfter::new(1);
        let error = TreeConversion::new(&TreeFmt::Auto, &TreeOutputFmt::Newick)
            .with_progress(&reporter)
            .convert(&files, output.path())
            .unwrap_err();
        assert!(error.downcast_ref::<progress::Cancelled>().is_some());
        assert_eq!(reporter.finished.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn test_output_tree_path() {
        let out = Path::new("out");
        assert_eq!(
            output_tree_path(
                out,
                Path::new("dir/gene.nwk.gz"),
                TreeOutputFmt::Nexus,
                false
            ),
            Path::new("out").join("gene.nex")
        );
        assert_eq!(
            output_tree_path(out, Path::new("dir/gene.nex"), TreeOutputFmt::Newick, true),
            Path::new("out").join("gene.nwk.gz")
        );
        assert_eq!(
            output_tree_path(
                out,
                Path::new("dir/net.enwk"),
                TreeOutputFmt::ExtendedNewick,
                false
            ),
            Path::new("out").join("net.enwk")
        );
    }
}
