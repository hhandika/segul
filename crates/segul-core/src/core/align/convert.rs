//! Convert sequence format from one to another.
use std::path::{Path, PathBuf};

use colored::Colorize;
use rayon::prelude::*;

use crate::core::OutputPrint;
use crate::helper::sequence::SeqParser;
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix};
use crate::helper::{files, utils};
use crate::report::spec::{AlignConvertParams, TaskKind, canonical};
use crate::report::stats::TaskCounts;
use crate::report::{DirSnapshot, ReportContext, RunReport, RunReportBuilder, report_path};
use crate::writer::sequences::SeqWriter;

impl OutputPrint for AlignmentConversion<'_> {}

/// Convert sequence format
pub struct AlignmentConversion<'a> {
    input_fmt: &'a InputFmt,
    output_fmt: &'a OutputFmt,
    datatype: &'a DataType,
    sort: bool,
    report_context: ReportContext,
}

impl<'a> AlignmentConversion<'a> {
    pub fn new(
        input_fmt: &'a InputFmt,
        output_fmt: &'a OutputFmt,
        datatype: &'a DataType,
        sort: bool,
    ) -> Self {
        Self {
            input_fmt,
            output_fmt,
            datatype,
            sort,
            report_context: ReportContext::default(),
        }
    }

    /// Attach frontend provenance and SHA-256 checksum settings for `report.json`.
    pub fn with_report_context(mut self, context: ReportContext) -> Self {
        self.report_context = context;
        self
    }

    pub fn convert(&self, files: &[PathBuf], output: &Path) {
        self.convert_with_report(files, output)
            .expect("Failed converting alignments");
    }

    /// Convert alignments and write a unified `report.json`.
    pub fn convert_with_report(
        &self,
        files: &[PathBuf],
        output: &Path,
    ) -> anyhow::Result<RunReport> {
        let mut sorted = files.to_vec();
        alphanumeric_sort::sort_path_slice(&mut sorted);
        let params = AlignConvertParams {
            input_format: canonical::input_fmt(self.input_fmt).to_string(),
            datatype: canonical::datatype(self.datatype).to_string(),
            output_format: canonical::output_fmt(self.output_fmt).to_string(),
            sort: self.sort,
        };
        let snapshot = DirSnapshot::capture(output);
        let builder = RunReportBuilder::new(TaskKind::AlignConvert, &self.report_context)
            .primary_inputs(&sorted)
            .parameters(&params)?;
        let spin = utils::set_spinner();
        spin.set_message("Converting sequence format...");
        sorted.par_iter().for_each(|file| {
            let output_fname = files::create_output_fname(output, file, self.output_fmt);
            self.convert_any(file, &output_fname);
        });
        spin.finish_with_message("Finished converting sequence format!\n");
        self.print_output_info(output);
        let written = snapshot.written_files(output);
        let counts = TaskCounts {
            files_read: sorted.len(),
            files_written: written.len(),
            items_processed: Some(written.len()),
            items_skipped: None,
        };
        let report_file = report_path(output, None);
        builder
            .outputs(output, written.into_iter().map(|p| ("alignment", p)))
            .statistics(&counts)?
            .write_to(&report_file)
    }

    fn convert_any(&self, input: &Path, output: &Path) {
        if self.sort {
            self.convert_sorted(input, output);
        } else {
            self.convert_unsorted(input, output);
        }
    }

    fn convert_unsorted(&self, input: &Path, output: &Path) {
        let (matrix, header) = self.get_sequence(input);
        self.write_results(&matrix, header, output);
    }

    fn convert_sorted(&self, input: &Path, output: &Path) {
        let (mut matrix, header) = self.get_sequence(input);
        matrix.sort_keys();
        self.write_results(&matrix, header, output)
    }

    fn get_sequence(&self, input: &Path) -> (SeqMatrix, Header) {
        let seq = SeqParser::new(input, self.datatype);
        seq.parse(self.input_fmt)
    }

    fn write_results(&self, matrix: &SeqMatrix, header: Header, output: &Path) {
        let mut convert = SeqWriter::new(output, matrix, &header);
        convert
            .write_sequence(self.output_fmt)
            .expect("Failed writing output files");
    }

    fn print_output_info(&self, output: &Path) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Output dir", output.display());
        self.print_output_fmt(self.output_fmt);
    }
}
