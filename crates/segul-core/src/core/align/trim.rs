//! Remove sites given a threshold of missing data.

use std::{
    fs,
    path::{Path, PathBuf},
    sync::mpsc,
};

use colored::Colorize;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::{
    core::OutputPrint,
    helper::{
        files,
        sequence::SeqParser,
        types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix},
        utils,
    },
    report::{
        DirSnapshot, ReportContext, RunReport, RunReportBuilder,
        paths::relative_posix,
        report_path,
        spec::{AlignTrimCriterion, AlignTrimParams, TaskKind, canonical},
        stats::{TrimEntry, TrimStats},
    },
    stats::sequence::Sites,
    writer::sequences::SeqWriter,
};

pub enum TrimmingParameters {
    /// Trim based on a threshold of missing data
    MissingData(f64),
    /// Trim based on a threshold of Parsimony Informative sites (PIS)
    /// PIS is the number of sites that have at least two different states
    /// and at least two of the states have a minimum frequency of 2
    ParsInf(usize),
    /// No trimming parameters
    None,
}

/// Trim alignment data structure
pub struct AlignmentTrimming<'a> {
    /// List of input files
    pub input_files: &'a [PathBuf],
    /// Input format of alignment files
    pub input_fmt: &'a InputFmt,
    /// Data type of sequences
    pub datatype: &'a DataType,
    /// Output directory
    pub output_dir: &'a Path,
    /// Output format of sequence files
    pub output_fmt: &'a OutputFmt,
    /// Trimming parameters
    pub params: &'a TrimmingParameters,
    report_context: ReportContext,
}

impl Default for AlignmentTrimming<'_> {
    fn default() -> Self {
        Self {
            input_files: &[],
            input_fmt: &InputFmt::Fasta,
            datatype: &DataType::Dna,
            output_dir: Path::new(""),
            output_fmt: &OutputFmt::Fasta,
            params: &TrimmingParameters::MissingData(0.1),
            report_context: ReportContext::default(),
        }
    }
}

impl OutputPrint for AlignmentTrimming<'_> {}

impl<'a> AlignmentTrimming<'a> {
    /// Create a new AlignmentTrimming instance
    pub fn new(
        input_files: &'a [PathBuf],
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        output_dir: &'a Path,
        output_fmt: &'a OutputFmt,
        params: &'a TrimmingParameters,
    ) -> Self {
        Self {
            input_files,
            input_fmt,
            datatype,
            output_dir,
            output_fmt,
            params,
            report_context: ReportContext::default(),
        }
    }

    /// Attach frontend provenance and SHA-256 checksum settings for `report.json`.
    pub fn with_report_context(mut self, context: ReportContext) -> Self {
        self.report_context = context;
        self
    }

    /// Trim alignment files based on the given parameters
    pub fn trim(&self) {
        self.trim_with_report().expect("Failed trimming alignments");
    }

    /// Trim alignment files and write a unified `report.json`.
    pub fn trim_with_report(&self) -> anyhow::Result<RunReport> {
        let mut sorted = self.input_files.to_vec();
        alphanumeric_sort::sort_path_slice(&mut sorted);
        let criterion = match self.params {
            TrimmingParameters::MissingData(max_missing) => AlignTrimCriterion::MissingData {
                max_missing: *max_missing,
            },
            TrimmingParameters::ParsInf(min_pars_inf) => AlignTrimCriterion::ParsInf {
                min_pars_inf: *min_pars_inf,
            },
            TrimmingParameters::None => {
                log::warn!("No trimming parameters provided. Skipping trimming...");
                std::process::exit(1);
            }
        };
        let params = AlignTrimParams {
            input_format: canonical::input_fmt(self.input_fmt).to_string(),
            datatype: canonical::datatype(self.datatype).to_string(),
            output_format: canonical::output_fmt(self.output_fmt).to_string(),
            criterion,
        };
        let snapshot = DirSnapshot::capture(self.output_dir);
        let builder = RunReportBuilder::new(TaskKind::AlignTrim, &self.report_context)
            .primary_inputs(&sorted)
            .parameters(&params)?;
        let spinner = utils::set_spinner();
        spinner.set_message("Trimming alignment files...");
        spinner.set_message("Writing summary...");
        let mut summary = self.trim_sites();
        summary.sort_by(|a, b| {
            alphanumeric_sort::compare_path(&a.parent_path, &b.parent_path)
                .then_with(|| alphanumeric_sort::compare_str(&a.file_name, &b.file_name))
        });
        self.write_summary(&summary);
        spinner.finish_with_message("Alignment files trimmed successfully\n");
        self.print_output_info();

        let base_dir = crate::report::paths::common_base_dir(&sorted);
        let entries = summary
            .iter()
            .map(|s| TrimEntry {
                path: relative_posix(&s.parent_path.join(&s.file_name), &base_dir),
                initial_sites: s.site_count_before,
                final_sites: s.site_count_after,
                removed_sites: s.site_removed,
            })
            .collect::<Vec<_>>();
        let stats = TrimStats {
            files: entries.len(),
            entries,
        };
        let written = snapshot.written_files(self.output_dir);
        let output_pairs = written
            .into_iter()
            .map(|p| {
                let role = if p.file_name().and_then(|n| n.to_str()) == Some("trimming_summary.csv")
                {
                    "summary"
                } else {
                    "alignment"
                };
                (role, p)
            })
            .collect::<Vec<_>>();
        let report_file = report_path(self.output_dir, None);
        builder
            .outputs(self.output_dir, output_pairs)
            .statistics(&stats)?
            .write_to(&report_file)
    }

    fn trim_sites(&self) -> Vec<TrimmingSummary> {
        match self.params {
            TrimmingParameters::MissingData(threshold) => self.par_trim_missing_data(*threshold),
            TrimmingParameters::ParsInf(threshold) => self.par_trim_informative_sites(*threshold),
            TrimmingParameters::None => {
                log::warn!("No trimming parameters provided. Skipping trimming...");
                std::process::exit(1);
            }
        }
    }

    fn par_trim_informative_sites(&self, threshold: usize) -> Vec<TrimmingSummary> {
        let (tx, rx) = mpsc::channel();
        self.input_files.par_iter().for_each_with(tx, |tx, file| {
            let (matrix, header) = self.parse_alignment(file);
            let mut summary = TrimmingSummary::new(file);
            match self.remove_site_with_informative(&matrix, threshold) {
                Some(new_matrix) => {
                    let nchar = self.write_output(&new_matrix, file);
                    summary.add_summary(header.nchar, nchar);
                }
                None => {
                    summary.add_summary(header.nchar, 0);
                }
            }
            tx.send(summary).expect("Failed to send summary");
        });
        rx.iter().collect()
    }

    fn par_trim_missing_data(&self, threshold: f64) -> Vec<TrimmingSummary> {
        let (tx, rx) = mpsc::channel();
        self.input_files.par_iter().for_each_with(tx, |tx, file| {
            let mut summary = TrimmingSummary::new(file);
            let (matrix, header) = self.parse_alignment(file);
            match self.remove_site_with_missing_data(&matrix, threshold) {
                Some(new_matrix) => {
                    let nchar = self.write_output(&new_matrix, file);
                    summary.add_summary(header.nchar, nchar);
                }
                None => {
                    let nchar = self.write_output(&matrix, file);
                    summary.add_summary(header.nchar, nchar);
                }
            }
            tx.send(summary).expect("Failed to send summary");
        });

        rx.iter().collect()
    }

    fn parse_alignment(&self, file: &Path) -> (SeqMatrix, Header) {
        SeqParser::new(file, self.datatype).get_alignment(self.input_fmt)
    }

    fn remove_site_with_informative(
        &self,
        matrix: &SeqMatrix,
        threshold: usize,
    ) -> Option<SeqMatrix> {
        let site_pos = Sites::default()
            .get_site_with_pars_informative(matrix, self.datatype, threshold)
            .iter()
            .map(|(i, _)| *i)
            .collect::<Vec<usize>>();

        if site_pos.is_empty() {
            None
        } else {
            Some(self.generate_new_matrix(matrix, site_pos))
        }
    }

    fn remove_site_with_missing_data(
        &self,
        matrix: &SeqMatrix,
        threshold: f64,
    ) -> Option<SeqMatrix> {
        let site_pos = Sites::default().get_site_without_missing_data(matrix, threshold);
        if site_pos.is_empty() {
            None
        } else {
            Some(self.generate_new_matrix(matrix, site_pos))
        }
    }

    fn generate_new_matrix(&self, matrix: &SeqMatrix, mut site_pos: Vec<usize>) -> SeqMatrix {
        site_pos.sort_unstable();
        let mut new_matrix = SeqMatrix::with_capacity(matrix.len(), matrix.len() * site_pos.len());
        for (k, seq) in matrix.iter() {
            new_matrix.push_with(k.clone(), |buf| {
                for &pos in &site_pos {
                    if let Some(&b) = seq.get(pos) {
                        buf.push(b);
                    }
                }
            });
        }

        new_matrix
    }

    // Write output and return the nchar (number of sites) in NEXUS terms
    fn write_output(&self, matrix: &SeqMatrix, file: &Path) -> usize {
        let alignment_dir = self.output_dir.join("trimmed_alignments");
        let output_path = files::create_output_fname(&alignment_dir, file, self.output_fmt);
        let mut header = Header::new();
        header.from_seq_matrix(matrix, true);
        let mut writer = SeqWriter::new(&output_path, matrix, &header);
        writer
            .write_sequence(self.output_fmt)
            .expect("Failed to write output file");
        header.nchar
    }

    fn write_summary(&self, summary: &[TrimmingSummary]) {
        let output_path = self
            .output_dir
            .join("trimming_summary")
            .with_extension("csv");
        fs::create_dir_all(self.output_dir).expect("Failed to create output directory");
        let mut writer = csv::Writer::from_path(output_path).expect("Failed to create CSV writer");
        summary.iter().for_each(|s| {
            writer.serialize(s).expect("Failed to write summary");
        });
    }

    fn print_output_info(&self) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Directory", self.output_dir.display());
        self.print_output_fmt(self.output_fmt);
    }
}

#[derive(Debug, Serialize, Default, Deserialize)]
struct TrimmingSummary {
    /// Parent path
    parent_path: PathBuf,
    /// File name
    file_name: String,
    /// Number of sites before trimming
    site_count_before: usize,
    /// Number of sites after trimming
    site_count_after: usize,
    /// Number of sites removed
    site_removed: usize,
}

impl TrimmingSummary {
    /// Create a new TrimmingSummary instance
    fn new(path: &Path) -> Self {
        Self {
            parent_path: path
                .parent()
                .expect("Failed to get parent path")
                .to_path_buf(),
            file_name: path
                .file_name()
                .expect("Failed to get file name")
                .to_str()
                .expect("Failed to convert to string")
                .to_string(),
            site_count_before: 0,
            site_count_after: 0,
            site_removed: 0,
        }
    }

    /// Generate a summary of the trimming process
    fn add_summary(&mut self, before: usize, after: usize) {
        if before < after {
            panic!("Number of sites after trimming is greater than before");
        }
        if before == after {
            self.site_count_before = before;
            self.site_count_after = after;
            self.site_removed = 0;
        } else {
            self.site_count_before = before;
            self.site_count_after = after;
            self.site_removed = before - after;
        }
    }
}

#[cfg(test)]
mod tests {
    use tempfile::TempDir;

    use super::*;

    macro_rules! init_trimming {
        ($input:expr_2021, $output:expr_2021, $params:expr_2021) => {
            AlignmentTrimming::new(
                $input,
                &InputFmt::Auto,
                &DataType::Dna,
                &$output.path(),
                &OutputFmt::Fasta,
                $params,
            )
        };
    }

    const INPUT_PATH: &str = "tests/files/trimming.fas";

    #[test]
    fn test_trim_missing_data() {
        let input_files = vec![PathBuf::from(INPUT_PATH)];
        let output_dir = TempDir::with_prefix("test").expect("Failed to create temp dir");
        let params = TrimmingParameters::MissingData(0.4);
        let align_trim = init_trimming!(&input_files, output_dir, &params);
        let summary = align_trim.trim_sites();
        let (matrix, header) = align_trim.parse_alignment(&input_files[0]);
        let site = Sites::default();
        let site_missing = site.get_site_without_missing_data(&matrix, 0.6);
        let index_site = site.index_site_with_missing_data(&matrix);
        assert_eq!(index_site.len(), 8);
        assert_eq!(index_site.get(&0).unwrap().len(), 4);
        assert_eq!(summary.len(), 1);
        assert_eq!(summary[0].site_count_before, 8);
        assert_eq!(summary[0].site_count_after, 7);
        assert_eq!(matrix.len(), 4);
        assert_eq!(header.nchar, 8);
        assert_eq!(site_missing.len(), 7);
    }

    #[test]
    fn test_trimming_results() {
        let input_files = vec![PathBuf::from(INPUT_PATH)];
        let output_dir = TempDir::with_prefix("test").expect("Failed to create temp dir");
        let params = TrimmingParameters::MissingData(0.4);
        let align_trim = init_trimming!(&input_files, output_dir, &params);
        align_trim.trim();
        let output_files = output_dir.path().read_dir().unwrap();
        assert_eq!(output_files.count(), 3);
    }
}
