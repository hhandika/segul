//! Generate alignment summary statistics
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;

use colored::Colorize;
use indexmap::IndexSet;
use rayon::prelude::*;

use crate::helper::finder::IDs;
use crate::helper::sequence::SeqParser;
use crate::helper::types::{DataType, InputFmt, TaxonRecords};
use crate::helper::utils;
use crate::report::spec::{AlignSummaryParams, TaskKind, canonical};
use crate::report::stats::TaskCounts;
use crate::report::{DirSnapshot, ReportContext, RunReport, RunReportBuilder, report_path};
use crate::stats::align_report::AlignReport;
use crate::stats::sequence::{CharMatrix, CharSummary, Completeness, SiteSummary, Sites, Taxa};
use crate::writer::summary::{CsvWriter, SummaryWriter};

/// Files written by [`AlignmentSummary::summarize_all`].
#[derive(Debug, Clone)]
pub struct AlignSummaryFiles {
    /// Per-locus statistics (`locus_summary.csv`).
    pub locus: PathBuf,
    /// Per-taxon statistics (`taxon_summary.csv`).
    pub taxon: PathBuf,
    /// Text summary (`alignment_summary.txt`).
    pub text: PathBuf,
    /// Unified provenance and statistics JSON (`report.json`).
    pub report: PathBuf,
}

pub struct AlignmentSummary<'a> {
    input_fmt: &'a InputFmt,
    output: &'a Path,
    datatype: &'a DataType,
    ntax: usize,
    interval: usize,
    report_context: ReportContext,
}

impl<'a> AlignmentSummary<'a> {
    pub fn new(
        input_fmt: &'a InputFmt,
        output: &'a Path,
        interval: usize,
        datatype: &'a DataType,
    ) -> Self {
        Self {
            input_fmt,
            output,
            ntax: 0,
            interval,
            datatype,
            report_context: ReportContext::default(),
        }
    }

    /// Attach frontend provenance and SHA-256 checksum settings for `report.json`.
    pub fn with_report_context(mut self, context: ReportContext) -> Self {
        self.report_context = context;
        self
    }

    /// Write the text, taxon, and locus summaries. Returns the written files.
    pub fn summarize_all(&mut self, files: &[PathBuf], prefix: Option<&str>) -> AlignSummaryFiles {
        self.check_datatype();
        let params = AlignSummaryParams {
            input_format: canonical::input_fmt(self.input_fmt).to_string(),
            datatype: canonical::datatype(self.datatype).to_string(),
            interval: self.interval,
            per_locus: false,
            prefix: prefix.map(str::to_string),
        };
        let mut builder = RunReportBuilder::new(TaskKind::AlignSummary, &self.report_context)
            .primary_inputs(files)
            .parameters(&params)
            .expect("Failed serializing alignment summary parameters");
        let spin = utils::set_spinner();
        spin.set_message("Indexing alignments...");
        let ids = self.get_id(files);
        self.ntax = ids.len();
        spin.set_message("Computing alignment stats...");
        let mut stats: Vec<(Sites, CharMatrix, Taxa)> = self.par_get_stats(files);
        stats.sort_by(|a, b| alphanumeric_sort::compare_path(&a.0.path, &b.0.path));
        let (sites, dna, complete) = self.summarize_char_matrix(&stats);
        let taxon_records = self.summarize_taxa(&ids, &stats);
        spin.finish_with_message("Finished computing summary stats!\n");
        let sum = SummaryWriter::new(&sites, &dna, &complete, self.datatype);
        let text = sum
            .write(self.output, prefix)
            .expect("Failed writing to stdout");
        let csv = CsvWriter::new(self.output, prefix, self.datatype);
        let taxon = csv
            .write_taxon_summary(&taxon_records)
            .expect("Failed writing to csv");
        let locus = csv
            .write_locus_summary(&stats)
            .expect("Failed writing a per locus csv file");
        let report = report_path(self.output, prefix);
        let align_report = AlignReport::new(
            &sites,
            &dna,
            &complete,
            &stats,
            &taxon_records,
            self.datatype,
        );
        builder = builder
            .outputs(
                self.output,
                [
                    ("locus_summary", locus.clone()),
                    ("taxon_summary", taxon.clone()),
                    ("text_summary", text.clone()),
                ],
            )
            .statistics(&align_report)
            .expect("Failed serializing alignment report");
        builder
            .write_to(&report)
            .expect("Failed writing the alignment report");
        log::info!("{:18}: {}", "Report", report.display());
        AlignSummaryFiles {
            locus,
            taxon,
            text,
            report,
        }
    }

    pub fn summarize_locus(&mut self, files: &[PathBuf], prefix: Option<&str>) {
        self.summarize_locus_with_report(files, prefix)
            .expect("Failed writing per-locus summary report");
    }

    /// Write per-locus summaries and record a unified `report.json`.
    pub fn summarize_locus_with_report(
        &mut self,
        files: &[PathBuf],
        prefix: Option<&str>,
    ) -> anyhow::Result<RunReport> {
        self.check_datatype();
        let params = AlignSummaryParams {
            input_format: canonical::input_fmt(self.input_fmt).to_string(),
            datatype: canonical::datatype(self.datatype).to_string(),
            interval: self.interval,
            per_locus: true,
            prefix: prefix.map(str::to_string),
        };
        let snapshot = DirSnapshot::capture(self.output);
        let builder = RunReportBuilder::new(TaskKind::AlignSummary, &self.report_context)
            .primary_inputs(files)
            .parameters(&params)?;
        let spin = utils::set_spinner();
        spin.set_message("Computing per locus summary...");
        files.iter().for_each(|file| {
            let (matrix, _) = SeqParser::new(file, self.datatype).get_alignment(self.input_fmt);
            let mut taxa = Taxa::new();
            taxa.summarize_taxa(&matrix, self.datatype);
            let csv = CsvWriter::new(self.output, prefix, self.datatype);
            csv.write_per_locus_summary(file, &taxa)
                .expect("Failed writing a taxon stats file");
        });
        spin.finish_with_message("Finished computing per locus summary!\n");
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Output dir", self.output.display());
        let written = snapshot.written_files(self.output);
        let counts = TaskCounts {
            files_read: files.len(),
            files_written: written.len(),
            items_processed: Some(files.len()),
            items_skipped: None,
        };
        let report = report_path(self.output, prefix);
        let run_report = builder
            .outputs(
                self.output,
                written.into_iter().map(|p| ("per_locus_summary", p)),
            )
            .statistics(&counts)?
            .write_to(&report)?;
        log::info!("{:18}: {}", "Report", report.display());
        Ok(run_report)
    }

    fn summarize_taxa(
        &self,
        ids: &IndexSet<String>,
        stats: &[(Sites, CharMatrix, Taxa)],
    ) -> BTreeMap<String, TaxonRecords> {
        let mut taxon_summary: BTreeMap<String, TaxonRecords> = BTreeMap::new();
        stats.iter().for_each(|(_, _, taxa)| {
            ids.iter().for_each(|id| {
                if let Some(char_counts) = taxa.records.get(id) {
                    match taxon_summary.get_mut(id) {
                        Some(taxon) => {
                            char_counts.chars.iter().for_each(|(c, count)| {
                                *taxon.char_counts.entry(*c).or_insert(0) += count;
                            });
                            taxon.locus_counts += 1;
                            taxon.total_chars += char_counts.total_chars;
                            taxon.missing_data += char_counts.missing_data;
                            if DataType::Dna == *self.datatype {
                                taxon.gc_count += char_counts.gc_count;
                                taxon.at_count += char_counts.at_count;
                                taxon.nucleotides += char_counts.nucleotides;
                            }
                        }
                        None => {
                            let mut taxon = TaxonRecords::new();
                            taxon.char_counts = char_counts.chars.clone();
                            taxon.locus_counts = 1;
                            taxon.total_chars = char_counts.total_chars;
                            taxon.missing_data = char_counts.missing_data;
                            if DataType::Dna == *self.datatype {
                                taxon.gc_count = char_counts.gc_count;
                                taxon.at_count = char_counts.at_count;
                                taxon.nucleotides = char_counts.nucleotides;
                            }
                            taxon_summary.insert(id.to_string(), taxon);
                        }
                    }
                }
            });
        });

        taxon_summary
    }

    fn get_id(&mut self, files: &[PathBuf]) -> IndexSet<String> {
        IDs::new(files, self.input_fmt, self.datatype).id_unique()
    }

    fn par_get_stats(&self, files: &[PathBuf]) -> Vec<(Sites, CharMatrix, Taxa)> {
        let (send, rec) = channel();
        files.par_iter().for_each_with(send, |s, file| {
            s.send(self.get_stats(file)).unwrap();
        });
        rec.iter().collect()
    }

    fn check_datatype(&mut self) {
        if let DataType::Ignore = self.datatype {
            self.datatype = &DataType::Dna
        }
    }

    fn get_stats(&self, path: &Path) -> (Sites, CharMatrix, Taxa) {
        let aln = SeqParser::new(path, self.datatype);
        let (matrix, header) = aln.get_alignment(self.input_fmt);
        let mut taxa = Taxa::new();
        taxa.summarize_taxa(&matrix, self.datatype);
        let mut dna = CharMatrix::new();
        dna.count_from_taxa(&taxa, &header, self.datatype);
        let mut sites = Sites::new(path);
        sites.get_stats(&matrix, self.datatype);

        (sites, dna, taxa)
    }

    fn summarize_char_matrix(
        &self,
        stats: &[(Sites, CharMatrix, Taxa)],
    ) -> (SiteSummary, CharSummary, Completeness) {
        let (sites, dna): (Vec<Sites>, Vec<CharMatrix>) =
            stats.par_iter().map(|p| (p.0.clone(), p.1.clone())).unzip();
        let mut sum_sites = SiteSummary::new();
        sum_sites.summarize(&sites);
        let mut sum_dna = CharSummary::new();
        sum_dna.summarize(&dna, self.datatype);
        let mut mat_comp = Completeness::new(&self.ntax, self.interval);
        mat_comp.matrix_completeness(&dna);
        (sum_sites, sum_dna, mat_comp)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_summary_report() {
        let output = TempDir::with_prefix("test").expect("Failed to create temp dir");
        let files: Vec<PathBuf> = ["gene_1.nex", "gene_2.nex", "gene_3.nex"]
            .iter()
            .map(|f| PathBuf::from("tests/files/alignments").join(f))
            .collect();
        let mut summary = AlignmentSummary::new(&InputFmt::Nexus, output.path(), 5, &DataType::Dna);
        let written = summary.summarize_all(&files, Some("test"));
        assert_eq!(written.report, output.path().join("test_report.json"));
        let run = RunReport::from_path(&written.report).unwrap();
        assert_eq!(run.schema_version, 2);
        assert_eq!(run.kind, TaskKind::AlignSummary);
        assert_eq!(run.inputs.primary.len(), 3);
        assert_eq!(run.outputs.files.len(), 3);
        let report = AlignReport::from_path(&written.report).unwrap();
        assert_eq!(report.datatype, "dna");
        assert_eq!(report.summary.loci, 3);
        assert_eq!(report.loci.len(), 3);
        assert_eq!(report.loci[0].name, "gene_1");
        assert_eq!(report.taxa.len(), report.summary.taxa);
        let sites: usize = report.loci.iter().map(|l| l.sites).sum();
        assert_eq!(sites, report.summary.sites);
        let missing: usize = report.loci.iter().map(|l| l.missing).sum();
        assert_eq!(missing, report.summary.missing);
        assert!(report.summary.gc.is_some());
        assert_eq!(report.completeness[0].percent, 100);
        // Taxa share the characters of the loci between them.
        let chars: usize = report.taxa.iter().map(|t| t.chars).sum();
        assert_eq!(chars, report.summary.chars);
    }
}
