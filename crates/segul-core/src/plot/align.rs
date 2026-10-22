//! Alignment summary plots.
//!
//! Locus charts read `locus_summary.csv` and taxon charts read
//! `taxon_summary.csv`, both written by `segul align summary`.
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};

use super::PlotPreset;
use super::batch::PlotWriter;
use super::preset::{bar_spec, histogram_spec, plot_presets, scatter_spec};
use super::spec::{BoxStats, ChartKind, PlotSpec, PlotStyle, Series};
use crate::core::align::summarize::AlignSummaryFiles;
use crate::parser::summary_csv::SummaryCsv;

const LOCUS_SUMMARY: &str = "locus_summary";
const TAXON_SUMMARY: &str = "taxon_summary";
// Matches the default interval of the text summary.
const COMPLETENESS_INTERVAL: usize = 5;

plot_presets! {
    /// Charts available for alignment summaries.
    AlignPlot {
        SiteHistogram => ("site-hist", "Histogram of alignment length (locus_summary.csv)"),
        PisHistogram => ("pis-hist", "Histogram of parsimony informative sites (locus_summary.csv)"),
        MissingHistogram => ("missing-hist", "Histogram of proportion of missing data (locus_summary.csv)"),
        PisVsLength => ("pis-vs-length", "Parsimony informative sites vs alignment length (locus_summary.csv)"),
        VarVsLength => ("var-vs-length", "Variable sites vs alignment length (locus_summary.csv)"),
        SiteProportions => ("site-proportions", "Box plot of site proportions across loci (locus_summary.csv)"),
        Completeness => ("completeness", "Loci retained by minimum taxon percentage (locus_summary.csv)"),
        TaxonLoci => ("taxon-loci", "Number of loci per taxon (taxon_summary.csv)"),
        TaxonMissing => ("taxon-missing", "Proportion of missing data per taxon (taxon_summary.csv)"),
        TaxonGc => ("taxon-gc", "GC content per taxon, DNA only (taxon_summary.csv)"),
    }
}

/// Summary file read by an alignment chart.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignSource {
    Locus,
    Taxon,
}

impl AlignPlot {
    pub fn source(&self) -> AlignSource {
        match self {
            Self::TaxonLoci | Self::TaxonMissing | Self::TaxonGc => AlignSource::Taxon,
            _ => AlignSource::Locus,
        }
    }
}

/// Write every alignment summary chart. Returns the written files.
/// GC content is skipped for amino acid summaries, which have no GC column.
pub fn plot_alignment_summary(files: &AlignSummaryFiles, writer: &PlotWriter) -> Vec<PathBuf> {
    let has_gc = SummaryCsv::from_path(&files.taxon)
        .map(|csv| csv.has_column("gc_content"))
        .unwrap_or(false);
    AlignPlot::all()
        .iter()
        .filter(|preset| has_gc || **preset != AlignPlot::TaxonGc)
        .filter_map(|preset| {
            let input = match preset.source() {
                AlignSource::Locus => &files.locus,
                AlignSource::Taxon => &files.taxon,
            };
            let spec = AlignmentPlot::new(input, writer.style()).build(preset);
            writer.write(preset.name(), spec)
        })
        .collect()
}

/// Builds alignment summary charts.
pub struct AlignmentPlot<'a> {
    input: &'a Path,
    style: &'a PlotStyle,
}

impl<'a> AlignmentPlot<'a> {
    pub fn new(input: &'a Path, style: &'a PlotStyle) -> Self {
        Self { input, style }
    }

    pub fn build(&self, preset: &AlignPlot) -> Result<PlotSpec> {
        let csv = SummaryCsv::from_path(self.input)?;
        if csv.is_empty() {
            bail!("No rows in {}", self.input.display());
        }
        let spec = match preset {
            AlignPlot::SiteHistogram => histogram_spec(
                "Alignment length",
                "Alignment length (sites)",
                &csv.column_f64("site_count")?,
                self.style,
            ),
            AlignPlot::PisHistogram => histogram_spec(
                "Parsimony informative sites",
                "Parsimony informative sites",
                &csv.column_f64("parsimony_informative_sites")?,
                self.style,
            ),
            AlignPlot::MissingHistogram => histogram_spec(
                "Missing data",
                "Proportion of missing data",
                &csv.column_f64("proportion_missing_data")?,
                self.style,
            ),
            AlignPlot::PisVsLength => scatter_spec(
                "Parsimony informative sites vs alignment length",
                ("Alignment length (sites)", "Parsimony informative sites"),
                csv.column_f64("site_count")?,
                csv.column_f64("parsimony_informative_sites")?,
                self.style,
            ),
            AlignPlot::VarVsLength => scatter_spec(
                "Variable sites vs alignment length",
                ("Alignment length (sites)", "Variable sites"),
                csv.column_f64("site_count")?,
                csv.column_f64("variable_sites")?,
                self.style,
            ),
            AlignPlot::SiteProportions => self.site_proportions(&csv)?,
            AlignPlot::Completeness => self.completeness(&csv)?,
            AlignPlot::TaxonLoci => {
                self.taxon_bar(&csv, "locus_counts", "Loci per taxon", "Number of loci")?
            }
            AlignPlot::TaxonMissing => self.taxon_bar(
                &csv,
                "proportion_missing_data",
                "Missing data per taxon",
                "Proportion of missing data",
            )?,
            AlignPlot::TaxonGc => {
                self.taxon_bar(&csv, "gc_content", "GC content per taxon", "GC content")?
            }
        };
        Ok(spec)
    }

    fn site_proportions(&self, csv: &SummaryCsv) -> Result<PlotSpec> {
        let columns = [
            ("proportion_cons_sites", "Conserved"),
            ("proportion_var_sites", "Variable"),
            ("proportion_pars_inf_sites", "Parsimony informative"),
            ("proportion_missing_data", "Missing data"),
        ];
        let mut spec = PlotSpec::new(
            "Site proportions across loci",
            ChartKind::BoxPlot,
            "",
            "Proportion",
        );
        for (column, label) in columns {
            if let Some(stats) = BoxStats::from_values(label, &csv.column_f64(column)?) {
                spec.boxes.push(stats);
            }
        }
        spec.style = self.style.clone();
        Ok(spec)
    }

    /// Matrix completeness: loci that contain at least a given percentage of taxa.
    fn completeness(&self, csv: &SummaryCsv) -> Result<PlotSpec> {
        let ntax = csv.column_f64("ntax")?;
        let total_taxa = self.total_taxa(&ntax);
        let percents: Vec<usize> = (1..=100 / COMPLETENESS_INTERVAL)
            .map(|i| i * COMPLETENESS_INTERVAL)
            .collect();
        let loci = percents
            .iter()
            .map(|p| {
                let min_taxa = (total_taxa as f64 * *p as f64 / 100.0).floor();
                ntax.iter().filter(|n| **n >= min_taxa).count() as f64
            })
            .collect();
        let mut spec = PlotSpec::new(
            "Matrix completeness",
            ChartKind::Line,
            "Minimum taxa (%)",
            "Number of loci",
        );
        let percents = percents.into_iter().map(|p| p as f64).collect();
        spec.series.push(Series::xy("Loci", percents, loci));
        spec.style = self.style.clone();
        Ok(spec)
    }

    /// Total taxa from the matching taxon summary, falling back to the
    /// largest locus when the taxon summary is not next to the input.
    fn total_taxa(&self, ntax: &[f64]) -> usize {
        let max_ntax = ntax
            .iter()
            .filter(|n| n.is_finite())
            .fold(0.0, |acc: f64, n| acc.max(*n)) as usize;
        match self.taxon_summary_path() {
            Some(path) => match SummaryCsv::from_path(&path) {
                Ok(taxa) => taxa.len().max(max_ntax),
                Err(_) => max_ntax,
            },
            None => {
                log::warn!(
                    "Taxon summary not found next to {}. \
                    Using the largest locus as the total number of taxa.",
                    self.input.display()
                );
                max_ntax
            }
        }
    }

    fn taxon_summary_path(&self) -> Option<PathBuf> {
        let file_name = self.input.file_name()?.to_str()?;
        if !file_name.contains(LOCUS_SUMMARY) {
            return None;
        }
        let path = self
            .input
            .with_file_name(file_name.replace(LOCUS_SUMMARY, TAXON_SUMMARY));
        path.is_file().then_some(path)
    }

    fn taxon_bar(
        &self,
        csv: &SummaryCsv,
        column: &str,
        title: &str,
        y_label: &str,
    ) -> Result<PlotSpec> {
        let taxa = csv.column_str("taxon")?;
        let values = csv.column_f64(column)?;
        Ok(bar_spec(
            title,
            ("Taxon", y_label),
            taxa,
            vec![Series::categorical(y_label, values)],
            self.style,
        ))
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::plot::fonts::PlotFonts;
    use crate::plot::spec::PlotFormat;

    const LOCUS: &str = "tests/files/plot/locus_summary.csv";
    const TAXON: &str = "tests/files/plot/taxon_summary.csv";

    fn build(input: &str, preset: AlignPlot) -> PlotSpec {
        let style = PlotStyle::default();
        AlignmentPlot::new(Path::new(input), &style)
            .build(&preset)
            .unwrap()
    }

    #[test]
    fn test_all_presets_build() {
        for preset in AlignPlot::all() {
            let input = match preset.source() {
                AlignSource::Taxon => TAXON,
                AlignSource::Locus => LOCUS,
            };
            let spec = build(input, *preset);
            spec.validate().unwrap();
        }
    }

    #[test]
    fn test_scatter_values() {
        let spec = build(LOCUS, AlignPlot::VarVsLength);
        assert_eq!(spec.series[0].x, vec![6.0, 8.0, 6.0, 6.0]);
        assert_eq!(spec.series[0].y, vec![0.0, 1.0, 1.0, 2.0]);
    }

    #[test]
    fn test_completeness_uses_taxon_summary() {
        // Three taxa in total, so 100% needs 3 taxa and no locus qualifies.
        let spec = build(LOCUS, AlignPlot::Completeness);
        let series = &spec.series[0];
        assert_eq!(series.x.last(), Some(&100.0));
        assert_eq!(series.y.last(), Some(&0.0));
        // 65% needs floor(1.95) = 1 taxon.
        let idx = series.x.iter().position(|x| *x == 65.0).unwrap();
        assert_eq!(series.y[idx], 4.0);
        // 70% needs floor(2.1) = 2 taxa.
        assert_eq!(series.y[idx + 1], 3.0);
    }

    #[test]
    fn test_taxon_bar() {
        let spec = build(TAXON, AlignPlot::TaxonLoci);
        assert_eq!(spec.categories, vec!["ABCD", "ABCE", "ABCF"]);
        assert_eq!(spec.series[0].y, vec![4.0, 2.0, 1.0]);
    }

    #[test]
    fn test_plot_alignment_summary() {
        let dir = tempfile::TempDir::with_prefix("plot").unwrap();
        let style = PlotStyle::default();
        let fonts = PlotFonts::new();
        let writer = PlotWriter::new(dir.path(), None, &PlotFormat::Svg, &style, &fonts);
        let files = AlignSummaryFiles {
            locus: PathBuf::from(LOCUS),
            taxon: PathBuf::from(TAXON),
            text: PathBuf::new(),
        };
        let written = plot_alignment_summary(&files, &writer);
        assert_eq!(written.len(), AlignPlot::all().len());
        assert!(dir.path().join("completeness.svg").is_file());
        assert!(dir.path().join("taxon-gc.svg").is_file());
    }

    #[test]
    fn test_wrong_input() {
        let style = PlotStyle::default();
        let result = AlignmentPlot::new(Path::new(TAXON), &style).build(&AlignPlot::PisHistogram);
        assert!(result.is_err());
    }
}
