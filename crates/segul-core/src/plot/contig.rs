//! Contig summary plots from `contig-summary.csv`, written by `segul contig summary`.
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};

use super::PlotPreset;
use super::batch::PlotWriter;
use super::preset::{bar_spec, plot_presets, scatter_spec};
use super::spec::{PlotSpec, PlotStyle, Series};
use crate::parser::summary_csv::SummaryCsv;

plot_presets! {
    /// Charts available for contig summaries.
    ContigPlot {
        Nstats => ("nstats", "N50, N75, and N90 per sample"),
        ContigCount => ("contig-count", "Number of contigs per sample"),
        LongContigs => ("long-contigs", "Contigs of at least 750, 1000, and 1500 bp per sample"),
        GcContent => ("gc-content", "GC content per sample"),
        CountVsN50 => ("count-vs-n50", "N50 vs number of contigs"),
    }
}

/// Write every contig summary chart. Returns the written files.
pub fn plot_contig_summary(summary: &Path, writer: &PlotWriter) -> Vec<PathBuf> {
    let plotter = ContigPlotter::new(summary, writer.style());
    ContigPlot::all()
        .iter()
        .filter_map(|preset| writer.write(preset.name(), plotter.build(preset)))
        .collect()
}

/// Builds contig summary charts.
pub struct ContigPlotter<'a> {
    input: &'a Path,
    style: &'a PlotStyle,
}

impl<'a> ContigPlotter<'a> {
    pub fn new(input: &'a Path, style: &'a PlotStyle) -> Self {
        Self { input, style }
    }

    pub fn build(&self, preset: &ContigPlot) -> Result<PlotSpec> {
        let csv = SummaryCsv::from_path(self.input)?;
        if csv.is_empty() {
            bail!("No rows in {}", self.input.display());
        }
        let samples = csv.column_str("contig_name")?;
        let spec = match preset {
            ContigPlot::Nstats => self.grouped_bar(
                &csv,
                samples,
                &[("N50", "N50"), ("N75", "N75"), ("N90", "N90")],
                "Contig N statistics",
                "Length (bp)",
            )?,
            ContigPlot::ContigCount => self.grouped_bar(
                &csv,
                samples,
                &[("contig_count", "Contigs")],
                "Number of contigs per sample",
                "Number of contigs",
            )?,
            ContigPlot::LongContigs => self.grouped_bar(
                &csv,
                samples,
                &[
                    ("contig_750", "≥750 bp"),
                    ("contig_1000", "≥1000 bp"),
                    ("contig_1500", "≥1500 bp"),
                ],
                "Long contigs per sample",
                "Number of contigs",
            )?,
            ContigPlot::GcContent => self.grouped_bar(
                &csv,
                samples,
                &[("GC_content", "GC content")],
                "GC content per sample",
                "GC content",
            )?,
            ContigPlot::CountVsN50 => scatter_spec(
                "N50 vs number of contigs",
                ("Number of contigs", "N50 (bp)"),
                csv.column_f64("contig_count")?,
                csv.column_f64("N50")?,
                self.style,
            ),
        };
        Ok(spec)
    }

    fn grouped_bar(
        &self,
        csv: &SummaryCsv,
        samples: Vec<String>,
        columns: &[(&str, &str)],
        title: &str,
        y_label: &str,
    ) -> Result<PlotSpec> {
        let series = columns
            .iter()
            .map(|(column, name)| Ok(Series::categorical(name, csv.column_f64(column)?)))
            .collect::<Result<Vec<Series>>>()?;
        Ok(bar_spec(
            title,
            ("Sample", y_label),
            samples,
            series,
            self.style,
        ))
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::plot::fonts::PlotFonts;
    use crate::plot::spec::PlotFormat;

    const SUMMARY: &str = "tests/files/plot/contig-summary.csv";

    #[test]
    fn test_all_presets_build() {
        let style = PlotStyle::default();
        let plotter = ContigPlotter::new(Path::new(SUMMARY), &style);
        for preset in ContigPlot::all() {
            plotter.build(preset).unwrap().validate().unwrap();
        }
    }

    #[test]
    fn test_plot_contig_summary() {
        let dir = tempfile::TempDir::with_prefix("plot").unwrap();
        let style = PlotStyle::default();
        let fonts = PlotFonts::new();
        let writer = PlotWriter::new(dir.path(), Some("run"), &PlotFormat::Svg, &style, &fonts);
        let written = plot_contig_summary(Path::new(SUMMARY), &writer);
        assert_eq!(written.len(), ContigPlot::all().len());
        assert!(dir.path().join("run_nstats.svg").is_file());
    }

    #[test]
    fn test_nstats() {
        let style = PlotStyle::default();
        let spec = ContigPlotter::new(Path::new(SUMMARY), &style)
            .build(&ContigPlot::Nstats)
            .unwrap();
        assert_eq!(spec.categories, vec!["contigs2", "contigs1"]);
        assert_eq!(spec.series.len(), 3);
        assert_eq!(spec.series[0].y, vec![316.0, 258.0]);
    }
}
