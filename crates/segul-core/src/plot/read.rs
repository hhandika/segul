//! Raw read summary plots.
//!
//! Sample charts read `default-read-summary.csv`. Per-position charts read
//! `read-pos-summary.zip` (or one CSV extracted from it), written by
//! `segul read summary --mode complete`.
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};

use super::PlotPreset;
use super::batch::PlotWriter;
use super::preset::{bar_spec, plot_presets, scatter_spec};
use super::spec::{ChartKind, PlotSpec, PlotStyle, Series};
use crate::core::read::summarize::ReadSummaryFiles;
use crate::helper::types::SummaryMode;
use crate::parser::summary_csv::SummaryCsv;

const BASES: [&str; 5] = ["G", "C", "A", "T", "N"];
const POSITION_BASES: [&str; 4] = ["G", "C", "A", "T"];
// Above this many samples, an all-sample quality chart is too crowded to read.
const MAX_OVERVIEW_SAMPLES: usize = 8;

plot_presets! {
    /// Charts available for raw read summaries.
    ReadPlot {
        ReadCount => ("read-count", "Number of reads per sample (default-read-summary.csv)"),
        MeanLength => ("mean-length", "Mean read length per sample (default-read-summary.csv)"),
        MeanQuality => ("mean-qscore", "Mean quality score per sample (default-read-summary.csv)"),
        GcContent => ("gc-content", "GC content per sample (default-read-summary.csv)"),
        BaseComposition => ("base-composition", "Base composition per sample (default-read-summary.csv)"),
        GcVsReads => ("gc-vs-reads", "GC content vs number of reads (default-read-summary.csv)"),
        PositionQuality => ("pos-qscore", "Quality score by read position (read-pos-summary.zip)"),
        PositionComposition => ("pos-composition", "Base composition by read position for one sample (read-pos-summary.zip)"),
    }
}

/// Summary file read by a raw read chart.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadSource {
    /// Per-sample summary (`default-read-summary.csv`).
    Sample,
    /// Per-position archive (`read-pos-summary.zip`).
    Position,
}

impl ReadPlot {
    pub fn source(&self) -> ReadSource {
        match self {
            Self::PositionQuality | Self::PositionComposition => ReadSource::Position,
            _ => ReadSource::Sample,
        }
    }
}

/// Write every chart that applies to a read summary. Returns the written files.
///
/// Minimal summaries only have read counts. Complete summaries also get
/// quality and base composition by position for each sample, named
/// `pos-qscore_<sample>` and `pos-composition_<sample>`, plus an
/// all-sample `pos-qscore` chart when there are few samples.
pub fn plot_read_summary(files: &ReadSummaryFiles, writer: &PlotWriter) -> Vec<PathBuf> {
    let presets: Vec<ReadPlot> = match files.mode {
        SummaryMode::Minimal => vec![ReadPlot::ReadCount],
        _ => ReadPlot::all()
            .iter()
            .filter(|p| p.source() == ReadSource::Sample)
            .copied()
            .collect(),
    };
    let plotter = ReadPlotter::new(&files.summary, None, writer.style());
    let mut written: Vec<PathBuf> = presets
        .iter()
        .filter_map(|preset| writer.write(preset.name(), plotter.build(preset)))
        .collect();
    if let Some(positions) = &files.positions {
        written.extend(plot_positions(positions, writer));
    }
    written
}

fn plot_positions(archive: &Path, writer: &PlotWriter) -> Vec<PathBuf> {
    let samples: Vec<String> = match SummaryCsv::from_zip(archive) {
        Ok(tables) => tables
            .into_iter()
            .filter(|(_, table)| !table.is_empty())
            .map(|(name, _)| name)
            .collect(),
        Err(error) => {
            log::warn!("Skipped per-position charts: {:#}", error);
            return Vec::new();
        }
    };
    let mut written = Vec::new();
    if samples.len() > 1 && samples.len() <= MAX_OVERVIEW_SAMPLES {
        let plotter = ReadPlotter::new(archive, None, writer.style());
        let preset = ReadPlot::PositionQuality;
        written.extend(writer.write(preset.name(), plotter.build(&preset)));
    }
    for sample in &samples {
        let plotter = ReadPlotter::new(archive, Some(sample), writer.style());
        for preset in [ReadPlot::PositionQuality, ReadPlot::PositionComposition] {
            let name = format!("{}_{}", preset.name(), sample);
            written.extend(writer.write(&name, plotter.build(&preset)));
        }
    }
    written
}

/// Builds raw read summary charts.
pub struct ReadPlotter<'a> {
    input: &'a Path,
    /// Sample to plot from a per-position archive. All samples when `None`.
    sample: Option<&'a str>,
    style: &'a PlotStyle,
}

impl<'a> ReadPlotter<'a> {
    pub fn new(input: &'a Path, sample: Option<&'a str>, style: &'a PlotStyle) -> Self {
        Self {
            input,
            sample,
            style,
        }
    }

    pub fn build(&self, preset: &ReadPlot) -> Result<PlotSpec> {
        match preset {
            ReadPlot::PositionQuality => return self.position_quality(),
            ReadPlot::PositionComposition => return self.position_composition(),
            _ => (),
        }
        let csv = SummaryCsv::from_path(self.input)?;
        if csv.is_empty() {
            bail!("No rows in {}", self.input.display());
        }
        let samples = csv.column_str("file_name")?;
        let spec = match preset {
            ReadPlot::ReadCount => {
                self.sample_bar(&csv, samples, "read_count", "Number of reads")?
            }
            ReadPlot::MeanLength => {
                self.sample_bar(&csv, samples, "mean_read_length", "Mean read length")?
            }
            ReadPlot::MeanQuality => {
                self.sample_bar(&csv, samples, "mean_QScore", "Mean quality score")?
            }
            ReadPlot::GcContent => self.sample_bar(&csv, samples, "GC_content", "GC content")?,
            ReadPlot::BaseComposition => self.base_composition(&csv, samples)?,
            ReadPlot::GcVsReads => scatter_spec(
                "GC content vs number of reads",
                ("Number of reads", "GC content"),
                csv.column_f64("read_count")?,
                csv.column_f64("GC_content")?,
                self.style,
            ),
            ReadPlot::PositionQuality | ReadPlot::PositionComposition => unreachable!(),
        };
        Ok(spec)
    }

    fn sample_bar(
        &self,
        csv: &SummaryCsv,
        samples: Vec<String>,
        column: &str,
        label: &str,
    ) -> Result<PlotSpec> {
        Ok(bar_spec(
            &format!("{} per sample", label),
            ("Sample", label),
            samples,
            vec![Series::categorical(label, csv.column_f64(column)?)],
            self.style,
        ))
    }

    fn base_composition(&self, csv: &SummaryCsv, samples: Vec<String>) -> Result<PlotSpec> {
        let totals = csv.column_f64("base_count")?;
        let mut spec = bar_spec(
            "Base composition per sample",
            ("Sample", "Proportion of bases"),
            samples,
            Vec::new(),
            self.style,
        );
        spec.kind = ChartKind::StackedBar;
        for base in BASES {
            let proportions = csv
                .column_f64(base)?
                .iter()
                .zip(&totals)
                .map(|(count, total)| count / total)
                .collect();
            spec.series.push(Series::categorical(base, proportions));
        }
        Ok(spec)
    }

    fn position_quality(&self) -> Result<PlotSpec> {
        let tables = self.position_tables()?;
        let mut spec = PlotSpec::new(
            "Quality score by read position",
            ChartKind::Line,
            "Position in read (bp)",
            "Quality score",
        );
        let single = tables.len() == 1;
        for (name, table) in &tables {
            let mut series = Series::xy(
                name,
                table.column_f64("index")?,
                table.column_f64("mean_QScore")?,
            );
            // A min-max band is only readable for one sample.
            if single {
                series.lower = Some(table.column_f64("min_QScore")?);
                series.upper = Some(table.column_f64("max_QScore")?);
            }
            spec.series.push(series);
        }
        spec.style = self.style.clone();
        Ok(spec)
    }

    fn position_composition(&self) -> Result<PlotSpec> {
        let tables = self.position_tables()?;
        if tables.len() > 1 {
            let names: Vec<&str> = tables.iter().map(|(n, _)| n.as_str()).collect();
            bail!(
                "Base composition plots one sample at a time. Choose one of: {}",
                names.join(", ")
            );
        }
        let (name, table) = &tables[0];
        let mut spec = PlotSpec::new(
            &format!("Base composition by read position: {}", name),
            ChartKind::Line,
            "Position in read (bp)",
            "Proportion of bases",
        );
        let positions = table.column_f64("index")?;
        for base in POSITION_BASES {
            let column = format!("proportion_{}", base);
            spec.series.push(Series::xy(
                base,
                positions.clone(),
                table.column_f64(&column)?,
            ));
        }
        spec.style = self.style.clone();
        Ok(spec)
    }

    /// Per-position tables with at least one row, filtered by sample.
    fn position_tables(&self) -> Result<Vec<(String, SummaryCsv)>> {
        let is_zip = self
            .input
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("zip"));
        let tables = if is_zip {
            SummaryCsv::from_zip(self.input)?
        } else {
            let name = self
                .input
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or_default()
                .to_string();
            vec![(name, SummaryCsv::from_path(self.input)?)]
        };
        let tables: Vec<(String, SummaryCsv)> =
            tables.into_iter().filter(|(_, t)| !t.is_empty()).collect();
        match self.sample {
            Some(sample) => {
                let names: Vec<String> = tables.iter().map(|(n, _)| n.clone()).collect();
                match tables.into_iter().find(|(n, _)| n == sample) {
                    Some(table) => Ok(vec![table]),
                    None => bail!(
                        "Sample '{}' not found. Available samples: {}",
                        sample,
                        names.join(", ")
                    ),
                }
            }
            None if tables.is_empty() => bail!("No reads in {}", self.input.display()),
            None => Ok(tables),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::plot::fonts::PlotFonts;
    use crate::plot::spec::PlotFormat;

    const SUMMARY: &str = "tests/files/plot/default-read-summary.csv";
    const POSITIONS: &str = "tests/files/plot/read-pos-summary.zip";

    fn build(input: &str, sample: Option<&str>, preset: ReadPlot) -> Result<PlotSpec> {
        let style = PlotStyle::default();
        ReadPlotter::new(Path::new(input), sample, &style).build(&preset)
    }

    #[test]
    fn test_sample_presets_build() {
        for preset in ReadPlot::all()
            .iter()
            .filter(|p| p.source() == ReadSource::Sample)
        {
            build(SUMMARY, None, *preset).unwrap().validate().unwrap();
        }
    }

    #[test]
    fn test_base_composition() {
        let spec = build(SUMMARY, None, ReadPlot::BaseComposition).unwrap();
        assert_eq!(spec.series.len(), 5);
        // Empty files have no bases.
        assert!(spec.series[0].y[0].is_nan());
        // read_1.fastq: 21 G out of 72 bases.
        assert_eq!(spec.series[0].y[2], 21.0 / 72.0);
    }

    #[test]
    fn test_position_quality_all_samples() {
        let spec = build(POSITIONS, None, ReadPlot::PositionQuality).unwrap();
        // Empty samples are skipped.
        assert_eq!(spec.series.len(), 2);
        assert!(spec.series[0].lower.is_none());
        spec.validate().unwrap();
    }

    #[test]
    fn test_position_quality_one_sample() {
        let spec = build(POSITIONS, Some("read_1"), ReadPlot::PositionQuality).unwrap();
        assert_eq!(spec.series.len(), 1);
        assert!(spec.series[0].upper.is_some());
    }

    #[test]
    fn test_position_composition_needs_sample() {
        assert!(build(POSITIONS, None, ReadPlot::PositionComposition).is_err());
        let spec = build(POSITIONS, Some("read_2"), ReadPlot::PositionComposition).unwrap();
        assert_eq!(spec.series.len(), 4);
        assert!(build(POSITIONS, Some("missing"), ReadPlot::PositionComposition).is_err());
    }

    fn plot_files(mode: SummaryMode, positions: Option<&str>) -> Vec<String> {
        let dir = tempfile::TempDir::with_prefix("plot").unwrap();
        let style = PlotStyle::default();
        let fonts = PlotFonts::new();
        let writer = PlotWriter::new(dir.path(), None, &PlotFormat::Svg, &style, &fonts);
        let files = ReadSummaryFiles {
            summary: PathBuf::from(SUMMARY),
            mode,
            positions: positions.map(PathBuf::from),
        };
        plot_read_summary(&files, &writer)
            .iter()
            .map(|p| p.file_name().unwrap().to_string_lossy().to_string())
            .collect()
    }

    #[test]
    fn test_plot_read_summary_minimal() {
        assert_eq!(
            plot_files(SummaryMode::Minimal, None),
            vec!["read-count.svg"]
        );
    }

    #[test]
    fn test_plot_read_summary_complete() {
        let names = plot_files(SummaryMode::Complete, Some(POSITIONS));
        assert!(names.contains(&"base-composition.svg".to_string()));
        assert!(names.contains(&"pos-qscore.svg".to_string()));
        assert!(names.contains(&"pos-qscore_read_1.svg".to_string()));
        assert!(names.contains(&"pos-composition_read_2.svg".to_string()));
        // Empty samples get no per-position charts.
        assert!(!names.iter().any(|n| n.contains("read2_1")));
    }
}
