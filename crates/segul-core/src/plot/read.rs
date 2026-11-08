//! Raw read summary plots.
//!
//! Sample charts read `default-read-summary.csv`. Per-position charts read
//! `read-pos-summary.zip` (or one CSV extracted from it). Distribution charts
//! read `read-report.json`. The last two are written by
//! `segul read summary --mode complete`.
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};

use super::PlotPreset;
use super::batch::PlotWriter;
use super::preset::{bar_spec, plot_presets, scatter_spec};
use super::spec::{BoxStats, ChartKind, PlotSpec, PlotStyle, Series};
use crate::core::read::summarize::ReadSummaryFiles;
use crate::helper::types::SummaryMode;
use crate::parser::summary_csv::SummaryCsv;
use crate::stats::fastq::report::{ReadReport, SampleReadReport};

const BASES: [&str; 5] = ["G", "C", "A", "T", "N"];
const POSITION_BASES: [&str; 4] = ["G", "C", "A", "T"];
// Above this many samples, an all-sample quality chart is too crowded to read.
const MAX_OVERVIEW_SAMPLES: usize = 8;
// Box plots with more positions than this show every n-th position.
const MAX_POSITION_BOXES: usize = 75;

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
        Q30Percent => ("q30-pct", "Bases with quality 30 or higher per sample (default-read-summary.csv)"),
        PositionQualityBox => ("pos-qscore-box", "Quality score distribution by read position for one sample (read-report.json)"),
        LengthDistribution => ("length-dist", "Read length distribution (read-report.json)"),
        ReadQualityDistribution => ("read-qscore-dist", "Mean quality score per read distribution (read-report.json)"),
        GcDistribution => ("gc-dist", "GC content per read with the theoretical distribution for one sample (read-report.json)"),
        PositionNContent => ("pos-n-content", "N content by read position (read-report.json)"),
        AdapterContent => ("adapter-content", "Cumulative adapter content by read position for one sample (read-report.json)"),
        Duplication => ("duplication", "Sequence duplication levels for one sample (read-report.json)"),
    }
}

/// Summary file read by a raw read chart.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadSource {
    /// Per-sample summary (`default-read-summary.csv`).
    Sample,
    /// Per-position archive (`read-pos-summary.zip`).
    Position,
    /// Complete statistics (`read-report.json`).
    Report,
}

impl ReadPlot {
    /// Whether the chart shows one sample at a time.
    pub fn is_per_sample(&self) -> bool {
        matches!(
            self,
            Self::PositionComposition
                | Self::PositionQualityBox
                | Self::GcDistribution
                | Self::AdapterContent
                | Self::Duplication
        )
    }

    pub fn source(&self) -> ReadSource {
        match self {
            Self::PositionQuality | Self::PositionComposition => ReadSource::Position,
            Self::PositionQualityBox
            | Self::LengthDistribution
            | Self::ReadQualityDistribution
            | Self::GcDistribution
            | Self::PositionNContent
            | Self::AdapterContent
            | Self::Duplication => ReadSource::Report,
            _ => ReadSource::Sample,
        }
    }
}

/// Write every chart that applies to a read summary. Returns the written files.
///
/// Minimal summaries only have read counts. Complete summaries also get
/// quality and base composition by position for each sample, named
/// `pos-qscore_<sample>` and `pos-composition_<sample>`, plus an
/// all-sample `pos-qscore` chart when there are few samples. When the run
/// wrote `read-report.json`, distribution charts are added the same way.
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
    if let Some(report) = &files.report {
        written.extend(plot_report(report, writer));
    }
    written
}

fn plot_report(path: &Path, writer: &PlotWriter) -> Vec<PathBuf> {
    let samples: Vec<String> = match ReadReport::from_path(path) {
        Ok(report) => report
            .samples
            .into_iter()
            .filter(|s| s.metrics.read_count > 0)
            .map(|s| s.sample)
            .collect(),
        Err(error) => {
            log::warn!("Skipped read distribution charts: {:#}", error);
            return Vec::new();
        }
    };
    let presets: Vec<ReadPlot> = ReadPlot::all()
        .iter()
        .filter(|p| p.source() == ReadSource::Report)
        .copied()
        .collect();
    let mut written = Vec::new();
    if !samples.is_empty() && samples.len() <= MAX_OVERVIEW_SAMPLES {
        let plotter = ReadPlotter::new(path, None, writer.style());
        for preset in presets.iter().filter(|p| !p.is_per_sample()) {
            written.extend(writer.write(preset.name(), plotter.build(preset)));
        }
    }
    for sample in &samples {
        let plotter = ReadPlotter::new(path, Some(sample), writer.style());
        for preset in presets.iter().filter(|p| p.is_per_sample()) {
            let name = format!("{}_{}", preset.name(), sample);
            written.extend(writer.write(&name, plotter.build(preset)));
        }
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
            preset if preset.source() == ReadSource::Report => return self.report_chart(preset),
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
            ReadPlot::Q30Percent => {
                self.sample_bar(&csv, samples, "Q30_pct", "Bases with Q30 or higher (%)")?
            }
            ReadPlot::BaseComposition => self.base_composition(&csv, samples)?,
            ReadPlot::GcVsReads => scatter_spec(
                "GC content vs number of reads",
                ("Number of reads", "GC content"),
                csv.column_f64("read_count")?,
                csv.column_f64("GC_content")?,
                self.style,
            ),
            _ => unreachable!("Chart reads another summary file"),
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

    fn report_chart(&self, preset: &ReadPlot) -> Result<PlotSpec> {
        let samples = self.report_samples(preset.is_per_sample())?;
        let mut spec = match preset {
            ReadPlot::PositionQualityBox => position_quality_box(&samples[0]),
            ReadPlot::LengthDistribution => overlay(
                "Read length distribution",
                ("Read length (bp)", "Number of reads"),
                &samples,
                |s| {
                    s.length_distribution
                        .iter()
                        .map(|l| (l.length as f64, l.reads as f64))
                        .collect()
                },
            ),
            ReadPlot::ReadQualityDistribution => overlay(
                "Mean quality score per read",
                ("Mean quality score", "Number of reads"),
                &samples,
                |s| indexed(&s.read_quality_distribution),
            ),
            ReadPlot::GcDistribution => gc_distribution(&samples[0]),
            ReadPlot::PositionNContent => overlay(
                "N content by read position",
                ("Position in read (bp)", "N content (%)"),
                &samples,
                |s| {
                    s.positions
                        .iter()
                        .map(|p| {
                            let share = p.n_count as f64 / p.base_count().max(1) as f64;
                            (p.start as f64, 100.0 * share)
                        })
                        .collect()
                },
            ),
            ReadPlot::AdapterContent => adapter_content(&samples[0])?,
            ReadPlot::Duplication => duplication(&samples[0])?,
            _ => unreachable!("Chart reads another summary file"),
        };
        spec.style = self.style.clone();
        Ok(spec)
    }

    /// Samples with reads in `read-report.json`, filtered by sample.
    fn report_samples(&self, single: bool) -> Result<Vec<SampleReadReport>> {
        let report = ReadReport::from_path(self.input)?;
        let names: Vec<String> = report.samples.iter().map(|s| s.sample.clone()).collect();
        let samples: Vec<SampleReadReport> = match self.sample {
            Some(sample) => match report.samples.into_iter().find(|s| s.sample == sample) {
                Some(found) => vec![found],
                None => bail!(
                    "Sample '{}' not found. Available samples: {}",
                    sample,
                    names.join(", ")
                ),
            },
            None => report
                .samples
                .into_iter()
                .filter(|s| s.metrics.read_count > 0)
                .collect(),
        };
        if samples.is_empty() {
            bail!("No reads in {}", self.input.display());
        }
        if single && samples.len() > 1 {
            bail!(
                "This chart plots one sample at a time. Choose one of: {}",
                names.join(", ")
            );
        }
        Ok(samples)
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

/// Line chart with one series per sample.
fn overlay(
    title: &str,
    labels: (&str, &str),
    samples: &[SampleReadReport],
    points: impl Fn(&SampleReadReport) -> Vec<(f64, f64)>,
) -> PlotSpec {
    let mut spec = PlotSpec::new(title, ChartKind::Line, labels.0, labels.1);
    spec.series = samples
        .iter()
        .map(|s| {
            let (x, y) = points(s).into_iter().unzip();
            Series::xy(&s.sample, x, y)
        })
        .collect();
    spec
}

fn indexed(counts: &[u64]) -> Vec<(f64, f64)> {
    counts
        .iter()
        .enumerate()
        .map(|(i, c)| (i as f64, *c as f64))
        .collect()
}

/// Box per position with 10th and 90th percentile whiskers, like FastQC.
fn position_quality_box(sample: &SampleReadReport) -> PlotSpec {
    let mut spec = PlotSpec::new(
        &format!("Quality score by read position: {}", sample.sample),
        ChartKind::BoxPlot,
        "Position in read (bp)",
        "Quality score",
    );
    let step = sample.positions.len().div_ceil(MAX_POSITION_BOXES).max(1);
    spec.boxes = sample
        .positions
        .iter()
        .step_by(step)
        .map(|p| BoxStats {
            label: p.start.to_string(),
            lower_whisker: p.p10_qscore as f64,
            q1: p.q1_qscore as f64,
            median: p.median_qscore as f64,
            q3: p.q3_qscore as f64,
            upper_whisker: p.p90_qscore as f64,
            outliers: Vec::new(),
        })
        .collect();
    spec
}

fn gc_distribution(sample: &SampleReadReport) -> PlotSpec {
    let mut spec = PlotSpec::new(
        &format!("GC content per read: {}", sample.sample),
        ChartKind::Line,
        "GC content (%)",
        "Number of reads",
    );
    let x: Vec<f64> = (0..sample.gc_distribution.len())
        .map(|i| i as f64)
        .collect();
    let y = sample.gc_distribution.clone();
    spec.series.push(Series::xy("Observed", x.clone(), y));
    spec.series
        .push(Series::xy("Theoretical", x, sample.gc_theoretical.clone()));
    spec
}

fn adapter_content(sample: &SampleReadReport) -> Result<PlotSpec> {
    if sample.adapters.is_empty() {
        bail!("No adapter content for {}", sample.sample);
    }
    let mut spec = PlotSpec::new(
        &format!("Adapter content: {}", sample.sample),
        ChartKind::Line,
        "Position in read (bp)",
        "Reads with adapter (%)",
    );
    let positions: Vec<f64> = sample.positions.iter().map(|p| p.start as f64).collect();
    spec.series = sample
        .adapters
        .iter()
        .map(|a| Series::xy(&a.name, positions.clone(), a.cumulative_pct.clone()))
        .collect();
    Ok(spec)
}

fn duplication(sample: &SampleReadReport) -> Result<PlotSpec> {
    let Some(duplication) = &sample.duplication else {
        bail!("No duplication estimate for {}", sample.sample);
    };
    let estimate = if duplication.exact { "" } else { ", estimated" };
    let spec = bar_spec(
        &format!(
            "Sequence duplication: {} ({:.1}% duplicate{})",
            sample.sample, duplication.duplicate_pct, estimate
        ),
        ("Copies of a sequence", "Share (%)"),
        duplication.levels.iter().map(|l| l.label.clone()).collect(),
        vec![
            Series::categorical(
                "Distinct sequences",
                duplication.levels.iter().map(|l| l.distinct_pct).collect(),
            ),
            Series::categorical(
                "Reads",
                duplication.levels.iter().map(|l| l.reads_pct).collect(),
            ),
        ],
        &PlotStyle::default(),
    );
    Ok(spec)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::plot::fonts::PlotFonts;
    use crate::plot::spec::PlotFormat;

    const SUMMARY: &str = "tests/files/plot/default-read-summary.csv";
    const POSITIONS: &str = "tests/files/plot/read-pos-summary.zip";
    const REPORT: &str = "tests/files/plot/read-report.json";

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
            distributions: None,
            report: positions.map(|_| PathBuf::from(REPORT)),
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
        assert!(names.contains(&"length-dist.svg".to_string()));
        assert!(names.contains(&"pos-qscore-box_read_1.svg".to_string()));
        assert!(names.contains(&"duplication_read_2.svg".to_string()));
    }

    #[test]
    fn test_report_presets_build() {
        for preset in ReadPlot::all()
            .iter()
            .filter(|p| p.source() == ReadSource::Report)
        {
            let sample = preset.is_per_sample().then_some("read_1");
            build(REPORT, sample, *preset).unwrap().validate().unwrap();
        }
    }

    #[test]
    fn test_report_per_sample_needs_sample() {
        assert!(build(REPORT, None, ReadPlot::PositionQualityBox).is_err());
        let spec = build(REPORT, Some("read_1"), ReadPlot::PositionQualityBox).unwrap();
        assert_eq!(spec.boxes.len(), 36);
        let spec = build(REPORT, Some("read_1"), ReadPlot::GcDistribution).unwrap();
        assert_eq!(spec.series.len(), 2);
    }
}
