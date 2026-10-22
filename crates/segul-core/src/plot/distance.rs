//! Pairwise distance plots.
//!
//! Reads the outputs of `segul align distance`: the per-locus summary
//! (`<prefix>_summary.csv`), a distance matrix (`per_locus/<locus>.csv` or
//! `<prefix>_group_mean.csv`), and the long table (`<prefix>_pairwise_long.csv`).
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};
use indexmap::IndexMap;

use super::PlotPreset;
use super::batch::PlotWriter;
use super::preset::{bar_spec, histogram_spec, plot_presets};
use super::spec::{BoxStats, ChartKind, HeatmapData, PlotSpec, PlotStyle, Series};
use crate::core::align::distance::DistanceFiles;
use crate::parser::summary_csv::SummaryCsv;

// Heatmaps are written for at most this many loci, in locus order.
const MAX_MATRIX_PLOTS: usize = 50;

plot_presets! {
    /// Charts available for pairwise distances.
    DistancePlot {
        LocusMean => ("locus-mean", "Mean pairwise distance per locus (distance_summary.csv)"),
        MeanHistogram => ("mean-hist", "Histogram of mean pairwise distance per locus (distance_summary.csv)"),
        Matrix => ("matrix", "Heatmap of a distance matrix (per_locus/<locus>.csv or group_mean.csv)"),
        PairwiseHistogram => ("pairwise-hist", "Histogram of all pairwise distances (pairwise_long.csv)"),
        LocusBox => ("locus-box", "Box plot of pairwise distances per locus (pairwise_long.csv)"),
    }
}

/// Write every chart that applies to the distance outputs. Returns the written files.
///
/// Writes the per-locus summary charts, one `matrix_<locus>` heatmap for each
/// locus with at least two samples (up to 50 loci), the `group-mean` heatmap
/// with grouping, and the pairwise charts with the long format.
pub fn plot_distance(files: &DistanceFiles, writer: &PlotWriter) -> Vec<PathBuf> {
    let style = writer.style();
    let summary = DistancePlotter::new(&files.summary, style);
    let mut written: Vec<PathBuf> = [DistancePlot::LocusMean, DistancePlot::MeanHistogram]
        .iter()
        .filter_map(|preset| writer.write(preset.name(), summary.build(preset)))
        .collect();
    let matrices = matrix_files(&files.per_locus_dir);
    if matrices.len() > MAX_MATRIX_PLOTS {
        log::info!(
            "Writing distance heatmaps for the first {} of {} loci",
            MAX_MATRIX_PLOTS,
            matrices.len()
        );
    }
    for matrix in matrices.iter().take(MAX_MATRIX_PLOTS) {
        let spec = DistancePlotter::new(matrix, style).build(&DistancePlot::Matrix);
        // A single sample has no pairwise distance to show.
        if let Ok(spec) = &spec
            && spec.heatmap.as_ref().is_some_and(|h| h.rows.len() < 2)
        {
            continue;
        }
        let locus = matrix
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or_default();
        written.extend(writer.write(&format!("matrix_{}", locus), spec));
    }
    if let Some(group_mean) = &files.group_mean {
        let spec = DistancePlotter::new(group_mean, style)
            .build(&DistancePlot::Matrix)
            .map(|mut spec| {
                spec.title = "Group mean pairwise distance".to_string();
                spec
            });
        written.extend(writer.write("group-mean", spec));
    }
    if let Some(long) = &files.pairwise_long {
        let plotter = DistancePlotter::new(long, style);
        for preset in [DistancePlot::PairwiseHistogram, DistancePlot::LocusBox] {
            written.extend(writer.write(preset.name(), plotter.build(&preset)));
        }
    }
    written
}

/// Per-locus matrix files, in locus order.
fn matrix_files(dir: &Path) -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = match fs::read_dir(dir) {
        Ok(entries) => entries
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.extension().is_some_and(|e| e == "csv"))
            .collect(),
        Err(error) => {
            log::warn!("Skipped distance heatmaps for {}: {}", dir.display(), error);
            Vec::new()
        }
    };
    files.sort_by(|a, b| alphanumeric_sort::compare_path(a, b));
    files
}

/// Builds pairwise distance charts.
pub struct DistancePlotter<'a> {
    input: &'a Path,
    style: &'a PlotStyle,
}

impl<'a> DistancePlotter<'a> {
    pub fn new(input: &'a Path, style: &'a PlotStyle) -> Self {
        Self { input, style }
    }

    pub fn build(&self, preset: &DistancePlot) -> Result<PlotSpec> {
        let csv = SummaryCsv::from_path(self.input)?;
        if csv.is_empty() {
            bail!("No rows in {}", self.input.display());
        }
        let spec = match preset {
            DistancePlot::LocusMean => bar_spec(
                "Mean pairwise distance per locus",
                ("Locus", "Mean distance"),
                csv.column_str("locus")?,
                vec![Series::categorical("Mean", csv.column_f64("mean")?)],
                self.style,
            ),
            DistancePlot::MeanHistogram => histogram_spec(
                "Mean pairwise distance per locus",
                "Mean distance",
                &csv.column_f64("mean")?,
                self.style,
            ),
            DistancePlot::Matrix => self.matrix(&csv)?,
            DistancePlot::PairwiseHistogram => histogram_spec(
                "Pairwise distances",
                "Distance",
                &csv.column_f64("distance")?,
                self.style,
            ),
            DistancePlot::LocusBox => self.locus_box(&csv)?,
        };
        Ok(spec)
    }

    /// Heatmap of a square matrix whose first column holds the row names.
    /// Empty cells (unwritten triangle or blank diagonal) are missing.
    fn matrix(&self, csv: &SummaryCsv) -> Result<PlotSpec> {
        let cols: Vec<String> = csv.headers().iter().skip(1).cloned().collect();
        if cols.is_empty() {
            bail!("{} is not a distance matrix", self.input.display());
        }
        let rows = csv.column_at(0);
        let values = (1..=cols.len())
            .map(|i| {
                csv.column_at(i)
                    .iter()
                    .map(|c| c.parse::<f64>().unwrap_or(f64::NAN))
                    .collect::<Vec<f64>>()
            })
            .collect::<Vec<Vec<f64>>>();
        // Column-major to row-major.
        let cells = (0..rows.len())
            .flat_map(|r| values.iter().map(move |col| col[r]))
            .collect();
        let title = format!("Pairwise distance: {}", self.locus_name());
        let mut spec = PlotSpec::new(&title, ChartKind::Heatmap, "", "");
        spec.heatmap = Some(HeatmapData::new(rows, cols, cells)?);
        spec.style = self.style.clone();
        Ok(spec)
    }

    fn locus_box(&self, csv: &SummaryCsv) -> Result<PlotSpec> {
        let loci = csv.column_str("locus")?;
        let distances = csv.column_f64("distance")?;
        let mut groups: IndexMap<String, Vec<f64>> = IndexMap::new();
        loci.into_iter()
            .zip(distances)
            .for_each(|(locus, d)| groups.entry(locus).or_default().push(d));
        let mut spec = PlotSpec::new(
            "Pairwise distances per locus",
            ChartKind::BoxPlot,
            "Locus",
            "Distance",
        );
        spec.boxes = groups
            .iter()
            .filter_map(|(locus, values)| BoxStats::from_values(locus, values))
            .collect();
        spec.style = self.style.clone();
        Ok(spec)
    }

    fn locus_name(&self) -> String {
        self.input
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or_default()
            .to_string()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::plot::fonts::PlotFonts;
    use crate::plot::spec::PlotFormat;

    #[test]
    fn test_plot_distance() {
        let dir = tempfile::TempDir::with_prefix("plot").unwrap();
        let per_locus = dir.path().join("per_locus");
        fs::create_dir_all(&per_locus).unwrap();
        fs::copy(
            "tests/files/plot/distance_matrix.csv",
            per_locus.join("gene_2.csv"),
        )
        .unwrap();
        fs::write(per_locus.join("gene_1.csv"), ",ABCD\nABCD,0.000000\n").unwrap();
        let plots = dir.path().join("plots");
        let style = PlotStyle::default();
        let fonts = PlotFonts::new();
        let writer = PlotWriter::new(&plots, None, &PlotFormat::Svg, &style, &fonts);
        let files = DistanceFiles {
            summary: PathBuf::from("tests/files/plot/distance_summary.csv"),
            per_locus_dir: per_locus,
            group_mean: None,
            pairwise_long: Some(PathBuf::from("tests/files/plot/distance_pairwise_long.csv")),
        };
        let written = plot_distance(&files, &writer);
        assert_eq!(written.len(), 5);
        assert!(plots.join("matrix_gene_2.svg").is_file());
        assert!(!plots.join("matrix_gene_1.svg").exists());
        assert!(plots.join("locus-box.svg").is_file());
    }

    fn build(input: &str, preset: DistancePlot) -> PlotSpec {
        let style = PlotStyle::default();
        let spec = DistancePlotter::new(Path::new(input), &style)
            .build(&preset)
            .unwrap();
        spec.validate().unwrap();
        spec
    }

    #[test]
    fn test_locus_mean() {
        let spec = build(
            "tests/files/plot/distance_summary.csv",
            DistancePlot::LocusMean,
        );
        assert_eq!(spec.categories.len(), 4);
        assert!(spec.series[0].y[0].is_nan());
        build(
            "tests/files/plot/distance_summary.csv",
            DistancePlot::MeanHistogram,
        );
    }

    #[test]
    fn test_matrix() {
        let spec = build("tests/files/plot/distance_matrix.csv", DistancePlot::Matrix);
        let heatmap = spec.heatmap.unwrap();
        assert_eq!(heatmap.rows, vec!["ABCD", "ABCF", "ABCE"]);
        assert_eq!(heatmap.cols, heatmap.rows);
        assert_eq!(heatmap.get(0, 1), Some(0.125));
        assert_eq!(heatmap.get(1, 1), Some(0.0));
    }

    #[test]
    fn test_pairwise() {
        let input = "tests/files/plot/distance_pairwise_long.csv";
        let hist = build(input, DistancePlot::PairwiseHistogram);
        assert_eq!(hist.series[0].y.iter().sum::<f64>(), 6.0);
        let boxes = build(input, DistancePlot::LocusBox);
        let labels: Vec<&str> = boxes.boxes.iter().map(|b| b.label.as_str()).collect();
        assert_eq!(labels, vec!["gene_1", "gene_2", "gene_3", "gene_4"]);
    }
}
