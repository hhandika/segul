//! Renderer-agnostic chart description.
//!
//! A [`PlotSpec`] holds everything needed to draw a chart: the chart kind,
//! axis titles, precomputed series, and the color style. The same spec is
//! rendered to SVG, PNG, or PDF by [`PlotRenderer`](super::render::PlotRenderer)
//! and can be passed as-is to other front ends (for example, Flutter through
//! flutter_rust_bridge) so they draw exactly the same values.
use std::path::Path;

use anyhow::{Result, bail};
use serde::{Deserialize, Serialize};

use super::color::{Colormap, CvdKind, Palette, Rgb, simulate_cvd};

/// Default chart width in pixels.
pub const DEFAULT_WIDTH: u32 = 960;
/// Default chart height in pixels.
pub const DEFAULT_HEIGHT: u32 = 600;
/// Default PNG scale factor (2x for sharp output on high-density screens).
pub const DEFAULT_PNG_SCALE: f32 = 2.0;

/// Supported chart kinds.
///
/// How a kind reads the data in [`PlotSpec`]:
/// - `Scatter` and `Line`: numeric `x` and `y` in each series. `Line` draws
///   a shaded band when `lower` and `upper` are set.
/// - `Bar` and `StackedBar`: `y` values aligned to `categories`.
///   Several series are drawn side by side (`Bar`) or stacked (`StackedBar`).
/// - `Histogram`: one series with `x` holding the bin edges
///   (one more than the counts) and `y` holding the counts.
/// - `Heatmap`: values in `heatmap`.
/// - `BoxPlot`: precomputed boxes in `boxes`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ChartKind {
    Scatter,
    Line,
    Bar,
    StackedBar,
    Histogram,
    Heatmap,
    BoxPlot,
}

/// Output file formats for rendered plots.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PlotFormat {
    Svg,
    Png,
    Pdf,
}

impl PlotFormat {
    /// Infer the plot format from a file extension (svg, png, or pdf).
    pub fn from_path(path: &Path) -> Result<Self> {
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or_default()
            .to_lowercase();
        ext.parse()
    }

    /// File extension for the format, without the leading dot.
    pub fn extension(&self) -> &'static str {
        match self {
            Self::Svg => "svg",
            Self::Png => "png",
            Self::Pdf => "pdf",
        }
    }
}

impl std::fmt::Display for PlotFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.extension())
    }
}

impl std::str::FromStr for PlotFormat {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "svg" => Ok(Self::Svg),
            "png" => Ok(Self::Png),
            "pdf" => Ok(Self::Pdf),
            _ => bail!("Unsupported plot format '{}'. Use svg, png, or pdf", s),
        }
    }
}

/// Axis title.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Axis {
    pub label: String,
}

impl Axis {
    pub fn new(label: &str) -> Self {
        Self {
            label: label.to_string(),
        }
    }
}

/// One data series.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Series {
    /// Series name, shown in the legend when there is more than one series.
    pub name: String,
    /// Numeric x values. Empty for categorical charts.
    #[serde(with = "nan_as_null")]
    pub x: Vec<f64>,
    /// Numeric y values. NaN marks a missing value.
    #[serde(with = "nan_as_null")]
    pub y: Vec<f64>,
    /// Lower band bound for line charts.
    #[serde(with = "optional_nan_as_null")]
    pub lower: Option<Vec<f64>>,
    /// Upper band bound for line charts.
    #[serde(with = "optional_nan_as_null")]
    pub upper: Option<Vec<f64>>,
}

impl Series {
    /// Series with numeric x and y values.
    pub fn xy(name: &str, x: Vec<f64>, y: Vec<f64>) -> Self {
        Self {
            name: name.to_string(),
            x,
            y,
            ..Default::default()
        }
    }

    /// Series of y values aligned to the chart categories.
    pub fn categorical(name: &str, y: Vec<f64>) -> Self {
        Self {
            name: name.to_string(),
            y,
            ..Default::default()
        }
    }

    /// Histogram series with equal-width bins.
    ///
    /// Non-finite values are skipped. When `bins` is zero, the bin count
    /// follows Sturges' rule.
    pub fn histogram(name: &str, values: &[f64], bins: usize) -> Self {
        let values: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
        if values.is_empty() {
            return Self::xy(name, Vec::new(), Vec::new());
        }
        let bins = if bins == 0 {
            ((values.len() as f64).log2().ceil() as usize + 1).max(1)
        } else {
            bins
        };
        let (min, max) = finite_range(&values).expect("Failed computing histogram range");
        // A single unique value still gets a visible bin.
        let (min, max) = if min == max {
            (min - 0.5, max + 0.5)
        } else {
            (min, max)
        };
        let width = (max - min) / bins as f64;
        let edges: Vec<f64> = (0..=bins).map(|i| min + width * i as f64).collect();
        let mut counts = vec![0.0; bins];
        values.iter().for_each(|v| {
            let idx = (((v - min) / width) as usize).min(bins - 1);
            counts[idx] += 1.0;
        });
        Self::xy(name, edges, counts)
    }
}

/// Values for a heatmap, stored row-major. `None` marks a missing cell.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct HeatmapData {
    pub rows: Vec<String>,
    pub cols: Vec<String>,
    pub values: Vec<Option<f64>>,
}

impl HeatmapData {
    /// Build a heatmap from row-major values. Non-finite values become missing cells.
    pub fn new(rows: Vec<String>, cols: Vec<String>, values: Vec<f64>) -> Result<Self> {
        if rows.len() * cols.len() != values.len() {
            bail!(
                "Heatmap has {} values but {} rows x {} columns",
                values.len(),
                rows.len(),
                cols.len()
            );
        }
        let values = values
            .into_iter()
            .map(|v| if v.is_finite() { Some(v) } else { None })
            .collect();
        Ok(Self { rows, cols, values })
    }

    pub fn get(&self, row: usize, col: usize) -> Option<f64> {
        self.values[row * self.cols.len() + col]
    }

    /// Minimum and maximum of the values.
    pub fn range(&self) -> Option<(f64, f64)> {
        let values: Vec<f64> = self.values.iter().flatten().copied().collect();
        finite_range(&values)
    }
}

/// Precomputed box plot statistics using Tukey's 1.5 IQR whiskers.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BoxStats {
    pub label: String,
    pub lower_whisker: f64,
    pub q1: f64,
    pub median: f64,
    pub q3: f64,
    pub upper_whisker: f64,
    pub outliers: Vec<f64>,
}

impl BoxStats {
    /// Compute box statistics from raw values. Non-finite values are skipped.
    /// Returns `None` when no finite value is left.
    pub fn from_values(label: &str, values: &[f64]) -> Option<Self> {
        let mut values: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
        if values.is_empty() {
            return None;
        }
        values.sort_by(|a, b| a.total_cmp(b));
        let q1 = quantile(&values, 0.25);
        let median = quantile(&values, 0.5);
        let q3 = quantile(&values, 0.75);
        let iqr = q3 - q1;
        let (low_fence, high_fence) = (q1 - 1.5 * iqr, q3 + 1.5 * iqr);
        let inside = values
            .iter()
            .filter(|v| **v >= low_fence && **v <= high_fence);
        let lower_whisker = inside.clone().copied().fold(f64::INFINITY, f64::min);
        let upper_whisker = inside.copied().fold(f64::NEG_INFINITY, f64::max);
        let outliers = values
            .iter()
            .copied()
            .filter(|v| *v < low_fence || *v > high_fence)
            .collect();
        Some(Self {
            label: label.to_string(),
            lower_whisker,
            q1,
            median,
            q3,
            upper_whisker,
            outliers,
        })
    }
}

/// Colors and size used to draw a plot.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PlotStyle {
    /// Categorical palette for series.
    pub palette: Palette,
    /// Continuous colormap for heatmaps.
    pub colormap: Colormap,
    /// Simulate a color vision deficiency on every color.
    pub cvd: Option<CvdKind>,
    pub width: u32,
    pub height: u32,
    /// Scale factor applied when rasterizing to PNG.
    pub png_scale: f32,
}

impl Default for PlotStyle {
    fn default() -> Self {
        Self {
            palette: Palette::OkabeIto,
            colormap: Colormap::Viridis,
            cvd: None,
            width: DEFAULT_WIDTH,
            height: DEFAULT_HEIGHT,
            png_scale: DEFAULT_PNG_SCALE,
        }
    }
}

impl PlotStyle {
    /// Series colors after applying the color vision deficiency simulation.
    pub fn series_colors(&self, n: usize) -> Vec<Rgb> {
        self.palette
            .colors(n)
            .into_iter()
            .map(|c| self.apply_cvd(c))
            .collect()
    }

    /// Colormap color at `t` in [0, 1] after the color vision deficiency simulation.
    pub fn colormap_at(&self, t: f64) -> Rgb {
        self.apply_cvd(self.colormap.at(t))
    }

    /// Evenly spaced colormap stops, for drawing gradients in other front ends.
    pub fn colormap_stops(&self, n: usize) -> Vec<Rgb> {
        let n = n.max(2);
        (0..n)
            .map(|i| self.colormap_at(i as f64 / (n - 1) as f64))
            .collect()
    }

    /// Apply the color vision deficiency simulation, if any.
    pub fn apply_cvd(&self, color: Rgb) -> Rgb {
        match &self.cvd {
            Some(kind) => simulate_cvd(color, kind, 1.0),
            None => color,
        }
    }
}

/// Complete description of a chart.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PlotSpec {
    pub title: String,
    pub kind: ChartKind,
    pub x_axis: Axis,
    pub y_axis: Axis,
    /// Category labels for bar charts and box plots.
    pub categories: Vec<String>,
    pub series: Vec<Series>,
    pub heatmap: Option<HeatmapData>,
    pub boxes: Vec<BoxStats>,
    pub style: PlotStyle,
}

impl PlotSpec {
    /// Empty spec of the given kind. Fill in the data fields afterward.
    pub fn new(title: &str, kind: ChartKind, x_label: &str, y_label: &str) -> Self {
        Self {
            title: title.to_string(),
            kind,
            x_axis: Axis::new(x_label),
            y_axis: Axis::new(y_label),
            categories: Vec::new(),
            series: Vec::new(),
            heatmap: None,
            boxes: Vec::new(),
            style: PlotStyle::default(),
        }
    }

    /// Serialize the spec to JSON, for front ends that do not link the crate.
    pub fn to_json(&self) -> Result<String> {
        Ok(serde_json::to_string(self)?)
    }

    /// Parse a spec from JSON.
    pub fn from_json(json: &str) -> Result<Self> {
        Ok(serde_json::from_str(json)?)
    }

    /// Check that the data fields match the chart kind.
    pub fn validate(&self) -> Result<()> {
        match self.kind {
            ChartKind::Scatter | ChartKind::Line => {
                if self.series.is_empty() {
                    bail!("No series to plot");
                }
                for s in &self.series {
                    if s.x.len() != s.y.len() {
                        bail!("Series '{}' has different x and y lengths", s.name);
                    }
                }
            }
            ChartKind::Bar | ChartKind::StackedBar => {
                if self.categories.is_empty() || self.series.is_empty() {
                    bail!("No categories to plot");
                }
                for s in &self.series {
                    if s.y.len() != self.categories.len() {
                        bail!(
                            "Series '{}' has {} values for {} categories",
                            s.name,
                            s.y.len(),
                            self.categories.len()
                        );
                    }
                }
            }
            ChartKind::Histogram => match self.series.first() {
                Some(s) if !s.y.is_empty() && s.x.len() == s.y.len() + 1 => (),
                _ => bail!("Histogram needs one series with bin edges and counts"),
            },
            ChartKind::Heatmap => match &self.heatmap {
                Some(h) if !h.rows.is_empty() && !h.cols.is_empty() => (),
                _ => bail!("No heatmap values to plot"),
            },
            ChartKind::BoxPlot => {
                if self.boxes.is_empty() {
                    bail!("No box plot values to plot");
                }
            }
        }
        Ok(())
    }
}

// JSON has no NaN, so missing values travel as null.
mod nan_as_null {
    use serde::{Deserialize, Deserializer, Serializer};

    pub fn serialize<S: Serializer>(values: &[f64], serializer: S) -> Result<S::Ok, S::Error> {
        serializer.collect_seq(values.iter().map(|v| v.is_finite().then_some(*v)))
    }

    pub fn deserialize<'de, D: Deserializer<'de>>(deserializer: D) -> Result<Vec<f64>, D::Error> {
        let values: Vec<Option<f64>> = Vec::deserialize(deserializer)?;
        Ok(values.into_iter().map(|v| v.unwrap_or(f64::NAN)).collect())
    }
}

mod optional_nan_as_null {
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    pub fn serialize<S: Serializer>(
        values: &Option<Vec<f64>>,
        serializer: S,
    ) -> Result<S::Ok, S::Error> {
        values
            .as_ref()
            .map(|v| {
                v.iter()
                    .map(|x| x.is_finite().then_some(*x))
                    .collect::<Vec<_>>()
            })
            .serialize(serializer)
    }

    pub fn deserialize<'de, D: Deserializer<'de>>(
        deserializer: D,
    ) -> Result<Option<Vec<f64>>, D::Error> {
        let values: Option<Vec<Option<f64>>> = Option::deserialize(deserializer)?;
        Ok(values.map(|v| v.into_iter().map(|x| x.unwrap_or(f64::NAN)).collect()))
    }
}

/// Minimum and maximum of the finite values.
pub(crate) fn finite_range(values: &[f64]) -> Option<(f64, f64)> {
    values
        .iter()
        .filter(|v| v.is_finite())
        .fold(None, |acc, &v| match acc {
            None => Some((v, v)),
            Some((min, max)) => Some((min.min(v), max.max(v))),
        })
}

/// Linear interpolation quantile (type 7) of sorted values.
fn quantile(sorted: &[f64], q: f64) -> f64 {
    let pos = (sorted.len() - 1) as f64 * q;
    let lower = pos.floor() as usize;
    let upper = pos.ceil() as usize;
    sorted[lower] + (sorted[upper] - sorted[lower]) * (pos - lower as f64)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_histogram_bins() {
        let values = [1.0, 2.0, 2.0, 3.0, 4.0, f64::NAN];
        let hist = Series::histogram("x", &values, 3);
        assert_eq!(hist.x, vec![1.0, 2.0, 3.0, 4.0]);
        assert_eq!(hist.y, vec![1.0, 2.0, 2.0]);
    }

    #[test]
    fn test_histogram_single_value() {
        let hist = Series::histogram("x", &[5.0, 5.0], 0);
        assert_eq!(hist.y.iter().sum::<f64>(), 2.0);
        assert_eq!(hist.x.len(), hist.y.len() + 1);
    }

    #[test]
    fn test_box_stats() {
        let values = [1.0, 2.0, 3.0, 4.0, 5.0, 100.0];
        let stats = BoxStats::from_values("x", &values).unwrap();
        assert_eq!(stats.q1, 2.25);
        assert_eq!(stats.median, 3.5);
        assert_eq!(stats.q3, 4.75);
        assert_eq!(stats.lower_whisker, 1.0);
        assert_eq!(stats.upper_whisker, 5.0);
        assert_eq!(stats.outliers, vec![100.0]);
    }

    #[test]
    fn test_box_stats_empty() {
        assert!(BoxStats::from_values("x", &[f64::NAN]).is_none());
    }

    #[test]
    fn test_plot_format() {
        assert_eq!(
            PlotFormat::from_path(Path::new("a/plot.PDF")).unwrap(),
            PlotFormat::Pdf
        );
        assert!(PlotFormat::from_path(Path::new("plot.jpg")).is_err());
    }

    #[test]
    fn test_spec_json_round_trip() {
        let mut spec = PlotSpec::new("t", ChartKind::Line, "x", "y");
        let mut series = Series::xy("s", vec![1.0, 2.0], vec![2.0, f64::NAN]);
        series.lower = Some(vec![1.0, f64::NAN]);
        spec.series.push(series);
        let json = spec.to_json().unwrap();
        let parsed = PlotSpec::from_json(&json).unwrap();
        assert_eq!(parsed.series[0].x, vec![1.0, 2.0]);
        assert!(parsed.series[0].y[1].is_nan());
        assert!(parsed.series[0].lower.as_ref().unwrap()[1].is_nan());
        assert_eq!(parsed.series[0].upper, None);
    }

    #[test]
    fn test_heatmap_missing_round_trip() {
        let mut spec = PlotSpec::new("t", ChartKind::Heatmap, "x", "y");
        let heatmap = HeatmapData::new(
            vec!["a".into()],
            vec!["b".into(), "c".into()],
            vec![1.0, f64::NAN],
        )
        .unwrap();
        assert_eq!(heatmap.get(0, 1), None);
        spec.heatmap = Some(heatmap);
        let json = spec.to_json().unwrap();
        assert_eq!(PlotSpec::from_json(&json).unwrap(), spec);
    }

    #[test]
    fn test_validate_bar_length() {
        let mut spec = PlotSpec::new("t", ChartKind::Bar, "x", "y");
        spec.categories = vec!["a".into(), "b".into()];
        spec.series.push(Series::categorical("s", vec![1.0]));
        assert!(spec.validate().is_err());
    }
}
