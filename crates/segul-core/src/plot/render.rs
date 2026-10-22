//! Render a [`PlotSpec`] to SVG, PNG, or PDF.
//!
//! Charts are drawn once as SVG with `plotters`. PNG (`resvg`) and PDF
//! (`svg2pdf`, behind the `pdf` feature) are converted from that SVG, using
//! the fonts given in [`PlotFonts`].
use std::fs;
use std::ops::Range;
use std::path::Path;

use anyhow::{Context, Result, anyhow};
use plotters::coord::Shift;
use plotters::coord::ranged1d::{DefaultFormatting, KeyPointHint};
use plotters::coord::types::RangedCoordf64;
use plotters::prelude::*;
use plotters::style::text_anchor::{HPos, Pos, VPos};
use resvg::tiny_skia;

use super::color::Rgb;
use super::fonts::PlotFonts;
use super::spec::{ChartKind, PlotFormat, PlotSpec, finite_range};

// Generic family in the SVG, so standalone SVG files use the viewer's
// sans-serif font. PlotFonts maps it to a real font for PNG and PDF.
const FONT: &str = "sans-serif";
const TITLE_SIZE: f64 = 22.0;
const AXIS_TITLE_SIZE: f64 = 17.0;
const LABEL_SIZE: f64 = 14.0;
const MARGIN: u32 = 16;
const Y_LABEL_AREA: u32 = 80;
const X_LABEL_AREA: u32 = 56;
// Approximate legend row height and padding in pixels.
const LEGEND_ROW: u32 = 20;
const LEGEND_PADDING: u32 = 24;
const TITLE_AREA: u32 = 40;
// Approximate label character width in pixels at LABEL_SIZE.
const CHAR_WIDTH: u32 = 8;
const MAX_LABEL_CHARS: usize = 24;
const COLORBAR_SPACE: u32 = 110;
const COLORBAR_WIDTH: i32 = 20;
// Heatmaps with more cells than this skip value annotations.
const MAX_ANNOTATED_CELLS: usize = 144;
const MISSING_CELL: RGBColor = RGBColor(220, 220, 220);
const GRID: RGBColor = RGBColor(225, 225, 225);
const AXIS: RGBColor = RGBColor(90, 90, 90);

type Area<'b> = DrawingArea<SVGBackend<'b>, Shift>;

/// Categorical axis: category `i` sits at `i`, with a tick on every category.
#[derive(Clone)]
struct CategoryCoord {
    count: usize,
}

impl Ranged for CategoryCoord {
    type FormatOption = DefaultFormatting;
    type ValueType = f64;

    fn map(&self, value: &f64, limit: (i32, i32)) -> i32 {
        let range = self.range();
        let t = (value - range.start) / (range.end - range.start);
        limit.0 + ((limit.1 - limit.0) as f64 * t).round() as i32
    }

    fn key_points<Hint: KeyPointHint>(&self, _hint: Hint) -> Vec<f64> {
        (0..self.count).map(|i| i as f64).collect()
    }

    fn range(&self) -> Range<f64> {
        -0.5..self.count as f64 - 0.5
    }
}

/// Draws a [`PlotSpec`] and writes it as SVG, PNG, or PDF.
///
/// # Example
/// ```rust
/// use segul_core::plot::fonts::PlotFonts;
/// use segul_core::plot::render::PlotRenderer;
/// use segul_core::plot::spec::{ChartKind, PlotFormat, PlotSpec, Series};
///
/// let mut spec = PlotSpec::new("Example", ChartKind::Scatter, "x", "y");
/// spec.series.push(Series::xy("points", vec![1.0, 2.0, 3.0], vec![2.0, 4.0, 3.0]));
/// // SVG needs no fonts. For PNG and PDF, pass a font from the app
/// // with `PlotFonts::from_data` or use `PlotFonts::system()`.
/// let fonts = PlotFonts::new();
/// let svg = PlotRenderer::new(&spec, &fonts)
///     .render_bytes(&PlotFormat::Svg)
///     .unwrap();
/// assert!(String::from_utf8(svg).unwrap().contains("<svg"));
/// ```
pub struct PlotRenderer<'a> {
    spec: &'a PlotSpec,
    fonts: &'a PlotFonts,
}

impl<'a> PlotRenderer<'a> {
    pub fn new(spec: &'a PlotSpec, fonts: &'a PlotFonts) -> Self {
        Self { spec, fonts }
    }

    /// Render the plot and write it to `output`, creating parent directories.
    pub fn render(&self, output: &Path, format: &PlotFormat) -> Result<()> {
        let bytes = self.render_bytes(format)?;
        if let Some(parent) = output.parent() {
            fs::create_dir_all(parent).with_context(|| {
                format!(
                    "Failed creating an output directory for {}",
                    output.display()
                )
            })?;
        }
        fs::write(output, bytes)
            .with_context(|| format!("Failed writing plot to {}", output.display()))
    }

    /// Render the plot to bytes in the given format.
    /// PDF needs the `pdf` feature.
    pub fn render_bytes(&self, format: &PlotFormat) -> Result<Vec<u8>> {
        let svg = self.to_svg()?;
        match format {
            PlotFormat::Svg => Ok(svg.into_bytes()),
            PlotFormat::Png => to_png(&self.parse_svg(&svg)?, self.spec.style.png_scale),
            PlotFormat::Pdf => to_pdf(&self.parse_svg(&svg)?),
        }
    }

    fn parse_svg(&self, svg: &str) -> Result<usvg::Tree> {
        let options = self.fonts.svg_options()?;
        usvg::Tree::from_str(svg, &options).context("Failed parsing the rendered SVG")
    }

    /// Render the plot as an SVG string.
    pub fn to_svg(&self) -> Result<String> {
        self.spec.validate()?;
        let size = (self.spec.style.width, self.spec.style.height);
        let mut svg = String::new();
        {
            let root = SVGBackend::with_string(&mut svg, size).into_drawing_area();
            root.fill(&WHITE).map_err(draw_error)?;
            match self.spec.kind {
                ChartKind::Scatter | ChartKind::Line => self.draw_xy(&root)?,
                ChartKind::Bar | ChartKind::StackedBar => self.draw_bars(&root)?,
                ChartKind::Histogram => self.draw_histogram(&root)?,
                ChartKind::Heatmap => self.draw_heatmap(&root)?,
                ChartKind::BoxPlot => self.draw_boxes(&root)?,
            }
            root.present().map_err(draw_error)?;
        }
        Ok(svg)
    }

    fn draw_xy(&self, root: &Area) -> Result<()> {
        let spec = self.spec;
        let xs: Vec<f64> = spec.series.iter().flat_map(|s| s.x.clone()).collect();
        let ys: Vec<f64> = spec
            .series
            .iter()
            .flat_map(|s| {
                let band = s.lower.iter().chain(s.upper.iter()).flatten();
                s.y.iter().chain(band).copied().collect::<Vec<f64>>()
            })
            .collect();
        let x_range = padded_range(&xs, false);
        let y_range = self.legend_headroom(padded_range(&ys, false));
        let mut chart = ChartBuilder::on(root)
            .caption(&spec.title, (FONT, TITLE_SIZE))
            .margin(MARGIN)
            .x_label_area_size(X_LABEL_AREA)
            .y_label_area_size(Y_LABEL_AREA)
            .build_cartesian_2d(x_range, y_range)
            .map_err(draw_error)?;
        chart
            .configure_mesh()
            .x_desc(&spec.x_axis.label)
            .y_desc(&spec.y_axis.label)
            .x_label_formatter(&format_tick)
            .y_label_formatter(&format_tick)
            .label_style((FONT, LABEL_SIZE))
            .axis_desc_style((FONT, AXIS_TITLE_SIZE))
            .axis_style(AXIS)
            .bold_line_style(GRID)
            .light_line_style(TRANSPARENT)
            .draw()
            .map_err(draw_error)?;

        let colors = spec.style.series_colors(spec.series.len());
        for (series, color) in spec.series.iter().zip(colors) {
            let color = to_rgb(color);
            let points: Vec<(f64, f64)> = series
                .x
                .iter()
                .zip(&series.y)
                .filter(|(x, y)| x.is_finite() && y.is_finite())
                .map(|(x, y)| (*x, *y))
                .collect();
            let anno = match spec.kind {
                ChartKind::Line => {
                    if let (Some(lower), Some(upper)) = (&series.lower, &series.upper) {
                        let band: Vec<(f64, f64)> = series
                            .x
                            .iter()
                            .zip(upper)
                            .chain(series.x.iter().zip(lower).rev())
                            .filter(|(x, y)| x.is_finite() && y.is_finite())
                            .map(|(x, y)| (*x, *y))
                            .collect();
                        chart
                            .draw_series(std::iter::once(Polygon::new(
                                band,
                                color.mix(0.2).filled(),
                            )))
                            .map_err(draw_error)?;
                    }
                    chart
                        .draw_series(LineSeries::new(points, color.stroke_width(2)))
                        .map_err(draw_error)?
                }
                _ => chart
                    .draw_series(
                        points
                            .into_iter()
                            .map(|p| Circle::new(p, 4, color.mix(0.8).filled())),
                    )
                    .map_err(draw_error)?,
            };
            anno.label(&series.name).legend(move |(x, y)| {
                Rectangle::new([(x, y - 5), (x + 12, y + 5)], color.filled())
            });
        }
        if spec.series.len() > 1 {
            chart
                .configure_series_labels()
                .position(SeriesLabelPosition::UpperRight)
                .background_style(WHITE.mix(0.85))
                .border_style(AXIS.mix(0.4))
                .label_font((FONT, LABEL_SIZE))
                .draw()
                .map_err(draw_error)?;
        }
        Ok(())
    }

    fn draw_bars(&self, root: &Area) -> Result<()> {
        let spec = self.spec;
        let stacked = spec.kind == ChartKind::StackedBar;
        let n_series = spec.series.len();
        let ys: Vec<f64> = if stacked {
            (0..spec.categories.len())
                .map(|i| {
                    spec.series
                        .iter()
                        .map(|s| s.y[i])
                        .filter(|y| y.is_finite())
                        .sum()
                })
                .collect()
        } else {
            spec.series.iter().flat_map(|s| s.y.clone()).collect()
        };
        let y_range = self.legend_headroom(padded_range(&ys, true));
        let y_start = y_range.start;
        let mut chart = self.categorical_chart(root, &spec.categories, y_range)?;

        let colors = spec.style.series_colors(n_series);
        let group_width = 0.8;
        let bar_width = if stacked {
            group_width
        } else {
            group_width / n_series as f64
        };
        let mut bases = vec![0.0; spec.categories.len()];
        for (j, (series, color)) in spec.series.iter().zip(colors).enumerate() {
            let color = to_rgb(color);
            let bars: Vec<[(f64, f64); 2]> = series
                .y
                .iter()
                .enumerate()
                .filter(|(_, y)| y.is_finite())
                .map(|(i, y)| {
                    let center = i as f64;
                    if stacked {
                        let base = bases[i];
                        bases[i] += y;
                        [
                            (center - bar_width / 2.0, base),
                            (center + bar_width / 2.0, base + y),
                        ]
                    } else {
                        let x0 = center - group_width / 2.0 + bar_width * j as f64;
                        [(x0, 0.0), (x0 + bar_width, *y)]
                    }
                })
                .collect();
            chart
                .draw_series(
                    bars.into_iter()
                        .map(|corners| Rectangle::new(corners, color.filled())),
                )
                .map_err(draw_error)?
                .label(&series.name)
                .legend(move |(x, y)| {
                    Rectangle::new([(x, y - 5), (x + 12, y + 5)], color.filled())
                });
        }
        if n_series > 1 {
            chart
                .configure_series_labels()
                .position(SeriesLabelPosition::UpperRight)
                .background_style(WHITE.mix(0.85))
                .border_style(AXIS.mix(0.4))
                .label_font((FONT, LABEL_SIZE))
                .draw()
                .map_err(draw_error)?;
        }
        self.draw_x_categories(root, &chart, &spec.categories, y_start)
    }

    fn draw_histogram(&self, root: &Area) -> Result<()> {
        let spec = self.spec;
        let series = &spec.series[0];
        let x_range = padded_range(&series.x, false);
        let y_range = padded_range(&series.y, true);
        let mut chart = ChartBuilder::on(root)
            .caption(&spec.title, (FONT, TITLE_SIZE))
            .margin(MARGIN)
            .x_label_area_size(X_LABEL_AREA)
            .y_label_area_size(Y_LABEL_AREA)
            .build_cartesian_2d(x_range, y_range)
            .map_err(draw_error)?;
        chart
            .configure_mesh()
            .disable_x_mesh()
            .x_desc(&spec.x_axis.label)
            .y_desc(&spec.y_axis.label)
            .x_label_formatter(&format_tick)
            .y_label_formatter(&format_tick)
            .label_style((FONT, LABEL_SIZE))
            .axis_desc_style((FONT, AXIS_TITLE_SIZE))
            .axis_style(AXIS)
            .bold_line_style(GRID)
            .light_line_style(TRANSPARENT)
            .draw()
            .map_err(draw_error)?;
        let color = to_rgb(spec.style.series_colors(1)[0]);
        let bins: Vec<[(f64, f64); 2]> = series
            .x
            .windows(2)
            .zip(&series.y)
            .map(|(edge, count)| [(edge[0], 0.0), (edge[1], *count)])
            .collect();
        chart
            .draw_series(
                bins.iter()
                    .map(|corners| Rectangle::new(*corners, color.filled())),
            )
            .map_err(draw_error)?;
        chart
            .draw_series(
                bins.iter()
                    .map(|corners| Rectangle::new(*corners, WHITE.stroke_width(1))),
            )
            .map_err(draw_error)?;
        Ok(())
    }

    fn draw_boxes(&self, root: &Area) -> Result<()> {
        let spec = self.spec;
        let categories: Vec<String> = spec.boxes.iter().map(|b| b.label.clone()).collect();
        let values: Vec<f64> = spec
            .boxes
            .iter()
            .flat_map(|b| {
                let mut v = vec![b.lower_whisker, b.upper_whisker];
                v.extend(&b.outliers);
                v
            })
            .collect();
        let y_range = padded_range(&values, false);
        let y_start = y_range.start;
        let mut chart = self.categorical_chart(root, &categories, y_range)?;
        let colors = spec.style.series_colors(spec.boxes.len());
        let half = 0.3;
        let cap = 0.15;
        for (i, (b, color)) in spec.boxes.iter().zip(colors).enumerate() {
            let color = to_rgb(color);
            let x = i as f64;
            let lines = vec![
                vec![(x, b.lower_whisker), (x, b.q1)],
                vec![(x, b.q3), (x, b.upper_whisker)],
                vec![(x - cap, b.lower_whisker), (x + cap, b.lower_whisker)],
                vec![(x - cap, b.upper_whisker), (x + cap, b.upper_whisker)],
            ];
            chart
                .draw_series(
                    lines
                        .into_iter()
                        .map(|l| PathElement::new(l, AXIS.stroke_width(1))),
                )
                .map_err(draw_error)?;
            let body = [(x - half, b.q1), (x + half, b.q3)];
            chart
                .draw_series([
                    Rectangle::new(body, color.mix(0.75).filled()),
                    Rectangle::new(body, AXIS.stroke_width(1)),
                ])
                .map_err(draw_error)?;
            chart
                .draw_series(std::iter::once(PathElement::new(
                    vec![(x - half, b.median), (x + half, b.median)],
                    BLACK.stroke_width(2),
                )))
                .map_err(draw_error)?;
            chart
                .draw_series(
                    b.outliers
                        .iter()
                        .map(|y| Circle::new((x, *y), 3, AXIS.stroke_width(1))),
                )
                .map_err(draw_error)?;
        }
        self.draw_x_categories(root, &chart, &categories, y_start)
    }

    fn draw_heatmap(&self, root: &Area) -> Result<()> {
        let spec = self.spec;
        let heatmap = spec
            .heatmap
            .as_ref()
            .expect("Failed reading heatmap values");
        let (n_rows, n_cols) = (heatmap.rows.len(), heatmap.cols.len());
        let (min, max) = heatmap.range().unwrap_or((0.0, 1.0));
        let span = if max > min { max - min } else { 1.0 };
        let x_label_area = self.x_label_area(&heatmap.cols);
        let y_label_area = label_width(&heatmap.rows) + AXIS_TITLE_SIZE as u32 * 2;
        let mut chart = ChartBuilder::on(root)
            .caption(&spec.title, (FONT, TITLE_SIZE))
            .margin(MARGIN)
            .margin_right(COLORBAR_SPACE)
            .x_label_area_size(x_label_area)
            .y_label_area_size(y_label_area)
            .build_cartesian_2d(category_coord(n_cols), category_coord(n_rows))
            .map_err(draw_error)?;
        chart
            .configure_mesh()
            .disable_mesh()
            .x_desc(&spec.x_axis.label)
            .y_desc(&spec.y_axis.label)
            .x_label_formatter(&|_| String::new())
            .y_label_formatter(&|_| String::new())
            .axis_desc_style((FONT, AXIS_TITLE_SIZE))
            .axis_style(AXIS)
            .draw()
            .map_err(draw_error)?;

        let annotate = n_rows * n_cols <= MAX_ANNOTATED_CELLS;
        for row in 0..n_rows {
            // First row at the top.
            let y = (n_rows - 1 - row) as f64;
            for col in 0..n_cols {
                let x = col as f64;
                let value = heatmap.get(row, col);
                let fill = value.map(|v| spec.style.colormap_at((v - min) / span));
                let cell = [(x - 0.5, y - 0.5), (x + 0.5, y + 0.5)];
                let fill_color = fill.map(to_rgb).unwrap_or(MISSING_CELL);
                chart
                    .draw_series([
                        Rectangle::new(cell, fill_color.filled()),
                        Rectangle::new(cell, WHITE.stroke_width(1)),
                    ])
                    .map_err(draw_error)?;
                if let (true, Some(v), Some(fill)) = (annotate, value, fill) {
                    let text_color = if fill.luminance() > 0.4 { BLACK } else { WHITE };
                    let style = TextStyle::from((FONT, LABEL_SIZE).into_font())
                        .color(&text_color)
                        .pos(Pos::new(HPos::Center, VPos::Center));
                    chart
                        .draw_series(std::iter::once(Text::new(format_tick(&v), (x, y), style)))
                        .map_err(draw_error)?;
                }
            }
        }

        let (x_start, y_start) = (-0.5, -0.5);
        self.draw_x_categories(root, &chart, &heatmap.cols, y_start)?;
        let style = TextStyle::from((FONT, LABEL_SIZE).into_font())
            .color(&BLACK)
            .pos(Pos::new(HPos::Right, VPos::Center));
        for (row, label) in heatmap.rows.iter().enumerate() {
            let y = (n_rows - 1 - row) as f64;
            let (px, py) = chart.backend_coord(&(x_start, y));
            root.draw(&Text::new(
                truncate_label(label),
                (px - 6, py),
                style.clone(),
            ))
            .map_err(draw_error)?;
        }
        self.draw_colorbar(root, &chart, min, max)
    }

    fn draw_colorbar(
        &self,
        root: &Area,
        chart: &ChartContext<SVGBackend, Cartesian2d<CategoryCoord, CategoryCoord>>,
        min: f64,
        max: f64,
    ) -> Result<()> {
        let (x_pixels, y_pixels) = chart.plotting_area().get_pixel_range();
        let x0 = x_pixels.end + 24;
        let (top, bottom) = (y_pixels.start, y_pixels.end);
        let steps = 64;
        let step = (bottom - top) as f64 / steps as f64;
        for k in 0..steps {
            // Highest value at the top.
            let t = 1.0 - (k as f64 + 0.5) / steps as f64;
            let color = to_rgb(self.spec.style.colormap_at(t));
            let y0 = top + (step * k as f64).floor() as i32;
            let y1 = top + (step * (k + 1) as f64).ceil() as i32;
            root.draw(&Rectangle::new(
                [(x0, y0), (x0 + COLORBAR_WIDTH, y1)],
                color.filled(),
            ))
            .map_err(draw_error)?;
        }
        root.draw(&Rectangle::new(
            [(x0, top), (x0 + COLORBAR_WIDTH, bottom)],
            AXIS.stroke_width(1),
        ))
        .map_err(draw_error)?;
        let style = TextStyle::from((FONT, LABEL_SIZE).into_font())
            .color(&BLACK)
            .pos(Pos::new(HPos::Left, VPos::Center));
        let ticks = if max > min { 5 } else { 1 };
        for i in 0..ticks {
            let t = if ticks > 1 {
                i as f64 / (ticks - 1) as f64
            } else {
                0.5
            };
            let value = min + (max - min) * t;
            let y = bottom - ((bottom - top) as f64 * t).round() as i32;
            let tick_end = x0 + COLORBAR_WIDTH + 4;
            root.draw(&PathElement::new(
                vec![(x0 + COLORBAR_WIDTH, y), (tick_end, y)],
                AXIS.stroke_width(1),
            ))
            .map_err(draw_error)?;
            root.draw(&Text::new(
                format_tick(&value),
                (tick_end + 4, y),
                style.clone(),
            ))
            .map_err(draw_error)?;
        }
        Ok(())
    }

    /// Chart with a categorical x axis (one slot per category) and numeric y axis.
    fn categorical_chart<'b>(
        &self,
        root: &Area<'b>,
        categories: &[String],
        y_range: Range<f64>,
    ) -> Result<ChartContext<'b, SVGBackend<'b>, Cartesian2d<CategoryCoord, RangedCoordf64>>> {
        let spec = self.spec;
        let mut chart = ChartBuilder::on(root)
            .caption(&spec.title, (FONT, TITLE_SIZE))
            .margin(MARGIN)
            .x_label_area_size(self.x_label_area(categories))
            .y_label_area_size(Y_LABEL_AREA)
            .build_cartesian_2d(category_coord(categories.len()), y_range)
            .map_err(draw_error)?;
        chart
            .configure_mesh()
            .disable_x_mesh()
            .x_desc(&spec.x_axis.label)
            .y_desc(&spec.y_axis.label)
            // Category labels are drawn separately so long labels can rotate.
            .x_label_formatter(&|_| String::new())
            .y_label_formatter(&format_tick)
            .label_style((FONT, LABEL_SIZE))
            .axis_desc_style((FONT, AXIS_TITLE_SIZE))
            .axis_style(AXIS)
            .bold_line_style(GRID)
            .light_line_style(TRANSPARENT)
            .draw()
            .map_err(draw_error)?;
        Ok(chart)
    }

    /// Draw category labels under the x axis, rotated when they would overlap.
    fn draw_x_categories<Y: Ranged<ValueType = f64>>(
        &self,
        root: &Area,
        chart: &ChartContext<SVGBackend, Cartesian2d<CategoryCoord, Y>>,
        categories: &[String],
        y_start: f64,
    ) -> Result<()> {
        let rotate = self.rotate_labels(categories);
        let base = TextStyle::from((FONT, LABEL_SIZE).into_font()).color(&BLACK);
        let style = if rotate {
            base.transform(FontTransform::Rotate270)
                .pos(Pos::new(HPos::Right, VPos::Center))
        } else {
            base.pos(Pos::new(HPos::Center, VPos::Top))
        };
        for (i, label) in categories.iter().enumerate() {
            let (px, py) = chart.backend_coord(&(i as f64, y_start));
            root.draw(&Text::new(
                truncate_label(label),
                (px, py + 8),
                style.clone(),
            ))
            .map_err(draw_error)?;
        }
        Ok(())
    }

    /// Extend the y range so the legend in the upper right corner sits above the data.
    fn legend_headroom(&self, range: Range<f64>) -> Range<f64> {
        let n_series = self.spec.series.len() as u32;
        if n_series < 2 {
            return range;
        }
        let plot_height = self
            .spec
            .style
            .height
            .saturating_sub(MARGIN * 2 + TITLE_AREA + X_LABEL_AREA)
            .max(1) as f64;
        let legend_height = (n_series * LEGEND_ROW + LEGEND_PADDING) as f64;
        let fraction = (legend_height / plot_height).min(0.5);
        let span = range.end - range.start;
        range.start..range.start + span / (1.0 - fraction)
    }

    fn rotate_labels(&self, categories: &[String]) -> bool {
        if categories.is_empty() {
            return false;
        }
        let plot_width = self
            .spec
            .style
            .width
            .saturating_sub(Y_LABEL_AREA + MARGIN * 2);
        let slot = plot_width / categories.len() as u32;
        label_width(categories) + 8 > slot
    }

    fn x_label_area(&self, categories: &[String]) -> u32 {
        if self.rotate_labels(categories) {
            let max = self.spec.style.height / 3;
            (label_width(categories) + AXIS_TITLE_SIZE as u32 * 2 + 8).min(max)
        } else {
            X_LABEL_AREA
        }
    }
}

fn category_coord(count: usize) -> CategoryCoord {
    CategoryCoord { count }
}

/// Data range padded by 5% on each side. `from_zero` anchors bars at zero.
fn padded_range(values: &[f64], from_zero: bool) -> Range<f64> {
    let (mut min, mut max) = finite_range(values).unwrap_or((0.0, 1.0));
    if from_zero {
        min = min.min(0.0);
        max = max.max(0.0);
    }
    if min == max {
        let pad = if min == 0.0 { 1.0 } else { min.abs() * 0.1 };
        return (min - pad)..(max + pad);
    }
    let pad = (max - min) * 0.05;
    let start = if from_zero && min == 0.0 {
        0.0
    } else {
        min - pad
    };
    let end = if from_zero && max == 0.0 {
        0.0
    } else {
        max + pad
    };
    start..end
}

/// Short tick label: plain decimals for typical values, scientific notation
/// for very large or very small ones.
fn format_tick(value: &f64) -> String {
    let abs = value.abs();
    if abs == 0.0 {
        return "0".to_string();
    }
    if !(1e-3..1e6).contains(&abs) {
        return format!("{:.1e}", value);
    }
    let text = format!("{:.3}", value);
    let text = text.trim_end_matches('0').trim_end_matches('.');
    if text == "-0" {
        "0".to_string()
    } else {
        text.to_string()
    }
}

/// Shorten long labels in the middle, keeping both ends.
/// Sample names often differ only at the end.
fn truncate_label(label: &str) -> String {
    let chars: Vec<char> = label.chars().collect();
    if chars.len() <= MAX_LABEL_CHARS {
        return label.to_string();
    }
    let head = (MAX_LABEL_CHARS - 1) / 2;
    let tail = MAX_LABEL_CHARS - 1 - head;
    let start: String = chars[..head].iter().collect();
    let end: String = chars[chars.len() - tail..].iter().collect();
    format!("{}…{}", start, end)
}

fn label_width(labels: &[String]) -> u32 {
    let chars = labels
        .iter()
        .map(|l| l.chars().count().min(MAX_LABEL_CHARS))
        .max()
        .unwrap_or(0);
    chars as u32 * CHAR_WIDTH
}

fn to_rgb(color: Rgb) -> RGBColor {
    RGBColor(color.r, color.g, color.b)
}

fn draw_error<E: std::error::Error + Send + Sync>(error: DrawingAreaErrorKind<E>) -> anyhow::Error {
    anyhow!("Failed drawing plot: {}", error)
}

fn to_png(tree: &usvg::Tree, scale: f32) -> Result<Vec<u8>> {
    let scale = if scale > 0.0 { scale } else { 1.0 };
    let size = tree
        .size()
        .to_int_size()
        .scale_by(scale)
        .context("Failed scaling the PNG output")?;
    let mut pixmap = tiny_skia::Pixmap::new(size.width(), size.height())
        .context("Failed allocating the PNG canvas")?;
    pixmap.fill(tiny_skia::Color::WHITE);
    resvg::render(
        tree,
        tiny_skia::Transform::from_scale(scale, scale),
        &mut pixmap.as_mut(),
    );
    pixmap.encode_png().context("Failed encoding PNG")
}

#[cfg(feature = "pdf")]
fn to_pdf(tree: &usvg::Tree) -> Result<Vec<u8>> {
    svg2pdf::to_pdf(
        tree,
        svg2pdf::ConversionOptions::default(),
        svg2pdf::PageOptions::default(),
    )
    .map_err(|e| anyhow!("Failed converting plot to PDF: {}", e))
}

#[cfg(not(feature = "pdf"))]
fn to_pdf(_tree: &usvg::Tree) -> Result<Vec<u8>> {
    anyhow::bail!("PDF output needs the `pdf` feature of segul-core")
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::plot::color::CvdKind;
    use crate::plot::spec::{BoxStats, HeatmapData, Series};

    fn scatter() -> PlotSpec {
        let mut spec = PlotSpec::new("Scatter", ChartKind::Scatter, "x", "y");
        spec.series
            .push(Series::xy("a", vec![1.0, 2.0, 3.0], vec![3.0, 1.0, 2.0]));
        spec.series
            .push(Series::xy("b", vec![1.0, 2.0], vec![2.0, 2.5]));
        spec
    }

    fn render_svg(spec: &PlotSpec) -> String {
        PlotRenderer::new(spec, &PlotFonts::new()).to_svg().unwrap()
    }

    #[test]
    fn test_render_svg_without_fonts() {
        let spec = scatter();
        let fonts = PlotFonts::new();
        let renderer = PlotRenderer::new(&spec, &fonts);
        let svg = renderer.render_bytes(&PlotFormat::Svg).unwrap();
        assert!(String::from_utf8(svg).unwrap().contains("<svg"));
        assert!(renderer.render_bytes(&PlotFormat::Png).is_err());
    }

    #[test]
    fn test_render_png_and_pdf() {
        let spec = scatter();
        let fonts = PlotFonts::system();
        let renderer = PlotRenderer::new(&spec, &fonts);
        let png = renderer.render_bytes(&PlotFormat::Png).unwrap();
        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
        let pdf = renderer.render_bytes(&PlotFormat::Pdf);
        if cfg!(feature = "pdf") {
            assert_eq!(&pdf.unwrap()[..5], b"%PDF-");
        } else {
            assert!(pdf.is_err());
        }
    }

    #[test]
    fn test_render_to_file() {
        let dir = tempfile::TempDir::with_prefix("plot").unwrap();
        let output = dir.path().join("nested").join("plot.png");
        let fonts = PlotFonts::system();
        PlotRenderer::new(&scatter(), &fonts)
            .render(&output, &PlotFormat::Png)
            .unwrap();
        assert!(output.is_file());
    }

    #[test]
    fn test_render_line_band() {
        let mut spec = PlotSpec::new("Line", ChartKind::Line, "x", "y");
        let mut series = Series::xy("q", vec![1.0, 2.0, 3.0], vec![30.0, 32.0, 31.0]);
        series.lower = Some(vec![20.0, 22.0, 21.0]);
        series.upper = Some(vec![40.0, 40.0, 40.0]);
        spec.series.push(series);
        assert!(render_svg(&spec).contains("<polygon"));
    }

    #[test]
    fn test_render_bars() {
        for kind in [ChartKind::Bar, ChartKind::StackedBar] {
            let mut spec = PlotSpec::new("Bars", kind, "sample", "count");
            spec.categories = vec!["s1".into(), "s2".into()];
            spec.series.push(Series::categorical("a", vec![1.0, 2.0]));
            spec.series
                .push(Series::categorical("b", vec![3.0, f64::NAN]));
            let svg = render_svg(&spec);
            assert!(svg.contains(">\ns1\n<"));
            assert!(svg.contains(">\ns2\n<"));
        }
    }

    #[test]
    fn test_render_histogram() {
        let mut spec = PlotSpec::new("Hist", ChartKind::Histogram, "x", "count");
        spec.series
            .push(Series::histogram("x", &[1.0, 2.0, 2.0, 5.0], 0));
        assert!(render_svg(&spec).contains("<rect"));
    }

    #[test]
    fn test_render_heatmap_with_cvd() {
        let mut spec = PlotSpec::new("Heatmap", ChartKind::Heatmap, "", "");
        spec.heatmap = Some(
            HeatmapData::new(
                vec!["r1".into(), "r2".into()],
                vec!["c1".into(), "c2".into()],
                vec![0.0, 0.5, f64::NAN, 1.0],
            )
            .unwrap(),
        );
        spec.style.cvd = Some(CvdKind::Achromatopsia);
        let svg = render_svg(&spec);
        assert!(svg.contains(">\nr2\n<"));
        assert!(svg.contains(">\n0.5\n<"));
    }

    #[test]
    fn test_render_boxes() {
        let mut spec = PlotSpec::new("Box", ChartKind::BoxPlot, "", "value");
        spec.boxes
            .push(BoxStats::from_values("a", &[1.0, 2.0, 3.0, 4.0, 50.0]).unwrap());
        spec.boxes
            .push(BoxStats::from_values("b", &[2.0, 2.0, 3.0]).unwrap());
        assert!(render_svg(&spec).contains("<circle"));
    }

    #[test]
    fn test_render_invalid_spec() {
        let spec = PlotSpec::new("Empty", ChartKind::Bar, "", "");
        let fonts = PlotFonts::new();
        assert!(PlotRenderer::new(&spec, &fonts).to_svg().is_err());
    }

    #[test]
    fn test_truncate_label() {
        assert_eq!(truncate_label("short"), "short");
        let label = "Genus_species_specimen_001_extra_long_name";
        let short = truncate_label(label);
        assert_eq!(short.chars().count(), MAX_LABEL_CHARS);
        assert!(short.starts_with("Genus_speci"));
        assert!(short.ends_with("long_name"));
    }

    #[test]
    fn test_format_tick() {
        assert_eq!(format_tick(&0.0), "0");
        assert_eq!(format_tick(&0.125), "0.125");
        assert_eq!(format_tick(&20.0), "20");
        assert_eq!(format_tick(&3.5952e-7), "3.6e-7");
        assert_eq!(format_tick(&2_500_000.0), "2.5e6");
    }
}
