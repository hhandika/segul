//! Charts for segul summary statistics.
//!
//! Plots are built in two steps:
//! 1. A preset reads a summary CSV written by segul and returns a
//!    [`PlotSpec`](spec::PlotSpec), a renderer-agnostic description of the chart.
//! 2. [`PlotRenderer`](render::PlotRenderer) draws the spec as SVG, PNG, or
//!    PDF (`pdf` feature). PNG and PDF use the fonts in
//!    [`PlotFonts`](fonts::PlotFonts), such as a font bundled with an app.
//!
//! Front ends that draw their own interactive charts (such as the segui app)
//! can use the spec directly and call the renderer only for export.
//! Default colors are colorblind-safe, and [`color::simulate_cvd`] previews
//! a plot as seen with a color vision deficiency.
//!
//! ## Example
//! ```rust
//! use std::path::Path;
//!
//! use segul_core::plot::PlotPreset;
//! use segul_core::plot::align::{AlignPlot, AlignmentPlot};
//! use segul_core::plot::fonts::PlotFonts;
//! use segul_core::plot::render::PlotRenderer;
//! use segul_core::plot::spec::{PlotFormat, PlotStyle};
//!
//! let input = Path::new("tests/files/plot/locus_summary.csv");
//! let style = PlotStyle::default();
//! let spec = AlignmentPlot::new(input, &style)
//!     .build(&AlignPlot::PisVsLength)
//!     .unwrap();
//! let fonts = PlotFonts::new();
//! let svg = PlotRenderer::new(&spec, &fonts)
//!     .render_bytes(&PlotFormat::Svg)
//!     .unwrap();
//! assert!(!svg.is_empty());
//! assert!(AlignPlot::all().contains(&"completeness".parse().unwrap()));
//! ```
pub mod batch;
pub mod color;
pub mod fonts;
pub mod render;
pub mod spec;

#[cfg(feature = "align")]
pub mod align;
#[cfg(feature = "contig")]
pub mod contig;
#[cfg(feature = "align")]
pub mod distance;
#[cfg(feature = "maf")]
pub mod maf;
#[cfg(feature = "read")]
pub mod read;

#[cfg(any(
    feature = "align",
    feature = "contig",
    feature = "maf",
    feature = "read"
))]
mod preset;

/// Named charts built from one kind of segul summary file.
pub trait PlotPreset: Sized + Copy + 'static {
    /// Every preset, in display order.
    fn all() -> &'static [Self];

    /// Name used on the command line, e.g. `pis-vs-length`.
    fn name(&self) -> &'static str;

    /// One-line description of the chart and its input file.
    fn description(&self) -> &'static str;
}
