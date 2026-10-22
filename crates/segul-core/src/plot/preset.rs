//! Shared helpers for the plot presets.
#[cfg(feature = "maf")]
use std::path::Path;

use super::spec::{ChartKind, PlotSpec, PlotStyle, Series};

/// Declare a preset enum with its command-line names and descriptions.
macro_rules! plot_presets {
    (
        $(#[$meta:meta])*
        $name:ident {
            $($(#[$vmeta:meta])* $variant:ident => ($id:literal, $desc:literal)),+ $(,)?
        }
    ) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
        pub enum $name {
            $($(#[$vmeta])* $variant),+
        }

        impl $crate::plot::PlotPreset for $name {
            fn all() -> &'static [Self] {
                &[$(Self::$variant),+]
            }

            fn name(&self) -> &'static str {
                match self {
                    $(Self::$variant => $id),+
                }
            }

            fn description(&self) -> &'static str {
                match self {
                    $(Self::$variant => $desc),+
                }
            }
        }

        impl std::fmt::Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(f, "{}", $crate::plot::PlotPreset::name(self))
            }
        }

        impl std::str::FromStr for $name {
            type Err = String;

            fn from_str(s: &str) -> Result<Self, Self::Err> {
                use $crate::plot::PlotPreset;
                Self::all()
                    .iter()
                    .find(|p| p.name() == s)
                    .copied()
                    .ok_or_else(|| {
                        let names: Vec<&str> = Self::all().iter().map(|p| p.name()).collect();
                        format!("{} is not a valid chart. Available charts: {}", s, names.join(", "))
                    })
            }
        }
    };
}

pub(super) use plot_presets;

/// File stem of a path stored in a summary file, used as a short label.
#[cfg(feature = "maf")]
pub(super) fn file_label(path: &str) -> String {
    Path::new(path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or(path)
        .to_string()
}

#[cfg(feature = "align")]
pub(super) fn histogram_spec(
    title: &str,
    x_label: &str,
    values: &[f64],
    style: &PlotStyle,
) -> PlotSpec {
    let mut spec = PlotSpec::new(title, ChartKind::Histogram, x_label, "Frequency");
    spec.series.push(Series::histogram(x_label, values, 0));
    spec.style = style.clone();
    spec
}

#[cfg(any(feature = "align", feature = "contig", feature = "read"))]
pub(super) fn scatter_spec(
    title: &str,
    axes: (&str, &str),
    x: Vec<f64>,
    y: Vec<f64>,
    style: &PlotStyle,
) -> PlotSpec {
    let mut spec = PlotSpec::new(title, ChartKind::Scatter, axes.0, axes.1);
    spec.series.push(Series::xy(title, x, y));
    spec.style = style.clone();
    spec
}

pub(super) fn bar_spec(
    title: &str,
    axes: (&str, &str),
    categories: Vec<String>,
    series: Vec<Series>,
    style: &PlotStyle,
) -> PlotSpec {
    let mut spec = PlotSpec::new(title, ChartKind::Bar, axes.0, axes.1);
    spec.categories = categories;
    spec.series = series;
    spec.style = style.clone();
    spec
}
