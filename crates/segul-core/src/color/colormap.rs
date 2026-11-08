//! Continuous colormaps for heatmaps.
use serde::{Deserialize, Serialize};

use super::Rgb;

/// Whether a colormap runs from low to high or diverges from a midpoint.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ColormapKind {
    Sequential,
    Diverging,
}

/// Continuous colormaps.
///
/// Viridis, magma, inferno, and plasma are from van der Walt and Smith
/// (2015, <https://bids.github.io/colormap/>, CC0). Cividis is from Nuñez
/// et al. (2018, <https://doi.org/10.1371/journal.pone.0199239>). These are
/// perceptually uniform and readable with color vision deficiencies.
///
/// The ColorBrewer schemes (<https://colorbrewer2.org>, Apache-2.0) are
/// interpolated from their nine- or eleven-class versions, and all seven
/// are marked colorblind safe there. The two ends of a diverging map have
/// similar lightness, so they look alike in grayscale.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Colormap {
    #[default]
    Viridis,
    Cividis,
    Magma,
    Inferno,
    Plasma,
    /// ColorBrewer Blues.
    Blues,
    /// ColorBrewer Greens.
    Greens,
    /// ColorBrewer YlGnBu.
    #[serde(rename = "ylgnbu")]
    YlGnBu,
    /// ColorBrewer YlOrRd.
    #[serde(rename = "ylorrd")]
    YlOrRd,
    /// ColorBrewer RdBu, diverging.
    #[serde(rename = "rdbu")]
    RdBu,
    /// ColorBrewer PuOr, diverging.
    #[serde(rename = "puor")]
    PuOr,
    /// ColorBrewer BrBG, diverging.
    #[serde(rename = "brbg")]
    BrBG,
}

impl Colormap {
    pub fn all() -> &'static [Self] {
        &[
            Self::Viridis,
            Self::Cividis,
            Self::Magma,
            Self::Inferno,
            Self::Plasma,
            Self::Blues,
            Self::Greens,
            Self::YlGnBu,
            Self::YlOrRd,
            Self::RdBu,
            Self::PuOr,
            Self::BrBG,
        ]
    }

    /// Name for menus and legends.
    pub fn label(&self) -> &'static str {
        match self {
            Self::Viridis => "Viridis",
            Self::Cividis => "Cividis",
            Self::Magma => "Magma",
            Self::Inferno => "Inferno",
            Self::Plasma => "Plasma",
            Self::Blues => "ColorBrewer Blues",
            Self::Greens => "ColorBrewer Greens",
            Self::YlGnBu => "ColorBrewer YlGnBu",
            Self::YlOrRd => "ColorBrewer YlOrRd",
            Self::RdBu => "ColorBrewer RdBu",
            Self::PuOr => "ColorBrewer PuOr",
            Self::BrBG => "ColorBrewer BrBG",
        }
    }

    pub fn kind(&self) -> ColormapKind {
        match self {
            Self::RdBu | Self::PuOr | Self::BrBG => ColormapKind::Diverging,
            _ => ColormapKind::Sequential,
        }
    }

    /// Color at `t`, clamped to [0, 1].
    pub fn at(&self, t: f64) -> Rgb {
        let t = if t.is_finite() {
            t.clamp(0.0, 1.0)
        } else {
            0.0
        };
        let c = self.gradient().eval_continuous(t);
        Rgb::new(c.r, c.g, c.b)
    }

    /// `n` evenly spaced colors from start to end (at least two).
    pub fn stops(&self, n: usize) -> Vec<Rgb> {
        let n = n.max(2);
        (0..n).map(|i| self.at(i as f64 / (n - 1) as f64)).collect()
    }

    fn gradient(&self) -> colorous::Gradient {
        match self {
            Self::Viridis => colorous::VIRIDIS,
            Self::Cividis => colorous::CIVIDIS,
            Self::Magma => colorous::MAGMA,
            Self::Inferno => colorous::INFERNO,
            Self::Plasma => colorous::PLASMA,
            Self::Blues => colorous::BLUES,
            Self::Greens => colorous::GREENS,
            Self::YlGnBu => colorous::YELLOW_GREEN_BLUE,
            Self::YlOrRd => colorous::YELLOW_ORANGE_RED,
            Self::RdBu => colorous::RED_BLUE,
            Self::PuOr => colorous::PURPLE_ORANGE,
            Self::BrBG => colorous::BROWN_GREEN,
        }
    }
}

impl std::fmt::Display for Colormap {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Viridis => write!(f, "viridis"),
            Self::Cividis => write!(f, "cividis"),
            Self::Magma => write!(f, "magma"),
            Self::Inferno => write!(f, "inferno"),
            Self::Plasma => write!(f, "plasma"),
            Self::Blues => write!(f, "blues"),
            Self::Greens => write!(f, "greens"),
            Self::YlGnBu => write!(f, "ylgnbu"),
            Self::YlOrRd => write!(f, "ylorrd"),
            Self::RdBu => write!(f, "rdbu"),
            Self::PuOr => write!(f, "puor"),
            Self::BrBG => write!(f, "brbg"),
        }
    }
}

impl std::str::FromStr for Colormap {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::all()
            .iter()
            .find(|c| c.to_string() == s)
            .copied()
            .ok_or_else(|| format!("{} is not a valid colormap", s))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_colormap_ends() {
        assert_eq!(Colormap::Viridis.at(0.0), Rgb::new(68, 1, 84));
        assert_eq!(Colormap::Viridis.at(2.0), Colormap::Viridis.at(1.0));
        assert_eq!(Colormap::Viridis.at(f64::NAN), Colormap::Viridis.at(0.0));
    }

    #[test]
    fn test_stops() {
        let stops = Colormap::Blues.stops(5);
        assert_eq!(stops.len(), 5);
        assert_eq!(stops[0], Colormap::Blues.at(0.0));
        assert_eq!(stops[4], Colormap::Blues.at(1.0));
        assert_eq!(Colormap::Magma.stops(0).len(), 2);
    }

    #[test]
    fn test_kind() {
        assert_eq!(Colormap::RdBu.kind(), ColormapKind::Diverging);
        assert_eq!(Colormap::Viridis.kind(), ColormapKind::Sequential);
    }

    #[test]
    fn test_parse_and_serde_names() {
        for colormap in Colormap::all() {
            assert_eq!(colormap.to_string().parse::<Colormap>().unwrap(), *colormap);
            let json = serde_json::to_string(colormap).unwrap();
            assert_eq!(json, format!("\"{}\"", colormap));
        }
        assert!("rainbow".parse::<Colormap>().is_err());
    }
}
