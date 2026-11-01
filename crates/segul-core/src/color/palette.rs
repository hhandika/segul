//! Categorical palettes for chart series.
use serde::{Deserialize, Serialize};

use super::Rgb;

/// Okabe and Ito (2008), <https://jfly.uni-koeln.de/color/>.
///
/// Reordered to start with the strongest colors, so charts with one or two
/// series get the most contrast. Black is replaced with gray so the last
/// color shows on dark backgrounds too. Same order as the segui app.
const OKABE_ITO: [u32; 8] = [
    0x0072B2, // blue
    0xE69F00, // orange
    0x009E73, // bluish green
    0xCC79A7, // reddish purple
    0x56B4E9, // sky blue
    0xD55E00, // vermillion
    0xF0E442, // yellow
    0x999999, // gray
];

/// Paul Tol's bright scheme, <https://personal.sron.nl/~pault/#sec:qualitative>.
const TOL_BRIGHT: [u32; 7] = [
    0x4477AA, 0xEE6677, 0x228833, 0xCCBB44, 0x66CCEE, 0xAA3377, 0xBBBBBB,
];

/// Paul Tol's muted scheme, <https://personal.sron.nl/~pault/#sec:qualitative>.
const TOL_MUTED: [u32; 9] = [
    0xCC6677, 0x332288, 0xDDCC77, 0x117733, 0x88CCEE, 0x882255, 0x44AA99, 0x999933, 0xAA4499,
];

/// Categorical palettes.
///
/// Okabe-Ito and the Paul Tol schemes are designed to stay distinct with
/// common color vision deficiencies. The ColorBrewer schemes
/// (<https://colorbrewer2.org>, Apache-2.0) are only colorblind safe for
/// their first few colors; see [`Palette::colorblind_safe_classes`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Palette {
    #[default]
    OkabeIto,
    TolBright,
    TolMuted,
    /// ColorBrewer Set2.
    Set2,
    /// ColorBrewer Dark2.
    Dark2,
    /// ColorBrewer Paired.
    Paired,
    /// ColorBrewer Set1. Not colorblind safe.
    Set1,
}

impl Palette {
    pub fn all() -> &'static [Self] {
        &[
            Self::OkabeIto,
            Self::TolBright,
            Self::TolMuted,
            Self::Set2,
            Self::Dark2,
            Self::Paired,
            Self::Set1,
        ]
    }

    /// Name for menus and legends.
    pub fn label(&self) -> &'static str {
        match self {
            Self::OkabeIto => "Okabe-Ito",
            Self::TolBright => "Tol bright",
            Self::TolMuted => "Tol muted",
            Self::Set2 => "ColorBrewer Set2",
            Self::Dark2 => "ColorBrewer Dark2",
            Self::Paired => "ColorBrewer Paired",
            Self::Set1 => "ColorBrewer Set1",
        }
    }

    /// Number of distinct colors before the palette repeats.
    pub fn len(&self) -> usize {
        self.base().len()
    }

    /// Always false: every palette has colors.
    pub fn is_empty(&self) -> bool {
        false
    }

    /// How many colors, from the start, stay distinct with color vision
    /// deficiencies. Tol and Okabe-Ito schemes are safe in full. ColorBrewer
    /// counts follow the colorblind filter on colorbrewer2.org.
    pub fn colorblind_safe_classes(&self) -> usize {
        match self {
            Self::OkabeIto | Self::TolBright | Self::TolMuted => self.len(),
            Self::Set2 | Self::Dark2 => 3,
            Self::Paired => 4,
            Self::Set1 => 0,
        }
    }

    /// Whether the first `n` colors stay distinct with color vision deficiencies.
    pub fn is_colorblind_safe(&self, n: usize) -> bool {
        n <= self.colorblind_safe_classes()
    }

    /// The first `n` palette colors, cycling when `n` is larger than the palette.
    pub fn colors(&self, n: usize) -> Vec<Rgb> {
        self.base().into_iter().cycle().take(n).collect()
    }

    fn base(&self) -> Vec<Rgb> {
        let from_hex = |hex: &[u32]| hex.iter().map(|&h| Rgb::from_u32(h)).collect();
        let from_brewer =
            |colors: &[colorous::Color]| colors.iter().map(|c| Rgb::new(c.r, c.g, c.b)).collect();
        match self {
            Self::OkabeIto => from_hex(&OKABE_ITO),
            Self::TolBright => from_hex(&TOL_BRIGHT),
            Self::TolMuted => from_hex(&TOL_MUTED),
            Self::Set2 => from_brewer(&colorous::SET2),
            Self::Dark2 => from_brewer(&colorous::DARK2),
            Self::Paired => from_brewer(&colorous::PAIRED),
            Self::Set1 => from_brewer(&colorous::SET1),
        }
    }
}

impl std::fmt::Display for Palette {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::OkabeIto => write!(f, "okabe-ito"),
            Self::TolBright => write!(f, "tol-bright"),
            Self::TolMuted => write!(f, "tol-muted"),
            Self::Set2 => write!(f, "set2"),
            Self::Dark2 => write!(f, "dark2"),
            Self::Paired => write!(f, "paired"),
            Self::Set1 => write!(f, "set1"),
        }
    }
}

impl std::str::FromStr for Palette {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "okabe-ito" => Ok(Self::OkabeIto),
            "tol-bright" => Ok(Self::TolBright),
            "tol-muted" => Ok(Self::TolMuted),
            "set2" => Ok(Self::Set2),
            "dark2" => Ok(Self::Dark2),
            "paired" => Ok(Self::Paired),
            "set1" => Ok(Self::Set1),
            _ => Err(format!("{} is not a valid palette", s)),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_palette_cycles() {
        let colors = Palette::OkabeIto.colors(10);
        assert_eq!(colors.len(), 10);
        assert_eq!(colors[0], colors[8]);
        assert_eq!(colors[0].to_hex(), "#0072B2");
    }

    #[test]
    fn test_brewer_colors() {
        // First Dark2 color on colorbrewer2.org.
        assert_eq!(Palette::Dark2.colors(1)[0].to_hex(), "#1B9E77");
        assert_eq!(Palette::Paired.len(), 12);
    }

    #[test]
    fn test_colorblind_safe() {
        assert!(Palette::OkabeIto.is_colorblind_safe(8));
        assert!(Palette::Dark2.is_colorblind_safe(3));
        assert!(!Palette::Dark2.is_colorblind_safe(4));
        assert!(!Palette::Set1.is_colorblind_safe(1));
    }

    #[test]
    fn test_parse_and_serde_names() {
        for palette in Palette::all() {
            assert_eq!(palette.to_string().parse::<Palette>().unwrap(), *palette);
            let json = serde_json::to_string(palette).unwrap();
            assert_eq!(json, format!("\"{}\"", palette));
        }
        assert!("rainbow".parse::<Palette>().is_err());
    }
}
