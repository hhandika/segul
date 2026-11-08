//! Color library shared by segul plots, the terminal viewers, and apps
//! built on segul-core (such as the segui GUI).
//!
//! Every color segul draws comes from this module, so charts, tables, and
//! alignment viewers look the same in every front end:
//!
//! - [`Palette`]: categorical palettes for chart series. ColorBrewer
//!   qualitative schemes (default Dark2), Okabe-Ito, and Paul Tol.
//! - [`Colormap`]: continuous colormaps for heatmaps. Viridis-family maps
//!   (default viridis) and ColorBrewer sequential and diverging schemes.
//! - [`residue::ResiduePalette`]: residue colors for alignment viewers.
//! - [`ui::ChartTokens`]: chart background, axis, grid, text, and status colors.
//! - WCAG 2.2 helpers ([`contrast_ratio`], [`ensure_contrast`],
//!   [`best_text_on`]) that every renderer uses so text reaches 4.5:1 and
//!   lines and marks reach 3:1 against their background.
//! - [`simulate_cvd`]: previews colors as seen with a color vision deficiency.
//!
//! ## Sources
//!
//! - Okabe, M. and Ito, K. (2008). Color Universal Design (CUD): How to make
//!   figures and presentations that are friendly to colorblind people.
//!   <https://jfly.uni-koeln.de/color/>
//! - Tol, P. (2021). Colour schemes. SRON technical note SRON/EPS/TN/09-002.
//!   <https://personal.sron.nl/~pault/>
//! - Brewer, C. A. ColorBrewer 2.0. <https://colorbrewer2.org>. Apache-2.0.
//!   Harrower, M. and Brewer, C. A. (2003). ColorBrewer.org: An online tool
//!   for selecting colour schemes for maps. The Cartographic Journal 40:27-37.
//!   <https://doi.org/10.1179/000870403235002042>
//! - van der Walt, S. and Smith, N. (2015). viridis, magma, inferno, and
//!   plasma colormaps for matplotlib. CC0. <https://bids.github.io/colormap/>
//! - Nuñez, J. R., Anderton, C. R., and Renslow, R. S. (2018). Optimizing
//!   colormaps with consideration for color vision deficiency to enable
//!   accurate interpretation of scientific data. PLOS ONE 13:e0199239.
//!   <https://doi.org/10.1371/journal.pone.0199239>
//! - Machado, G. M., Oliveira, M. M., and Fernandes, L. A. F. (2009). A
//!   physiologically-based model for simulation of color vision deficiency.
//!   IEEE TVCG 15:1291-1298. <https://doi.org/10.1109/TVCG.2009.113>
//! - Thompson, J. D. et al. (1997). The CLUSTAL_X windows interface. Nucleic
//!   Acids Research 25:4876-4882. <https://doi.org/10.1093/nar/25.24.4876>
//! - W3C (2023). Web Content Accessibility Guidelines (WCAG) 2.2, success
//!   criteria 1.4.3 (contrast, minimum) and 1.4.11 (non-text contrast).
//!   <https://www.w3.org/TR/WCAG22/>
//!
//! Gradient and ColorBrewer data are read from the `colorous` crate
//! (<https://github.com/dtolnay/colorous>), which takes them from d3-scale-chromatic.
//!
//! ## Example
//! ```rust
//! use segul_core::color::{Palette, Rgb, contrast_ratio, ensure_contrast, NON_TEXT};
//!
//! let white = Rgb::new(255, 255, 255);
//! let yellow = Palette::OkabeIto.colors(7)[6];
//! assert!(contrast_ratio(yellow, white) < NON_TEXT);
//! let line = ensure_contrast(yellow, white, NON_TEXT);
//! assert!(contrast_ratio(line, white) >= NON_TEXT);
//! ```
use anyhow::bail;
use serde::{Deserialize, Serialize};

mod colormap;
mod contrast;
mod cvd;
mod palette;
pub mod residue;
pub mod ui;

pub use colormap::{Colormap, ColormapKind};
pub use contrast::{
    AA_LARGE_TEXT, AA_TEXT, NON_TEXT, best_text_on, blend, contrast_ratio, ensure_contrast,
    relative_luminance,
};
pub use cvd::{CvdKind, simulate_cvd};
pub use palette::Palette;

// Rec. 709 luminance weights for linear RGB, as used by WCAG.
const LUMINANCE: [f64; 3] = [0.2126, 0.7152, 0.0722];

/// An sRGB color.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Rgb {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

impl Rgb {
    pub const BLACK: Self = Self::new(0, 0, 0);
    pub const WHITE: Self = Self::new(255, 255, 255);

    pub const fn new(r: u8, g: u8, b: u8) -> Self {
        Self { r, g, b }
    }

    /// Build from `0xRRGGBB`.
    pub const fn from_u32(rgb: u32) -> Self {
        Self::new((rgb >> 16) as u8, (rgb >> 8) as u8, rgb as u8)
    }

    /// Parse a `#RRGGBB` hex string.
    pub fn from_hex(hex: &str) -> anyhow::Result<Self> {
        let hex = hex.trim_start_matches('#');
        if hex.len() != 6 {
            bail!("Invalid hex color: {}", hex);
        }
        let channel = |i: usize| u8::from_str_radix(&hex[i..i + 2], 16);
        Ok(Self::new(channel(0)?, channel(2)?, channel(4)?))
    }

    /// Format as `#RRGGBB`.
    pub fn to_hex(&self) -> String {
        format!("#{:02X}{:02X}{:02X}", self.r, self.g, self.b)
    }

    /// Opaque `0xAARRGGBB`, the layout Flutter's `Color` takes.
    pub fn to_argb(&self) -> u32 {
        0xFF00_0000 | (self.r as u32) << 16 | (self.g as u32) << 8 | self.b as u32
    }

    /// Read `0xAARRGGBB`, ignoring alpha.
    pub fn from_argb(argb: u32) -> Self {
        Self::from_u32(argb & 0x00FF_FFFF)
    }

    /// WCAG relative luminance in [0, 1]. Same as [`relative_luminance`].
    pub fn luminance(&self) -> f64 {
        relative_luminance(*self)
    }

    fn to_linear(self) -> [f64; 3] {
        [
            srgb_to_linear(self.r),
            srgb_to_linear(self.g),
            srgb_to_linear(self.b),
        ]
    }
}

/// Light or dark appearance. Residue and chart colors have a value for each.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Theme {
    #[default]
    Light,
    Dark,
}

impl Theme {
    pub fn all() -> &'static [Self] {
        &[Self::Light, Self::Dark]
    }

    // Picks the light or dark value of a `[light, dark]` pair.
    fn pick(&self, pair: [u32; 2]) -> Rgb {
        match self {
            Self::Light => Rgb::from_u32(pair[0]),
            Self::Dark => Rgb::from_u32(pair[1]),
        }
    }
}

impl std::fmt::Display for Theme {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Light => write!(f, "light"),
            Self::Dark => write!(f, "dark"),
        }
    }
}

impl std::str::FromStr for Theme {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "light" => Ok(Self::Light),
            "dark" => Ok(Self::Dark),
            _ => Err(format!("{} is not a valid theme", s)),
        }
    }
}

fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn srgb_to_linear(channel: u8) -> f64 {
    let c = channel as f64 / 255.0;
    if c <= 0.04045 {
        c / 12.92
    } else {
        ((c + 0.055) / 1.055).powf(2.4)
    }
}

fn linear_to_srgb(channel: f64) -> u8 {
    let c = channel.clamp(0.0, 1.0);
    let encoded = if c <= 0.0031308 {
        c * 12.92
    } else {
        1.055 * c.powf(1.0 / 2.4) - 0.055
    };
    (encoded * 255.0).round() as u8
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_hex_round_trip() {
        let color = Rgb::from_hex("#56B4E9").unwrap();
        assert_eq!(color, Rgb::new(86, 180, 233));
        assert_eq!(color.to_hex(), "#56B4E9");
        assert!(Rgb::from_hex("#FFF").is_err());
    }

    #[test]
    fn test_argb_round_trip() {
        let color = Rgb::from_u32(0x0072B2);
        assert_eq!(color.to_argb(), 0xFF0072B2);
        assert_eq!(Rgb::from_argb(0x800072B2), color);
    }

    #[test]
    fn test_theme_names() {
        for theme in Theme::all() {
            assert_eq!(theme.to_string().parse::<Theme>().unwrap(), *theme);
        }
    }
}
