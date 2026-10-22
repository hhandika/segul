//! Colorblind-safe palettes and color vision deficiency simulation.
//!
//! Categorical palettes come from Okabe and Ito (2008) and Paul Tol.
//! Continuous colormaps (viridis, cividis, magma, inferno) are perceptually
//! uniform and stay readable for most color vision deficiencies.
//!
//! [`simulate_cvd`] shows how a color looks to someone with protanopia,
//! deuteranopia, tritanopia, or achromatopsia. It uses the Machado et al.
//! (2009) matrices on linear RGB.
use anyhow::bail;
use serde::{Deserialize, Serialize};

const OKABE_ITO: [&str; 8] = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000",
];

const TOL_BRIGHT: [&str; 7] = [
    "#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB",
];

const TOL_MUTED: [&str; 9] = [
    "#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255", "#44AA99", "#999933",
    "#AA4499",
];

// Machado, Oliveira, and Fernandes (2009), severity 1.0.
// https://doi.org/10.1109/TVCG.2009.113
const PROTANOPIA: [[f64; 3]; 3] = [
    [0.152286, 1.052583, -0.204868],
    [0.114503, 0.786281, 0.099216],
    [-0.003882, -0.048116, 1.051998],
];

const DEUTERANOPIA: [[f64; 3]; 3] = [
    [0.367322, 0.860646, -0.227968],
    [0.280085, 0.672501, 0.047413],
    [-0.011820, 0.042940, 0.968881],
];

const TRITANOPIA: [[f64; 3]; 3] = [
    [1.255528, -0.076749, -0.178779],
    [-0.078411, 0.930809, 0.147602],
    [0.004733, 0.691367, 0.303900],
];

// Rec. 709 luminance weights for linear RGB.
const LUMINANCE: [f64; 3] = [0.2126, 0.7152, 0.0722];

/// An sRGB color.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Rgb {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

impl Rgb {
    pub fn new(r: u8, g: u8, b: u8) -> Self {
        Self { r, g, b }
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

    /// Relative luminance in [0, 1]. Used to pick readable text over a fill.
    pub fn luminance(&self) -> f64 {
        let linear = [
            srgb_to_linear(self.r),
            srgb_to_linear(self.g),
            srgb_to_linear(self.b),
        ];
        dot(&LUMINANCE, &linear)
    }
}

/// Categorical palettes. All are distinguishable with common color vision deficiencies.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum Palette {
    #[default]
    OkabeIto,
    TolBright,
    TolMuted,
}

impl Palette {
    pub fn all() -> &'static [Self] {
        &[Self::OkabeIto, Self::TolBright, Self::TolMuted]
    }

    /// The first `n` palette colors, cycling when `n` is larger than the palette.
    pub fn colors(&self, n: usize) -> Vec<Rgb> {
        let hex: &[&str] = match self {
            Self::OkabeIto => &OKABE_ITO,
            Self::TolBright => &TOL_BRIGHT,
            Self::TolMuted => &TOL_MUTED,
        };
        hex.iter()
            .cycle()
            .take(n)
            .map(|h| Rgb::from_hex(h).expect("Failed parsing palette color"))
            .collect()
    }
}

impl std::fmt::Display for Palette {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::OkabeIto => write!(f, "okabe-ito"),
            Self::TolBright => write!(f, "tol-bright"),
            Self::TolMuted => write!(f, "tol-muted"),
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
            _ => Err(format!("{} is not a valid palette", s)),
        }
    }
}

/// Continuous colormaps for heatmaps.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum Colormap {
    #[default]
    Viridis,
    Cividis,
    Magma,
    Inferno,
}

impl Colormap {
    pub fn all() -> &'static [Self] {
        &[Self::Viridis, Self::Cividis, Self::Magma, Self::Inferno]
    }

    /// Color at `t`, clamped to [0, 1].
    pub fn at(&self, t: f64) -> Rgb {
        let gradient = match self {
            Self::Viridis => colorous::VIRIDIS,
            Self::Cividis => colorous::CIVIDIS,
            Self::Magma => colorous::MAGMA,
            Self::Inferno => colorous::INFERNO,
        };
        let t = if t.is_finite() {
            t.clamp(0.0, 1.0)
        } else {
            0.0
        };
        let c = gradient.eval_continuous(t);
        Rgb::new(c.r, c.g, c.b)
    }
}

impl std::fmt::Display for Colormap {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Viridis => write!(f, "viridis"),
            Self::Cividis => write!(f, "cividis"),
            Self::Magma => write!(f, "magma"),
            Self::Inferno => write!(f, "inferno"),
        }
    }
}

impl std::str::FromStr for Colormap {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "viridis" => Ok(Self::Viridis),
            "cividis" => Ok(Self::Cividis),
            "magma" => Ok(Self::Magma),
            "inferno" => Ok(Self::Inferno),
            _ => Err(format!("{} is not a valid colormap", s)),
        }
    }
}

/// Color vision deficiencies that can be simulated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CvdKind {
    /// Missing long-wavelength (red) cones.
    Protanopia,
    /// Missing medium-wavelength (green) cones.
    Deuteranopia,
    /// Missing short-wavelength (blue) cones.
    Tritanopia,
    /// No color vision (luminance only).
    Achromatopsia,
}

impl CvdKind {
    pub fn all() -> &'static [Self] {
        &[
            Self::Protanopia,
            Self::Deuteranopia,
            Self::Tritanopia,
            Self::Achromatopsia,
        ]
    }
}

impl std::fmt::Display for CvdKind {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Protanopia => write!(f, "protanopia"),
            Self::Deuteranopia => write!(f, "deuteranopia"),
            Self::Tritanopia => write!(f, "tritanopia"),
            Self::Achromatopsia => write!(f, "achromatopsia"),
        }
    }
}

impl std::str::FromStr for CvdKind {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "protanopia" => Ok(Self::Protanopia),
            "deuteranopia" => Ok(Self::Deuteranopia),
            "tritanopia" => Ok(Self::Tritanopia),
            "achromatopsia" => Ok(Self::Achromatopsia),
            _ => Err(format!("{} is not a valid color vision deficiency", s)),
        }
    }
}

/// Simulate how `color` looks with a color vision deficiency.
///
/// `severity` ranges from 0 (normal vision) to 1 (full deficiency).
/// Partial severities interpolate between the original and the fully
/// simulated color in linear RGB, which closely tracks Machado's
/// per-severity matrices.
pub fn simulate_cvd(color: Rgb, kind: &CvdKind, severity: f32) -> Rgb {
    let linear = [
        srgb_to_linear(color.r),
        srgb_to_linear(color.g),
        srgb_to_linear(color.b),
    ];
    let simulated = match kind {
        CvdKind::Protanopia => apply_matrix(&PROTANOPIA, &linear),
        CvdKind::Deuteranopia => apply_matrix(&DEUTERANOPIA, &linear),
        CvdKind::Tritanopia => apply_matrix(&TRITANOPIA, &linear),
        CvdKind::Achromatopsia => [dot(&LUMINANCE, &linear); 3],
    };
    let severity = severity.clamp(0.0, 1.0) as f64;
    let mix = |i: usize| linear[i] + (simulated[i] - linear[i]) * severity;
    Rgb::new(
        linear_to_srgb(mix(0)),
        linear_to_srgb(mix(1)),
        linear_to_srgb(mix(2)),
    )
}

fn apply_matrix(matrix: &[[f64; 3]; 3], rgb: &[f64; 3]) -> [f64; 3] {
    [
        dot(&matrix[0], rgb),
        dot(&matrix[1], rgb),
        dot(&matrix[2], rgb),
    ]
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
    fn test_palette_cycles() {
        let colors = Palette::OkabeIto.colors(10);
        assert_eq!(colors.len(), 10);
        assert_eq!(colors[0], colors[8]);
    }

    #[test]
    fn test_colormap_ends() {
        assert_eq!(Colormap::Viridis.at(0.0), Rgb::new(68, 1, 84));
        assert_eq!(Colormap::Viridis.at(2.0), Colormap::Viridis.at(1.0));
    }

    #[test]
    fn test_simulate_protanopia() {
        let red = Rgb::new(255, 0, 0);
        assert_eq!(
            simulate_cvd(red, &CvdKind::Protanopia, 1.0),
            Rgb::new(109, 95, 0)
        );
        let vermillion = Rgb::new(213, 94, 0);
        assert_eq!(
            simulate_cvd(vermillion, &CvdKind::Protanopia, 1.0),
            Rgb::new(129, 113, 0)
        );
    }

    #[test]
    fn test_simulate_deuteranopia() {
        let green = Rgb::new(0, 158, 115);
        assert_eq!(
            simulate_cvd(green, &CvdKind::Deuteranopia, 1.0),
            Rgb::new(138, 134, 118)
        );
    }

    #[test]
    fn test_simulate_achromatopsia_is_gray() {
        let color = simulate_cvd(Rgb::new(230, 159, 0), &CvdKind::Achromatopsia, 1.0);
        assert_eq!(color.r, color.g);
        assert_eq!(color.g, color.b);
    }

    #[test]
    fn test_simulate_zero_severity() {
        let color = Rgb::new(204, 121, 167);
        assert_eq!(simulate_cvd(color, &CvdKind::Tritanopia, 0.0), color);
    }

    #[test]
    fn test_parse_names() {
        assert_eq!("tol-muted".parse::<Palette>().unwrap(), Palette::TolMuted);
        assert_eq!("cividis".parse::<Colormap>().unwrap(), Colormap::Cividis);
        assert_eq!(
            "deuteranopia".parse::<CvdKind>().unwrap(),
            CvdKind::Deuteranopia
        );
        assert!("rainbow".parse::<Colormap>().is_err());
    }
}
