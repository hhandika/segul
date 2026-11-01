//! Color vision deficiency simulation.
//!
//! Uses the Machado, Oliveira, and Fernandes (2009) matrices at severity 1.0
//! on linear RGB. Source: <https://doi.org/10.1109/TVCG.2009.113>, matrices
//! from <https://www.inf.ufrgs.br/~oliveira/pubs_files/CVD_Simulation/CVD_Simulation.html>.
use serde::{Deserialize, Serialize};

use super::{LUMINANCE, Rgb, dot, linear_to_srgb};

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

/// Color vision deficiencies that can be simulated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
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

    /// Name for menus and legends.
    pub fn label(&self) -> &'static str {
        match self {
            Self::Protanopia => "Protanopia",
            Self::Deuteranopia => "Deuteranopia",
            Self::Tritanopia => "Tritanopia",
            Self::Achromatopsia => "Achromatopsia",
        }
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
    let linear = color.to_linear();
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

#[cfg(test)]
mod test {
    use super::*;

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
    fn test_parse_and_serde_names() {
        for kind in CvdKind::all() {
            assert_eq!(kind.to_string().parse::<CvdKind>().unwrap(), *kind);
            let json = serde_json::to_string(kind).unwrap();
            assert_eq!(json, format!("\"{}\"", kind));
        }
    }
}
