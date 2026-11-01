//! Residue colors, matching the SEGUI alignment viewer.

use ratatui::style::Color;

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum Theme {
    Dark,
    Light,
}

/// How residues are colored.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum ColorMode {
    /// Colored by residue.
    Base,
    /// One color for residues, another for gaps.
    Mono,
    /// Residues equal to the first sequence are shown as dots.
    Reference,
}

impl ColorMode {
    pub(super) fn next(self) -> Self {
        match self {
            ColorMode::Base => ColorMode::Mono,
            ColorMode::Mono => ColorMode::Reference,
            ColorMode::Reference => ColorMode::Base,
        }
    }

    pub(super) fn label(self) -> &'static str {
        match self {
            ColorMode::Base => "Base",
            ColorMode::Mono => "Mono",
            ColorMode::Reference => "Diff to first",
        }
    }
}

const DNA_GROUPS: &[(&str, &[u8], [u32; 2])] = &[
    ("A", b"A", [0x97C459, 0x3B6D11]),
    ("C", b"C", [0x85B7EB, 0x185FA5]),
    ("G", b"G", [0xFAC775, 0x854F0B]),
    ("T", b"TU", [0xF09595, 0xA32D2D]),
];

// Clustal X groups.
const AA_GROUPS: &[(&str, &[u8], [u32; 2])] = &[
    ("Hydrophobic", b"AILMFWV", [0x85B7EB, 0x185FA5]),
    ("Positive", b"KR", [0xF09595, 0xA32D2D]),
    ("Negative", b"DE", [0xED93B1, 0x993556]),
    ("Polar", b"NQST", [0x97C459, 0x3B6D11]),
    ("Aromatic", b"HY", [0x5DCAA5, 0x0F6E56]),
    ("G", b"G", [0xFAC775, 0x854F0B]),
    ("P", b"P", [0xFAEEDA, 0x633806]),
    ("C", b"C", [0xF5C4B3, 0x993C1D]),
];

// Light and dark values of the colors shared by both palettes.
const FALLBACK: [u32; 2] = [0xB4B2A9, 0x5F5E5A];
const GAP: [u32; 2] = [0xE6E4DC, 0x2C2C2A];
const ON_RESIDUE: [u32; 2] = [0x2C2C2A, 0xF1EFE8];
const MONO: [u32; 2] = [0xD3D1C7, 0x444441];
const IDENTITY: [u32; 2] = [0xF1EFE8, 0x242422];
const ON_IDENTITY: [u32; 2] = [0x888780, 0x8A8984];
const ACCENT: [u32; 2] = [0x0F6E56, 0x5DCAA5];
const TERTIARY: [u32; 2] = [0x993556, 0xED93B1];

/// Colors for one theme and data type.
pub(super) struct Palette {
    residues: [Option<Color>; 128],
    pub(super) legend: Vec<(&'static str, Color)>,
    pub(super) fallback: Color,
    pub(super) gap: Color,
    pub(super) on_residue: Color,
    pub(super) mono: Color,
    pub(super) identity: Color,
    pub(super) on_identity: Color,
    pub(super) accent: Color,
    pub(super) tertiary: Color,
    is_truecolor: bool,
}

impl Palette {
    pub(super) fn new(theme: Theme, is_dna: bool, is_truecolor: bool) -> Self {
        let pick = |hex: [u32; 2]| -> Color {
            let rgb = match theme {
                Theme::Light => hex[0],
                Theme::Dark => hex[1],
            };
            to_color(rgb, is_truecolor)
        };
        let groups = if is_dna { DNA_GROUPS } else { AA_GROUPS };
        let mut residues = [None; 128];
        let mut legend = Vec::with_capacity(groups.len());
        for (label, codes, hex) in groups {
            let color = pick(*hex);
            codes
                .iter()
                .for_each(|&c| residues[c as usize] = Some(color));
            legend.push((*label, color));
        }
        Self {
            residues,
            legend,
            fallback: pick(FALLBACK),
            gap: pick(GAP),
            on_residue: pick(ON_RESIDUE),
            mono: pick(MONO),
            identity: pick(IDENTITY),
            on_identity: pick(ON_IDENTITY),
            accent: pick(ACCENT),
            tertiary: pick(TERTIARY),
            is_truecolor,
        }
    }

    /// Background of a residue in base mode.
    pub(super) fn residue(&self, code: u8) -> Color {
        if code == b'-' || code == b'?' {
            return self.gap;
        }
        self.residues
            .get(code as usize)
            .copied()
            .flatten()
            .unwrap_or(self.fallback)
    }

    /// `color` mixed with the accent by `amount` (0 to 1), for selection and
    /// cursor washes. Indexed colors are returned unchanged.
    pub(super) fn wash(&self, color: Color, amount: f32) -> Color {
        match (color, self.accent) {
            (Color::Rgb(r, g, b), Color::Rgb(ar, ag, ab)) if self.is_truecolor => {
                let mix = |c: u8, a: u8| (c as f32 + (a as f32 - c as f32) * amount).round() as u8;
                Color::Rgb(mix(r, ar), mix(g, ag), mix(b, ab))
            }
            _ => color,
        }
    }
}

/// Whether the terminal reports 24-bit color support.
pub(super) fn is_truecolor() -> bool {
    std::env::var("COLORTERM")
        .map(|v| v.eq_ignore_ascii_case("truecolor") || v.eq_ignore_ascii_case("24bit"))
        .unwrap_or(false)
}

fn to_color(rgb: u32, is_truecolor: bool) -> Color {
    let (r, g, b) = ((rgb >> 16) as u8, (rgb >> 8) as u8, rgb as u8);
    if is_truecolor {
        Color::Rgb(r, g, b)
    } else {
        Color::Indexed(to_ansi256(r, g, b))
    }
}

// Nearest color in the xterm-256 color cube or gray ramp.
fn to_ansi256(r: u8, g: u8, b: u8) -> u8 {
    const LEVELS: [u8; 6] = [0, 95, 135, 175, 215, 255];
    let level = |c: u8| -> usize {
        LEVELS
            .iter()
            .enumerate()
            .min_by_key(|(_, l)| (**l as i32 - c as i32).abs())
            .map(|(i, _)| i)
            .unwrap_or(0)
    };
    let (ri, gi, bi) = (level(r), level(g), level(b));
    let cube = (LEVELS[ri], LEVELS[gi], LEVELS[bi]);
    let gray_index = ((r as u32 + g as u32 + b as u32) / 3).saturating_sub(8) / 10;
    let gray_index = gray_index.min(23) as u8;
    let gray = 8 + gray_index * 10;
    let distance = |(cr, cg, cb): (u8, u8, u8)| -> i32 {
        let d = |a: u8, b: u8| (a as i32 - b as i32).pow(2);
        d(cr, r) + d(cg, g) + d(cb, b)
    };
    if distance((gray, gray, gray)) < distance(cube) {
        232 + gray_index
    } else {
        16 + 36 * ri as u8 + 6 * gi as u8 + bi as u8
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_dna_palette() {
        let palette = Palette::new(Theme::Light, true, true);
        assert_eq!(palette.residue(b'A'), Color::Rgb(0x97, 0xC4, 0x59));
        assert_eq!(palette.residue(b'U'), palette.residue(b'T'));
        assert_eq!(palette.residue(b'-'), palette.gap);
        assert_eq!(palette.residue(b'N'), palette.fallback);
        assert_eq!(palette.legend.len(), 4);
    }

    #[test]
    fn test_aa_palette() {
        let palette = Palette::new(Theme::Dark, false, true);
        assert_eq!(palette.residue(b'K'), Color::Rgb(0xA3, 0x2D, 0x2D));
        assert_eq!(palette.residue(b'X'), palette.fallback);
    }

    #[test]
    fn test_ansi256() {
        assert_eq!(to_ansi256(0, 0, 0), 16);
        assert_eq!(to_ansi256(255, 255, 255), 231);
        assert_eq!(to_ansi256(128, 128, 128), 244);
    }
}
