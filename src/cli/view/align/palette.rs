//! Residue colors for the terminal. The colors come from
//! `segul_core::color::residue`, shared with the SEGUI alignment viewer.

use ratatui::style::Color;
use segul_core::color::Rgb;
use segul_core::color::residue::{ResiduePalette, ResidueScheme};
use segul_core::helper::types::DataType;

pub(super) use segul_core::color::Theme;

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

/// Residue colors from `segul_core::color::residue`, as terminal colors.
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
    selection: Color,
    is_truecolor: bool,
}

impl Palette {
    pub(super) fn new(
        scheme: ResidueScheme,
        theme: Theme,
        is_dna: bool,
        is_truecolor: bool,
    ) -> Self {
        let datatype = if is_dna { DataType::Dna } else { DataType::Aa };
        let source = ResiduePalette::new(&scheme, &theme, &datatype);
        let pick = |rgb: Rgb| to_color(rgb, is_truecolor);
        let mut residues = [None; 128];
        let mut legend = Vec::with_capacity(source.groups.len());
        for group in &source.groups {
            let color = pick(group.color);
            group
                .codes
                .iter()
                .for_each(|&c| residues[c as usize] = Some(color));
            legend.push((group.label, color));
        }
        Self {
            residues,
            legend,
            fallback: pick(source.fallback),
            gap: pick(source.gap),
            on_residue: pick(source.on_residue),
            mono: pick(source.mono),
            identity: pick(source.identity),
            on_identity: pick(source.on_identity),
            accent: pick(source.accent),
            tertiary: pick(source.tertiary),
            selection: pick(source.selection),
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

    /// `color` mixed with the selection color by `amount` (0 to 0.25), for
    /// selection and cursor washes. Indexed colors are returned unchanged.
    pub(super) fn wash(&self, color: Color, amount: f32) -> Color {
        match (color, self.selection) {
            (Color::Rgb(r, g, b), Color::Rgb(ar, ag, ab)) if self.is_truecolor => {
                let mix = |c: u8, a: u8| (c as f32 + (a as f32 - c as f32) * amount).round() as u8;
                Color::Rgb(mix(r, ar), mix(g, ag), mix(b, ab))
            }
            _ => color,
        }
    }
}

/// Whether the terminal reports 24-bit color support.
pub(in crate::cli::view) fn is_truecolor() -> bool {
    std::env::var("COLORTERM")
        .map(|v| v.eq_ignore_ascii_case("truecolor") || v.eq_ignore_ascii_case("24bit"))
        .unwrap_or(false)
}

pub(in crate::cli::view) fn to_color(rgb: Rgb, is_truecolor: bool) -> Color {
    let Rgb { r, g, b } = rgb;
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
        let palette = Palette::new(ResidueScheme::Clustal, Theme::Light, true, true);
        assert_eq!(palette.residue(b'A'), Color::Rgb(0x97, 0xC4, 0x59));
        assert_eq!(palette.residue(b'U'), palette.residue(b'T'));
        assert_eq!(palette.residue(b'-'), palette.gap);
        assert_eq!(palette.residue(b'N'), palette.fallback);
        assert_eq!(palette.legend.len(), 4);
    }

    #[test]
    fn test_aa_palette() {
        let palette = Palette::new(ResidueScheme::Clustal, Theme::Dark, false, true);
        assert_eq!(palette.residue(b'K'), Color::Rgb(0xA3, 0x2D, 0x2D));
        assert_eq!(palette.residue(b'X'), palette.fallback);
    }

    #[test]
    fn test_washed_letters_meet_wcag_aa() {
        use segul_core::color::{AA_TEXT, contrast_ratio};

        let to_rgb = |c: Color| match c {
            Color::Rgb(r, g, b) => Rgb::new(r, g, b),
            _ => unreachable!(),
        };
        for scheme in ResidueScheme::all() {
            for theme in Theme::all() {
                for is_dna in [true, false] {
                    let palette = Palette::new(*scheme, *theme, is_dna, true);
                    let cells = palette.legend.iter().map(|(_, c)| *c);
                    for cell in cells.chain([palette.fallback, palette.mono]) {
                        // Selected rows wash cells 25% toward the accent.
                        let washed = to_rgb(palette.wash(cell, 0.25));
                        let ratio = contrast_ratio(to_rgb(palette.on_residue), washed);
                        assert!(ratio >= AA_TEXT, "{:?} {:?}: {:.2}", scheme, theme, ratio);
                    }
                }
            }
        }
    }

    #[test]
    fn test_ansi256() {
        assert_eq!(to_ansi256(0, 0, 0), 16);
        assert_eq!(to_ansi256(255, 255, 255), 231);
        assert_eq!(to_ansi256(128, 128, 128), 244);
    }
}
