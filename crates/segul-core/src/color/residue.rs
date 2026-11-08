//! Residue colors for alignment viewers.
//!
//! Residues are drawn as colored cells with a letter on top. In the light
//! theme cells are pale with dark letters; in the dark theme cells are deep
//! with light letters. Every letter reaches the WCAG 2.2 AA text contrast
//! (4.5:1) on every cell it can be drawn on, which the tests check.
//!
//! Amino acids use the Clustal X groups (Thompson et al. 1997,
//! <https://doi.org/10.1093/nar/25.24.4876>).
use serde::{Deserialize, Serialize};

use super::{Rgb, Theme, blend, ensure_contrast};
use crate::helper::types::DataType;

/// How residues are colored.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ResidueScheme {
    /// Clustal X groups in the soft colors of the segui app.
    #[default]
    Clustal,
    /// Okabe-Ito (DNA) and Paul Tol muted (amino acid) hues, distinct with
    /// color vision deficiencies.
    ColorblindSafe,
}

impl ResidueScheme {
    pub fn all() -> &'static [Self] {
        &[Self::Clustal, Self::ColorblindSafe]
    }

    /// Name for menus and legends.
    pub fn label(&self) -> &'static str {
        match self {
            Self::Clustal => "Clustal",
            Self::ColorblindSafe => "Colorblind safe",
        }
    }
}

impl std::fmt::Display for ResidueScheme {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Clustal => write!(f, "clustal"),
            Self::ColorblindSafe => write!(f, "colorblind-safe"),
        }
    }
}

impl std::str::FromStr for ResidueScheme {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "clustal" => Ok(Self::Clustal),
            "colorblind-safe" => Ok(Self::ColorblindSafe),
            _ => Err(format!("{} is not a valid residue scheme", s)),
        }
    }
}

// Residue groups: label, residue codes, and the Clustal colors as
// `[light, dark]`, followed by the base hue of the colorblind-safe scheme.
type Group = (&'static str, &'static [u8], [u32; 2], u32);

const DNA_GROUPS: &[Group] = &[
    ("A", b"A", [0x97C459, 0x3B6D11], 0x009E73),
    ("C", b"C", [0x85B7EB, 0x185FA5], 0x0072B2),
    ("G", b"G", [0xFAC775, 0x854F0B], 0xF0E442),
    ("T", b"TU", [0xF09595, 0xA32D2D], 0xD55E00),
];

const AA_GROUPS: &[Group] = &[
    ("Hydrophobic", b"AILMFWV", [0x85B7EB, 0x185FA5], 0x332288),
    ("Positive", b"KR", [0xF09595, 0xA32D2D], 0xCC6677),
    ("Negative", b"DE", [0xED93B1, 0x993556], 0x882255),
    ("Polar", b"NQST", [0x97C459, 0x3B6D11], 0x117733),
    ("Aromatic", b"HY", [0x5DCAA5, 0x0F6E56], 0x44AA99),
    ("G", b"G", [0xFAC775, 0x854F0B], 0xDDCC77),
    ("P", b"P", [0xFAEEDA, 0x633806], 0x999933),
    ("C", b"C", [0xF5C4B3, 0x993C1D], 0xAA4499),
];

// Light and dark values of the colors shared by both schemes.
const FALLBACK: [u32; 2] = [0xB4B2A9, 0x5F5E5A];
const GAP: [u32; 2] = [0xE6E4DC, 0x2C2C2A];
const ON_RESIDUE: [u32; 2] = [0x2C2C2A, 0xF1EFE8];
const MONO: [u32; 2] = [0xD3D1C7, 0x444441];
const IDENTITY: [u32; 2] = [0xF1EFE8, 0x242422];
// Darker (light theme) and lighter (dark theme) than the earlier
// #888780 and #8A8984, which were 3.1:1 and 4.4:1, so they also pass on
// selected cells.
const ON_IDENTITY: [u32; 2] = [0x5F5E58, 0x9A9994];
const ACCENT: [u32; 2] = [0x0F6E56, 0x5DCAA5];
const TERTIARY: [u32; 2] = [0x993556, 0xED93B1];
// Selection washes, lighter than pale cells' letters and deeper than dark
// cells' letters. Washing toward the accent instead dropped letters to
// 3.3:1 (light) and 3.9:1 (dark).
const SELECTION: [u32; 2] = [0x5DCAA5, 0x085041];
// Largest share of the selection color a viewer mixes into a cell.
const MAX_WASH: f64 = 0.25;

// Letter contrast the colorblind-safe cells are tinted to. Above AA so
// selection washes in viewers keep letters readable.
const TINT_CONTRAST: f64 = 6.0;

/// A residue group, such as "Hydrophobic", and its cell color.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ResidueGroup {
    pub label: &'static str,
    pub codes: &'static [u8],
    pub color: Rgb,
}

/// Every color an alignment viewer needs for one scheme, theme, and data type.
#[derive(Debug, Clone, PartialEq)]
pub struct ResiduePalette {
    /// Residue groups in legend order.
    pub groups: Vec<ResidueGroup>,
    /// Cell of residues outside every group, such as `N` or `X`.
    pub fallback: Rgb,
    /// Cell of gaps (`-`) and missing data (`?`).
    pub gap: Rgb,
    /// Letter on residue, fallback, gap, and mono cells.
    pub on_residue: Rgb,
    /// Cell of every residue in single-color mode.
    pub mono: Rgb,
    /// Cell of a residue equal to the first sequence.
    pub identity: Rgb,
    /// Letter on identity cells.
    pub on_identity: Rgb,
    /// Selection, cursor, and highlight color.
    pub accent: Rgb,
    /// Secondary highlight, such as chunk marks.
    pub tertiary: Rgb,
    /// Color mixed into selected cells, at most [`ResiduePalette::MAX_WASH`].
    pub selection: Rgb,
    lookup: [Option<Rgb>; 128],
}

impl ResiduePalette {
    /// Largest share of [`selection`](Self::selection) to mix into a cell
    /// while letters keep 4.5:1.
    pub const MAX_WASH: f64 = MAX_WASH;

    pub fn new(scheme: &ResidueScheme, theme: &Theme, datatype: &DataType) -> Self {
        let on_residue = theme.pick(ON_RESIDUE);
        let groups = match datatype {
            DataType::Aa => AA_GROUPS,
            _ => DNA_GROUPS,
        };
        let groups: Vec<ResidueGroup> = groups
            .iter()
            .map(|(label, codes, clustal, base)| {
                let color = match scheme {
                    ResidueScheme::Clustal => theme.pick(*clustal),
                    ResidueScheme::ColorblindSafe => {
                        ensure_contrast(Rgb::from_u32(*base), on_residue, TINT_CONTRAST)
                    }
                };
                ResidueGroup {
                    label,
                    codes,
                    color,
                }
            })
            .collect();
        let mut lookup = [None; 128];
        groups.iter().for_each(|g| {
            g.codes
                .iter()
                .for_each(|&c| lookup[c as usize] = Some(g.color))
        });
        Self {
            groups,
            fallback: theme.pick(FALLBACK),
            gap: theme.pick(GAP),
            on_residue,
            mono: theme.pick(MONO),
            identity: theme.pick(IDENTITY),
            on_identity: theme.pick(ON_IDENTITY),
            accent: theme.pick(ACCENT),
            tertiary: theme.pick(TERTIARY),
            selection: theme.pick(SELECTION),
            lookup,
        }
    }

    /// Cell color of a residue code in base mode.
    pub fn residue(&self, code: u8) -> Rgb {
        if code == b'-' || code == b'?' {
            return self.gap;
        }
        self.lookup
            .get(code as usize)
            .copied()
            .flatten()
            .unwrap_or(self.fallback)
    }

    /// Every (letter, cell) pair the palette can draw, including selected
    /// cells, for contrast checks.
    pub fn text_pairs(&self) -> Vec<(Rgb, Rgb)> {
        let mut pairs: Vec<(Rgb, Rgb)> = self
            .groups
            .iter()
            .map(|g| (self.on_residue, g.color))
            .collect();
        pairs.extend([
            (self.on_residue, self.fallback),
            (self.on_residue, self.gap),
            (self.on_residue, self.mono),
            (self.on_identity, self.identity),
        ]);
        let selected: Vec<(Rgb, Rgb)> = pairs
            .iter()
            .map(|(text, cell)| (*text, blend(self.selection, *cell, MAX_WASH)))
            .collect();
        pairs.extend(selected);
        pairs
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::color::{AA_TEXT, contrast_ratio};

    fn every_palette() -> Vec<ResiduePalette> {
        let mut palettes = Vec::new();
        for scheme in ResidueScheme::all() {
            for theme in Theme::all() {
                for datatype in [DataType::Dna, DataType::Aa] {
                    palettes.push(ResiduePalette::new(scheme, theme, &datatype));
                }
            }
        }
        palettes
    }

    #[test]
    fn test_dna_palette() {
        let palette = ResiduePalette::new(&ResidueScheme::Clustal, &Theme::Light, &DataType::Dna);
        assert_eq!(palette.residue(b'A'), Rgb::from_u32(0x97C459));
        assert_eq!(palette.residue(b'U'), palette.residue(b'T'));
        assert_eq!(palette.residue(b'-'), palette.gap);
        assert_eq!(palette.residue(b'N'), palette.fallback);
        assert_eq!(palette.groups.len(), 4);
    }

    #[test]
    fn test_aa_palette() {
        let palette = ResiduePalette::new(&ResidueScheme::Clustal, &Theme::Dark, &DataType::Aa);
        assert_eq!(palette.residue(b'K'), Rgb::from_u32(0xA32D2D));
        assert_eq!(palette.residue(b'X'), palette.fallback);
        assert_eq!(palette.groups.len(), 8);
    }

    #[test]
    fn test_letters_meet_wcag_aa() {
        for palette in every_palette() {
            for (text, cell) in palette.text_pairs() {
                let ratio = contrast_ratio(text, cell);
                assert!(
                    ratio >= AA_TEXT,
                    "{} on {} is {:.2}",
                    text.to_hex(),
                    cell.to_hex(),
                    ratio
                );
            }
        }
    }

    #[test]
    fn test_accents_visible_on_terminal_background() {
        let light = ResiduePalette::new(&ResidueScheme::Clustal, &Theme::Light, &DataType::Dna);
        let dark = ResiduePalette::new(&ResidueScheme::Clustal, &Theme::Dark, &DataType::Dna);
        for color in [light.accent, light.tertiary] {
            assert!(contrast_ratio(color, Rgb::WHITE) >= AA_TEXT);
        }
        for color in [dark.accent, dark.tertiary] {
            assert!(contrast_ratio(color, Rgb::BLACK) >= AA_TEXT);
        }
    }

    #[test]
    fn test_colorblind_groups_differ() {
        for palette in every_palette() {
            let colors: Vec<Rgb> = palette.groups.iter().map(|g| g.color).collect();
            for (i, a) in colors.iter().enumerate() {
                for b in &colors[i + 1..] {
                    assert_ne!(a, b);
                }
            }
        }
    }

    #[test]
    fn test_scheme_names() {
        for scheme in ResidueScheme::all() {
            assert_eq!(
                scheme.to_string().parse::<ResidueScheme>().unwrap(),
                *scheme
            );
            let json = serde_json::to_string(scheme).unwrap();
            assert_eq!(json, format!("\"{}\"", scheme));
        }
    }
}
