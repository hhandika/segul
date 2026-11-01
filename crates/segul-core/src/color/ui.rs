//! Chart and table colors other than the data series: background, axes,
//! grid, text, and status colors.
//!
//! Text reaches 4.5:1 and axes reach 3:1 against the background (WCAG 2.2
//! SC 1.4.3 and 1.4.11). Grid lines are decorative and lighter on purpose:
//! SC 1.4.11 does not cover them because the axes carry the scale.
use serde::{Deserialize, Serialize};

use super::{NON_TEXT, Rgb, Theme, ensure_contrast};

/// Pass, warn, and fail colors for quality summaries, plus highlights.
///
/// Taken from Okabe-Ito so they stay distinct with color vision
/// deficiencies. Use [`ChartTokens::line`] before drawing them as lines or
/// dots, because orange and gray fall below 3:1 on white.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct StatusColors {
    pub pass: Rgb,
    pub warn: Rgb,
    pub fail: Rgb,
    /// The sample or item the user opened.
    pub open: Rgb,
    /// Outliers and other de-emphasized marks.
    pub outlier: Rgb,
}

impl Default for StatusColors {
    fn default() -> Self {
        Self {
            pass: Rgb::from_u32(0x009E73),
            warn: Rgb::from_u32(0xE69F00),
            fail: Rgb::from_u32(0xD55E00),
            open: Rgb::from_u32(0xCC79A7),
            outlier: Rgb::from_u32(0x999999),
        }
    }
}

/// Colors for chart chrome in one theme.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct ChartTokens {
    pub background: Rgb,
    /// Titles, labels, and tick text.
    pub text: Rgb,
    /// Axis lines, box plot whiskers, and outlines.
    pub axis: Rgb,
    /// Grid lines. Decorative.
    pub grid: Rgb,
    /// Heatmap cells with no value. Drawn with an `axis` outline.
    pub missing: Rgb,
    /// Legend box behind series names.
    pub legend_background: Rgb,
    pub status: StatusColors,
}

impl ChartTokens {
    pub fn new(theme: &Theme) -> Self {
        let pick = |pair: [u32; 2]| theme.pick(pair);
        Self {
            background: pick([0xFFFFFF, 0x121212]),
            text: pick([0x000000, 0xE8E8E8]),
            axis: pick([0x5A5A5A, 0xA0A0A0]),
            grid: pick([0xE1E1E1, 0x333333]),
            missing: pick([0xDCDCDC, 0x3A3A3A]),
            legend_background: pick([0xFFFFFF, 0x1E1E1E]),
            status: StatusColors::default(),
        }
    }

    /// `color` adjusted to reach 3:1 against the background, for lines,
    /// points, and outlines. Fills can keep the original color.
    pub fn line(&self, color: Rgb) -> Rgb {
        ensure_contrast(color, self.background, NON_TEXT)
    }
}

impl Default for ChartTokens {
    fn default() -> Self {
        Self::new(&Theme::Light)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::color::{AA_TEXT, contrast_ratio};

    #[test]
    fn test_tokens_meet_wcag() {
        for theme in Theme::all() {
            let tokens = ChartTokens::new(theme);
            assert!(contrast_ratio(tokens.text, tokens.background) >= AA_TEXT);
            assert!(contrast_ratio(tokens.text, tokens.legend_background) >= AA_TEXT);
            assert!(contrast_ratio(tokens.axis, tokens.background) >= NON_TEXT);
            assert!(contrast_ratio(tokens.axis, tokens.missing) >= NON_TEXT);
        }
    }

    #[test]
    fn test_status_lines_meet_wcag() {
        for theme in Theme::all() {
            let tokens = ChartTokens::new(theme);
            let s = tokens.status;
            for color in [s.pass, s.warn, s.fail, s.open, s.outlier] {
                let line = tokens.line(color);
                assert!(contrast_ratio(line, tokens.background) >= NON_TEXT);
            }
        }
    }
}
