//! WCAG 2.2 contrast helpers.
//!
//! Contrast ratio and relative luminance follow the WCAG 2.2 definitions:
//! <https://www.w3.org/TR/WCAG22/#dfn-contrast-ratio> and
//! <https://www.w3.org/TR/WCAG22/#dfn-relative-luminance>.
//! Thresholds come from success criteria 1.4.3 (text) and 1.4.11 (graphical
//! objects such as lines, markers, and bar outlines).
use super::{LUMINANCE, Rgb, dot};

/// Minimum contrast for normal text (WCAG 2.2 SC 1.4.3, level AA).
pub const AA_TEXT: f64 = 4.5;
/// Minimum contrast for large text, 18 pt or 14 pt bold (SC 1.4.3, level AA).
pub const AA_LARGE_TEXT: f64 = 3.0;
/// Minimum contrast for graphical objects needed to read a chart (SC 1.4.11).
pub const NON_TEXT: f64 = 3.0;

// Lightness step used by `ensure_contrast`.
const STEP: f64 = 0.02;

/// WCAG relative luminance in [0, 1].
pub fn relative_luminance(color: Rgb) -> f64 {
    dot(&LUMINANCE, &color.to_linear())
}

/// WCAG contrast ratio between two colors, from 1 to 21.
pub fn contrast_ratio(a: Rgb, b: Rgb) -> f64 {
    let (la, lb) = (relative_luminance(a), relative_luminance(b));
    let (light, dark) = if la > lb { (la, lb) } else { (lb, la) };
    (light + 0.05) / (dark + 0.05)
}

/// `fg` drawn with `alpha` (0 to 1) over `bg`, blended in sRGB like
/// browsers, Flutter, and plotters do.
pub fn blend(fg: Rgb, bg: Rgb, alpha: f64) -> Rgb {
    let alpha = alpha.clamp(0.0, 1.0);
    let mix = |f: u8, b: u8| (b as f64 + (f as f64 - b as f64) * alpha).round() as u8;
    Rgb::new(mix(fg.r, bg.r), mix(fg.g, bg.g), mix(fg.b, bg.b))
}

/// Black or white, whichever has more contrast with `bg`.
pub fn best_text_on(bg: Rgb) -> Rgb {
    if contrast_ratio(Rgb::BLACK, bg) >= contrast_ratio(Rgb::WHITE, bg) {
        Rgb::BLACK
    } else {
        Rgb::WHITE
    }
}

/// `color` darkened (on light backgrounds) or lightened (on dark ones) until
/// its contrast with `bg` reaches `min`.
///
/// Colors that already pass are returned unchanged. The color is mixed toward
/// black or white in small steps, which keeps its hue, so a darkened palette
/// color is still recognizable next to its fill. Returns black or white when
/// `min` cannot be reached otherwise.
pub fn ensure_contrast(color: Rgb, bg: Rgb, min: f64) -> Rgb {
    if contrast_ratio(color, bg) >= min {
        return color;
    }
    let target = best_text_on(bg);
    let mut amount = STEP;
    while amount < 1.0 {
        let adjusted = blend(target, color, amount);
        if contrast_ratio(adjusted, bg) >= min {
            return adjusted;
        }
        amount += STEP;
    }
    target
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::color::Palette;

    #[test]
    fn test_contrast_extremes() {
        assert_eq!(contrast_ratio(Rgb::BLACK, Rgb::WHITE), 21.0);
        assert_eq!(contrast_ratio(Rgb::WHITE, Rgb::WHITE), 1.0);
    }

    #[test]
    fn test_known_ratio() {
        // #767676 is the lightest gray that passes 4.5:1 on white.
        let gray = Rgb::from_u32(0x767676);
        assert!(contrast_ratio(gray, Rgb::WHITE) >= AA_TEXT);
        assert!(contrast_ratio(Rgb::from_u32(0x777777), Rgb::WHITE) < AA_TEXT);
    }

    #[test]
    fn test_best_text() {
        assert_eq!(best_text_on(Rgb::from_u32(0xF0E442)), Rgb::BLACK);
        assert_eq!(best_text_on(Rgb::from_u32(0x0072B2)), Rgb::WHITE);
    }

    #[test]
    fn test_blend() {
        let color = blend(Rgb::BLACK, Rgb::WHITE, 0.5);
        assert_eq!(color, Rgb::new(128, 128, 128));
        assert_eq!(blend(Rgb::BLACK, Rgb::WHITE, 0.0), Rgb::WHITE);
    }

    #[test]
    fn test_ensure_contrast_every_palette() {
        let backgrounds = [Rgb::WHITE, Rgb::BLACK, Rgb::from_u32(0x121212)];
        for palette in Palette::all() {
            for color in palette.colors(palette.len()) {
                for bg in backgrounds {
                    let line = ensure_contrast(color, bg, NON_TEXT);
                    assert!(
                        contrast_ratio(line, bg) >= NON_TEXT,
                        "{} {} on {}",
                        palette,
                        color.to_hex(),
                        bg.to_hex()
                    );
                }
            }
        }
    }

    #[test]
    fn test_ensure_contrast_keeps_passing_colors() {
        let blue = Rgb::from_u32(0x0072B2);
        assert_eq!(ensure_contrast(blue, Rgb::WHITE, NON_TEXT), blue);
    }
}
