//! Interactive viewers in the terminal.
pub(crate) mod align;
pub(crate) mod summary;
pub(crate) mod tree;

use ratatui::layout::Rect;

/// [value] with commas between groups of three digits.
fn thousands(value: i128) -> String {
    let digits = value.unsigned_abs().to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, ch) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(ch);
    }
    if value < 0 {
        out.insert(0, '-');
    }
    out
}

/// Counts with thousands separators.
fn count(value: usize) -> String {
    thousands(value as i128)
}

/// `text` cut to `width` characters, ending with an ellipsis when cut.
fn truncate(text: &str, width: usize) -> String {
    if text.chars().count() <= width {
        return text.to_string();
    }
    let mut cut: String = text.chars().take(width.saturating_sub(1)).collect();
    cut.push('…');
    cut
}

/// A `width` by `height` box in the middle of `area`.
fn centered(area: Rect, width: u16, height: u16) -> Rect {
    Rect::new(
        area.x + (area.width - width) / 2,
        area.y + (area.height - height) / 2,
        width,
        height,
    )
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_thousands() {
        assert_eq!(thousands(0), "0");
        assert_eq!(thousands(999), "999");
        assert_eq!(thousands(1_234_567), "1,234,567");
        assert_eq!(thousands(-1_000), "-1,000");
    }

    #[test]
    fn test_truncate() {
        assert_eq!(truncate("sample", 10), "sample");
        assert_eq!(truncate("sample_long", 6), "sampl…");
    }
}
