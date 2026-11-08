//! Interactive viewers in the terminal.
pub(crate) mod align;
pub(crate) mod summary;

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
}
