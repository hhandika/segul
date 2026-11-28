//! Bracket comment scanner and structured metadata parsers (BEAST, FigTree, NHX).

use anyhow::{Result, bail};

use super::model::{Annotation, AnnotationFormat, AnnotationValue, Rootedness, TreeMetadata};

/// Scan a `[...]` comment starting at `bytes[start] == b'['`.
///
/// Supports nested `[...]` comments and single/double quotes inside comments,
/// returning `(inner_bytes, index_after_closing_bracket)`.
pub fn scan_bracket_comment(bytes: &[u8], start: usize) -> Result<(&[u8], usize)> {
    if bytes.get(start) != Some(&b'[') {
        bail!("Expected '[' at byte offset {}", start);
    }
    let mut depth = 1usize;
    let mut in_single = false;
    let mut in_double = false;
    let mut i = start + 1;

    while i < bytes.len() {
        let b = bytes[i];
        if b == b'\'' && !in_double {
            in_single = !in_single;
        } else if b == b'"' && !in_single {
            in_double = !in_double;
        } else if !in_single && !in_double {
            if b == b'[' {
                depth += 1;
            } else if b == b']' {
                depth -= 1;
                if depth == 0 {
                    return Ok((&bytes[start + 1..i], i + 1));
                }
            }
        }
        i += 1;
    }

    bail!("Unclosed '[' comment starting at byte offset {}", start)
}

/// Interpret tree-level `[&R]`, `[&U]`, and `[&W <weight>]` comment directives.
///
/// Returns `true` if `inner` was a pure tree-level directive and was consumed
/// into `meta`.
pub fn parse_tree_comment_directives(inner: &[u8], meta: &mut TreeMetadata) -> bool {
    let trimmed = inner.trim_ascii();
    if !trimmed.starts_with(b"&") {
        return false;
    }
    let body = trimmed[1..].trim_ascii();
    if body.eq_ignore_ascii_case(b"r") {
        meta.rootedness = Rootedness::Rooted;
        return true;
    }
    if body.eq_ignore_ascii_case(b"u") {
        meta.rootedness = Rootedness::Unrooted;
        return true;
    }
    if body.len() >= 2 && (body[0] == b'w' || body[0] == b'W') && body[1].is_ascii_whitespace() {
        let weight_str = String::from_utf8_lossy(body[1..].trim_ascii());
        if let Some(w) = parse_weight(&weight_str) {
            meta.weight = Some(w);
            return true;
        }
    }
    false
}

/// Parse a node or branch bracket comment body into a structured [`Annotation`].
pub fn parse_comment_annotation(inner: &[u8]) -> Option<Annotation> {
    let trimmed = inner.trim_ascii();
    if trimmed.is_empty() {
        return None;
    }
    let raw = String::from_utf8_lossy(trimmed).into_owned();

    if trimmed.len() >= 6 && trimmed[..6].eq_ignore_ascii_case(b"&&NHX:") {
        let entries = parse_nhx_entries(&raw[6..]);
        return Some(Annotation {
            format: AnnotationFormat::Nhx,
            entries,
            raw,
        });
    }

    if let Some(rest) = raw.strip_prefix('&') {
        let entries = parse_beast_entries(rest.trim());
        return Some(Annotation {
            format: AnnotationFormat::Beast,
            entries,
            raw,
        });
    }

    Some(Annotation {
        format: AnnotationFormat::Raw,
        entries: Vec::new(),
        raw,
    })
}

fn parse_weight(text: &str) -> Option<f64> {
    if let Some((num, den)) = text.split_once('/') {
        let n: f64 = num.trim().parse().ok()?;
        let d: f64 = den.trim().parse().ok()?;
        if d != 0.0 && n.is_finite() && d.is_finite() {
            return Some(n / d);
        }
        return None;
    }
    text.parse::<f64>().ok().filter(|v| v.is_finite())
}

fn parse_nhx_entries(body: &str) -> Vec<(String, AnnotationValue)> {
    let mut entries = Vec::new();
    for segment in body.split(':') {
        let seg = segment.trim();
        if seg.is_empty() {
            continue;
        }
        if let Some((k, v)) = seg.split_once('=') {
            entries.push((k.trim().to_string(), parse_scalar_value(v.trim())));
        } else {
            entries.push((seg.to_string(), AnnotationValue::Bool(true)));
        }
    }
    entries
}

fn parse_beast_entries(body: &str) -> Vec<(String, AnnotationValue)> {
    let mut entries = Vec::new();
    for item in split_top_level_commas(body) {
        let trimmed = item.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some((k, v)) = trimmed.split_once('=') {
            entries.push((k.trim().to_string(), parse_annotation_value(v.trim())));
        } else {
            entries.push((trimmed.to_string(), AnnotationValue::Bool(true)));
        }
    }
    entries
}

fn parse_annotation_value(raw: &str) -> AnnotationValue {
    if raw.starts_with('{') && raw.ends_with('}') && raw.len() >= 2 {
        let inner = &raw[1..raw.len() - 1];
        let items = split_top_level_commas(inner)
            .into_iter()
            .filter(|s| !s.trim().is_empty())
            .map(|s| parse_annotation_value(s.trim()))
            .collect();
        return AnnotationValue::Array(items);
    }
    if ((raw.starts_with('"') && raw.ends_with('"'))
        || (raw.starts_with('\'') && raw.ends_with('\'')))
        && raw.len() >= 2
    {
        return AnnotationValue::Str(raw[1..raw.len() - 1].to_string());
    }
    if raw.eq_ignore_ascii_case("true") {
        return AnnotationValue::Bool(true);
    }
    if raw.eq_ignore_ascii_case("false") {
        return AnnotationValue::Bool(false);
    }
    parse_scalar_value(raw)
}

fn parse_scalar_value(raw: &str) -> AnnotationValue {
    if let Ok(i) = raw.parse::<i64>() {
        return AnnotationValue::Int(i);
    }
    if let Ok(f) = raw.parse::<f64>()
        && f.is_finite()
    {
        return AnnotationValue::Float(f);
    }
    AnnotationValue::Str(raw.to_string())
}

fn split_top_level_commas(input: &str) -> Vec<&str> {
    let bytes = input.as_bytes();
    let mut parts = Vec::new();
    let mut start = 0usize;
    let mut brace_depth = 0usize;
    let mut in_single = false;
    let mut in_double = false;

    for (i, &b) in bytes.iter().enumerate() {
        if b == b'\'' && !in_double {
            in_single = !in_single;
        } else if b == b'"' && !in_single {
            in_double = !in_double;
        } else if !in_single && !in_double {
            if b == b'{' {
                brace_depth += 1;
            } else if b == b'}' && brace_depth > 0 {
                brace_depth -= 1;
            } else if b == b',' && brace_depth == 0 {
                parts.push(&input[start..i]);
                start = i + 1;
            }
        }
    }
    parts.push(&input[start..]);
    parts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scan_nested_bracket_comment() {
        let input = b"[outer [inner ']'] end]:0.1";
        let (inner, next) = scan_bracket_comment(input, 0).unwrap();
        assert_eq!(inner, b"outer [inner ']'] end");
        assert_eq!(&input[next..], b":0.1");
    }

    #[test]
    fn test_scan_unclosed_comment_errors() {
        assert!(scan_bracket_comment(b"[unclosed [nested]", 0).is_err());
    }

    #[test]
    fn test_parse_beast_annotation() {
        let annot = parse_comment_annotation(
            br##"&height_95%_HPD={0.1,0.5},posterior=0.99,!color="#ff0000",label="a, b""##,
        )
        .unwrap();
        assert_eq!(annot.format, AnnotationFormat::Beast);
        assert_eq!(annot.entries.len(), 4);
        assert_eq!(
            annot.entries[0],
            (
                "height_95%_HPD".to_string(),
                AnnotationValue::Array(vec![
                    AnnotationValue::Float(0.1),
                    AnnotationValue::Float(0.5)
                ])
            )
        );
        assert_eq!(
            annot.entries[1],
            ("posterior".to_string(), AnnotationValue::Float(0.99))
        );
        assert_eq!(
            annot.entries[2],
            (
                "!color".to_string(),
                AnnotationValue::Str("#ff0000".to_string())
            )
        );
        assert_eq!(
            annot.entries[3],
            (
                "label".to_string(),
                AnnotationValue::Str("a, b".to_string())
            )
        );
    }

    #[test]
    fn test_parse_nhx_annotation() {
        let annot = parse_comment_annotation(b"&&NHX:S=Homo_sapiens:D=N:B=95").unwrap();
        assert_eq!(annot.format, AnnotationFormat::Nhx);
        assert_eq!(
            annot.entries,
            vec![
                (
                    "S".to_string(),
                    AnnotationValue::Str("Homo_sapiens".to_string())
                ),
                ("D".to_string(), AnnotationValue::Str("N".to_string())),
                ("B".to_string(), AnnotationValue::Int(95)),
            ]
        );
    }

    #[test]
    fn test_tree_directives() {
        let mut meta = TreeMetadata::default();
        assert!(parse_tree_comment_directives(b"&R", &mut meta));
        assert_eq!(meta.rootedness, Rootedness::Rooted);
        assert!(parse_tree_comment_directives(b"&W 1/4", &mut meta));
        assert_eq!(meta.weight, Some(0.25));
        assert!(!parse_tree_comment_directives(b"&posterior=0.9", &mut meta));
    }
}
