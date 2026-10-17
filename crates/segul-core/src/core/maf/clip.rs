//! Clip reference-oriented MAF blocks to a region.
use std::ops::Range;

/// Columns of a block whose reference positions fall in `[start, end)`.
///
/// `ref_text` is the forward-oriented reference row, and `block_start` is its
/// first reference position. A column with a reference base at position `p` is
/// kept when `start <= p < end`. A reference gap column, an insertion before the
/// next reference base `p`, is kept only when `start < p < end`, so insertions at
/// either edge of the region are dropped. The kept columns are always contiguous.
/// Returns `None` when no column is kept.
pub(crate) fn clip_columns(
    ref_text: &[u8],
    block_start: usize,
    start: usize,
    end: usize,
) -> Option<Range<usize>> {
    let mut position = block_start;
    let mut kept: Option<Range<usize>> = None;
    for (column, &base) in ref_text.iter().enumerate() {
        let is_kept = if base == b'-' {
            start < position && position < end
        } else {
            let in_region = start <= position && position < end;
            position += 1;
            in_region
        };
        if is_kept {
            match kept.as_mut() {
                Some(range) => range.end = column + 1,
                None => kept = Some(column..column + 1),
            }
        }
    }
    kept
}

#[cfg(test)]
mod tests {
    use super::*;

    // Reference positions: A10 C11 - - G12 T13 A14 - C15
    const REF: &[u8] = b"AC--GTA-C";

    fn clip(start: usize, end: usize) -> Option<String> {
        clip_columns(REF, 10, start, end)
            .map(|range| String::from_utf8(REF[range].to_vec()).unwrap())
    }

    #[test]
    fn test_clip_inside_block() {
        assert_eq!(clip(12, 15), Some("GTA".to_string()));
        assert_eq!(clip(11, 15), Some("C--GTA".to_string()));
        assert_eq!(clip(15, 16), Some("C".to_string()));
    }

    #[test]
    fn test_region_covers_block() {
        assert_eq!(clip(0, 100), Some("AC--GTA-C".to_string()));
        assert_eq!(clip(10, 16), Some("AC--GTA-C".to_string()));
    }

    #[test]
    fn test_no_overlap() {
        assert_eq!(clip(20, 30), None);
        assert_eq!(clip(0, 10), None);
    }

    #[test]
    fn test_edge_insertions_are_dropped() {
        // Insertion columns sit before position 12; a region starting at 12 drops them.
        assert_eq!(clip(12, 13), Some("G".to_string()));
        // A region ending at 12 also drops them.
        assert_eq!(clip(10, 12), Some("AC".to_string()));
    }
}
