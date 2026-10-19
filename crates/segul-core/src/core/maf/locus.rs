//! Merge clipped MAF blocks into one alignment per locus.
use indexmap::{IndexMap, IndexSet};

use crate::helper::{alphabet, types::SeqMatrix};

/// Character used for samples missing from a block.
const PAD: u8 = b'-';

/// Part of a block assigned to a locus.
pub(crate) struct LocusPiece {
    /// Forward reference interval covered by the piece.
    pub ref_start: usize,
    pub ref_end: usize,
    /// Position of the source file in the sorted input list.
    pub file_idx: usize,
    /// Position of the block in its file.
    pub block_idx: usize,
    /// Number of alignment columns.
    pub width: usize,
    pub rows: IndexMap<String, Vec<u8>>,
}

/// Statistics for a merged locus.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub(crate) struct LocusStats {
    pub blocks: usize,
    pub ref_bases: usize,
    pub columns: usize,
    pub samples: usize,
    /// Pieces dropped because they overlap an earlier piece.
    pub overlapping: usize,
}

/// Concatenate pieces in reference order and pad missing samples with gaps.
///
/// Pieces are sorted by reference position, then by file and block order, so the
/// result does not depend on input order. A piece that overlaps the previous kept
/// piece is dropped. The reference sample comes first, followed by other samples in
/// order of first appearance. Returns `None` when no sequence data remains.
pub(crate) fn build_locus(
    mut pieces: Vec<LocusPiece>,
    ref_sample: &str,
    reverse_complement: bool,
) -> (Option<SeqMatrix>, LocusStats) {
    let mut stats = LocusStats::default();
    pieces.sort_by_key(|p| (p.ref_start, p.ref_end, p.file_idx, p.block_idx));
    let mut kept: Vec<LocusPiece> = Vec::with_capacity(pieces.len());
    for piece in pieces {
        match kept.last() {
            Some(last) if piece.ref_start < last.ref_end => stats.overlapping += 1,
            _ => kept.push(piece),
        }
    }
    let mut samples: IndexSet<&str> = IndexSet::new();
    if kept.iter().any(|p| p.rows.contains_key(ref_sample)) {
        samples.insert(ref_sample);
    }
    for piece in &kept {
        samples.extend(piece.rows.keys().map(|s| s.as_str()));
    }
    stats.blocks = kept.len();
    stats.ref_bases = kept.iter().map(|p| p.ref_end - p.ref_start).sum();
    stats.columns = kept.iter().map(|p| p.width).sum();
    stats.samples = samples.len();
    if samples.is_empty() || stats.columns == 0 {
        return (None, stats);
    }
    let matrix = samples
        .iter()
        .map(|&sample| {
            let mut seq = Vec::with_capacity(stats.columns);
            for piece in &kept {
                match piece.rows.get(sample) {
                    Some(row) => seq.extend_from_slice(row),
                    None => seq.extend(std::iter::repeat_n(PAD, piece.width)),
                }
            }
            if reverse_complement {
                seq = alphabet::reverse_complement_dna(&seq);
            }
            (
                sample.to_string(),
                String::from_utf8_lossy(&seq).into_owned(),
            )
        })
        .collect();
    (Some(matrix), stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn piece(start: usize, file_idx: usize, rows: &[(&str, &str)]) -> LocusPiece {
        let width = rows.first().map(|(_, s)| s.len()).unwrap_or(0);
        LocusPiece {
            ref_start: start,
            ref_end: start + width,
            file_idx,
            block_idx: 0,
            width,
            rows: rows
                .iter()
                .map(|(k, v)| (k.to_string(), v.as_bytes().to_vec()))
                .collect(),
        }
    }

    #[test]
    fn test_padding_and_sample_order() {
        let pieces = vec![
            piece(10, 0, &[("ref", "GC"), ("c", "GG")]),
            piece(0, 0, &[("b", "AAA"), ("ref", "ACG")]),
        ];
        let (matrix, stats) = build_locus(pieces, "ref", false);
        let matrix = matrix.unwrap();
        assert_eq!(matrix.keys().collect::<Vec<_>>(), ["ref", "b", "c"]);
        assert_eq!(matrix["ref"], "ACGGC");
        assert_eq!(matrix["b"], "AAA--");
        assert_eq!(matrix["c"], "---GG");
        assert_eq!((stats.blocks, stats.columns, stats.samples), (2, 5, 3));
    }

    #[test]
    fn test_input_order_does_not_matter() {
        let a = || piece(0, 0, &[("ref", "AC"), ("x", "AA")]);
        let b = || piece(5, 1, &[("ref", "GT"), ("y", "GG")]);
        let first = build_locus(vec![a(), b()], "ref", false).0;
        let second = build_locus(vec![b(), a()], "ref", false).0;
        assert_eq!(first, second);
    }

    #[test]
    fn test_overlapping_piece_is_dropped() {
        let pieces = vec![
            piece(0, 1, &[("ref", "ACGT")]),
            piece(0, 0, &[("ref", "ACGA")]),
            piece(2, 2, &[("ref", "GT")]),
        ];
        let (matrix, stats) = build_locus(pieces, "ref", false);
        assert_eq!(matrix.unwrap()["ref"], "ACGA");
        assert_eq!(stats.overlapping, 2);
    }

    #[test]
    fn test_reverse_complement_and_empty() {
        let (matrix, _) = build_locus(vec![piece(0, 0, &[("ref", "AACg")])], "ref", true);
        assert_eq!(matrix.unwrap()["ref"], "cGTT");
        assert!(build_locus(Vec::new(), "ref", false).0.is_none());
    }
}
