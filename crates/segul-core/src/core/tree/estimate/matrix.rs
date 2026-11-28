//! Symmetric pairwise distance matrix construction for tree estimation.

use anyhow::{Result, bail};
use rayon::prelude::*;

use crate::helper::types::{DataType, SeqMatrix};
use crate::stats::distance::{self, BaseFreqs, DistanceModel, INVALID, PairCounts};

/// Symmetric pairwise distance matrix ready for distance-based tree estimation.
#[derive(Debug, Clone)]
pub struct DistanceMatrix {
    /// Ordered taxon names matching rows and columns.
    pub taxa: Vec<String>,
    /// Number of taxa (`>= 3`).
    pub n: usize,
    /// Row-major `n x n` symmetric distances with `0.0` on the diagonal.
    pub values: Vec<f64>,
    /// Count of off-diagonal pairs that were undefined or saturated and imputed.
    pub undefined_pairs: usize,
}

impl DistanceMatrix {
    /// Compute a symmetric distance matrix from an aligned [`SeqMatrix`].
    pub fn from_seq_matrix(
        matrix: &SeqMatrix,
        datatype: &DataType,
        model: &DistanceModel,
    ) -> Result<Self> {
        let n = matrix.len();
        if n < 3 {
            bail!(
                "At least 3 taxa are required for tree estimation (found {})",
                n
            );
        }
        if !model.is_valid_for(datatype) {
            bail!(
                "The {} model is not valid for {} data. \
                Use p, jc69, k2p, or tn93 for DNA, and p, jc69, or poisson for amino acid.",
                model,
                datatype
            );
        }
        let mut seq_len = None;
        for (id, seq) in matrix.iter() {
            if seq.is_empty() {
                bail!("Sequence '{}' is empty", id);
            }
            match seq_len {
                None => seq_len = Some(seq.len()),
                Some(len) if seq.len() != len => {
                    bail!(
                        "Sequences must be aligned with equal length (expected {}, '{}' has {})",
                        len,
                        id,
                        seq.len()
                    );
                }
                _ => {}
            }
        }

        let taxa: Vec<String> = matrix.keys().cloned().collect();
        let encoded: Vec<Vec<u8>> = matrix
            .values()
            .map(|seq| distance::encode_bytes(seq, datatype))
            .collect();

        Self::from_encoded_inner(&taxa, &encoded, datatype, model, None)
    }

    /// Compute a distance matrix from pre-encoded sequences using per-site integer weights.
    pub(crate) fn from_encoded_weighted(
        taxa: &[String],
        encoded: &[Vec<u8>],
        datatype: &DataType,
        model: &DistanceModel,
        weights: &[u32],
    ) -> Result<Self> {
        Self::from_encoded_inner(taxa, encoded, datatype, model, Some(weights))
    }

    /// Return the pairwise distance between taxa `i` and `j`.
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.values[i * self.n + j]
    }

    /// Set the distance between taxa `i` and `j` (single cell).
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, value: f64) {
        self.values[i * self.n + j] = value;
    }

    fn from_encoded_inner(
        taxa: &[String],
        encoded: &[Vec<u8>],
        datatype: &DataType,
        model: &DistanceModel,
        weights: Option<&[u32]>,
    ) -> Result<Self> {
        let n = taxa.len();
        if n < 3 {
            bail!(
                "At least 3 taxa are required for tree estimation (found {})",
                n
            );
        }
        if !model.is_valid_for(datatype) {
            bail!("The {} model is not valid for {} data", model, datatype);
        }

        let freqs = match model {
            DistanceModel::TamuraNei => match weights {
                Some(w) => weighted_base_freqs(encoded, w),
                None => BaseFreqs::from_encoded(encoded),
            },
            _ => BaseFreqs::default(),
        };

        let compute_row = |i: usize| -> Vec<Option<f64>> {
            (i + 1..n)
                .map(|j| {
                    let counts = match weights {
                        Some(w) => weighted_pair_counts(&encoded[i], &encoded[j], datatype, w),
                        None => PairCounts::new(&encoded[i], &encoded[j], datatype),
                    };
                    model.distance(&counts, datatype, &freqs)
                })
                .collect()
        };

        let upper_rows: Vec<Vec<Option<f64>>> = if weights.is_none() && n >= 16 {
            (0..n).into_par_iter().map(compute_row).collect()
        } else {
            (0..n).map(compute_row).collect()
        };

        let mut max_finite = 0.0f64;
        let mut has_finite = false;
        let mut undefined_pairs = 0usize;

        for row in &upper_rows {
            for cell in row {
                match cell {
                    Some(d) if d.is_finite() => {
                        has_finite = true;
                        if *d > max_finite {
                            max_finite = *d;
                        }
                    }
                    _ => undefined_pairs += 1,
                }
            }
        }

        if !has_finite {
            bail!(
                "No finite pairwise distances could be computed under the {} model",
                model
            );
        }

        let imputed = if max_finite > 0.0 {
            max_finite * 1.5
        } else {
            1.0
        };

        let mut values = vec![0.0f64; n * n];
        for (i, row) in upper_rows.into_iter().enumerate() {
            for (offset, cell) in row.into_iter().enumerate() {
                let j = i + 1 + offset;
                let d = match cell {
                    Some(v) if v.is_finite() => v,
                    _ => imputed,
                };
                values[i * n + j] = d;
                values[j * n + i] = d;
            }
        }

        Ok(Self {
            taxa: taxa.to_vec(),
            n,
            values,
            undefined_pairs,
        })
    }
}

fn weighted_base_freqs(encoded: &[Vec<u8>], weights: &[u32]) -> BaseFreqs {
    let mut counts = [0u64; 4];
    for seq in encoded {
        for (&b, &w) in seq.iter().zip(weights) {
            if w > 0 && b < 4 {
                counts[b as usize] += u64::from(w);
            }
        }
    }
    let total: u64 = counts.iter().sum();
    if total == 0 {
        return BaseFreqs::default();
    }
    let total = total as f64;
    BaseFreqs {
        a: counts[0] as f64 / total,
        c: counts[1] as f64 / total,
        g: counts[2] as f64 / total,
        t: counts[3] as f64 / total,
    }
}

fn weighted_pair_counts(
    seq_1: &[u8],
    seq_2: &[u8],
    datatype: &DataType,
    weights: &[u32],
) -> PairCounts {
    let mut counts = PairCounts::default();
    let is_dna = !matches!(datatype, DataType::Aa);
    for ((&a, &b), &w) in seq_1.iter().zip(seq_2).zip(weights) {
        if w == 0 || a == INVALID || b == INVALID {
            continue;
        }
        let w = w as usize;
        counts.sites += w;
        if a == b {
            continue;
        }
        counts.diffs += w;
        if is_dna {
            match (a.min(b), a.max(b)) {
                (0, 2) => counts.ts_purine += w,
                (1, 3) => counts.ts_pyrimidine += w,
                _ => counts.transversions += w,
            }
        }
    }
    counts
}

#[cfg(test)]
mod tests {
    use indexmap::IndexSet;

    use super::*;

    #[test]
    fn test_distance_matrix_rejects_fewer_than_three_taxa() {
        let mut matrix = SeqMatrix::new();
        matrix.insert("A".to_string(), b"ACGT");
        matrix.insert("B".to_string(), b"ACGA");
        let err =
            DistanceMatrix::from_seq_matrix(&matrix, &DataType::Dna, &DistanceModel::Uncorrected)
                .unwrap_err();
        assert!(err.to_string().contains("At least 3 taxa"));
    }

    #[test]
    fn test_distance_matrix_imputes_isolated_undefined_pair() {
        let mut ids = IndexSet::new();
        ids.insert("A".to_string());
        ids.insert("B".to_string());
        ids.insert("C".to_string());
        // A and C have no overlapping valid sites; A-B and B-C are defined.
        let mut matrix = SeqMatrix::new();
        matrix.insert("A".to_string(), b"AA--");
        matrix.insert("B".to_string(), b"ATAA");
        matrix.insert("C".to_string(), b"--AA");

        let dm =
            DistanceMatrix::from_seq_matrix(&matrix, &DataType::Dna, &DistanceModel::Uncorrected)
                .unwrap();
        assert_eq!(dm.undefined_pairs, 1);
        assert!((dm.get(0, 1) - 0.5).abs() < 1e-12);
        assert!((dm.get(1, 2) - 0.0).abs() < 1e-12);
        assert!((dm.get(0, 2) - 0.75).abs() < 1e-12);
    }

    #[test]
    fn test_weighted_matches_unweighted() {
        let mut matrix = SeqMatrix::new();
        matrix.insert("A".to_string(), b"AACGT");
        matrix.insert("B".to_string(), b"GATGT");
        matrix.insert("C".to_string(), b"AACGA");

        let dm1 =
            DistanceMatrix::from_seq_matrix(&matrix, &DataType::Dna, &DistanceModel::Kimura2P)
                .unwrap();
        let taxa: Vec<String> = matrix.keys().cloned().collect();
        let encoded: Vec<Vec<u8>> = matrix
            .values()
            .map(|s| distance::encode_bytes(s, &DataType::Dna))
            .collect();
        let weights = vec![1u32; 5];
        let dm2 = DistanceMatrix::from_encoded_weighted(
            &taxa,
            &encoded,
            &DataType::Dna,
            &DistanceModel::Kimura2P,
            &weights,
        )
        .unwrap();
        assert_eq!(dm1.values, dm2.values);
    }
}
