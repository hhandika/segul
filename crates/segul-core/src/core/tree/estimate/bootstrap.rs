//! Parallel non-parametric bootstrapping and bipartition support mapping.

use ahash::{AHashMap, AHashSet};
use anyhow::Result;
use rayon::prelude::*;

use crate::helper::types::{DataType, SeqMatrix};
use crate::parser::tree::{NodeId, Tree};
use crate::stats::distance;

use super::distance::{estimate_bionj, estimate_nj, estimate_upgma};
use super::matrix::DistanceMatrix;
use super::parsimony::{CompressedAlignment, ParsimonyOptions, estimate_parsimony_from_compressed};
use super::rng::SeededRng;
use super::{TreeEstimateOptions, TreeMethod};

/// Run parallel non-parametric bootstrap replicates and annotate `reference` internal nodes
/// with rounded bipartition support percentages (`0..=100`).
pub(crate) fn run_bootstrap(
    reference: &Tree,
    matrix: &SeqMatrix,
    datatype: &DataType,
    options: &TreeEstimateOptions,
) -> Result<Tree> {
    let n_reps = options.bootstrap;
    if n_reps == 0 {
        return Ok(reference.clone());
    }

    let n_taxa = matrix.len();
    let mut taxa = Vec::with_capacity(n_taxa);
    let mut taxon_index = AHashMap::with_capacity(n_taxa);
    for (idx, name) in matrix.keys().enumerate() {
        taxa.push(name.clone());
        taxon_index.insert(name.as_str(), idx);
    }

    let ref_splits = extract_bipartitions(reference, &taxon_index);
    if ref_splits.is_empty() {
        return Ok(reference.clone());
    }

    let replicate_splits: Vec<AHashSet<Vec<u64>>> = match options.method {
        TreeMethod::Nj | TreeMethod::BioNj | TreeMethod::Upgma => {
            let encoded: Vec<Vec<u8>> = matrix
                .values()
                .map(|seq| distance::encode_bytes(seq, datatype))
                .collect();
            let nchar = encoded.first().map(Vec::len).unwrap_or(0);
            (0..n_reps)
                .into_par_iter()
                .map(|rep| -> Result<AHashSet<Vec<u64>>> {
                    let mut rng = SeededRng::derive(options.seed, rep as u64 + 1);
                    let mut weights = vec![0u32; nchar];
                    for _ in 0..nchar {
                        weights[rng.gen_range(nchar)] += 1;
                    }
                    let dm = DistanceMatrix::from_encoded_weighted(
                        &taxa,
                        &encoded,
                        datatype,
                        &options.model,
                        &weights,
                    )?;
                    let rep_tree = match options.method {
                        TreeMethod::Nj => estimate_nj(&dm)?,
                        TreeMethod::BioNj => estimate_bionj(&dm)?,
                        TreeMethod::Upgma => estimate_upgma(&dm)?,
                        TreeMethod::Parsimony => unreachable!(),
                    };
                    Ok(extract_split_set(&rep_tree, &taxon_index))
                })
                .collect::<Result<Vec<_>>>()?
        }
        TreeMethod::Parsimony => {
            let aln =
                CompressedAlignment::from_seq_matrix(matrix, datatype, options.parsimony.gap_mode)?;
            let rep_opts = ParsimonyOptions {
                seed: options.seed,
                ..options.parsimony.clone()
            };
            (0..n_reps)
                .into_par_iter()
                .map(|rep| -> Result<AHashSet<Vec<u64>>> {
                    let mut rng = SeededRng::derive(options.seed, rep as u64 + 1);
                    let mut pat_weights = vec![0u32; aln.n_all];
                    for _ in 0..aln.nchar {
                        let col = rng.gen_range(aln.nchar);
                        pat_weights[aln.col_to_pattern[col]] += 1;
                    }
                    let rep_aln = aln.with_resampled_weights(pat_weights);
                    let mut local_opts = rep_opts.clone();
                    local_opts.seed = options.seed.wrapping_add(rep as u64 + 1);
                    let res = estimate_parsimony_from_compressed(&rep_aln, &local_opts, false)?;
                    Ok(extract_split_set(&res.tree, &taxon_index))
                })
                .collect::<Result<Vec<_>>>()?
        }
    };

    let mut counts = vec![0usize; ref_splits.len()];
    for rep_set in &replicate_splits {
        for (i, (_, split)) in ref_splits.iter().enumerate() {
            if rep_set.contains(split) {
                counts[i] += 1;
            }
        }
    }

    let mut annotated = reference.clone();
    for (i, (node_id, _)) in ref_splits.into_iter().enumerate() {
        let pct = ((counts[i] as f64 / n_reps as f64) * 100.0).round() as u32;
        annotated.set_label(node_id, &pct.to_string(), false);
    }

    Ok(annotated)
}

/// Extract canonical non-trivial bipartitions for every non-root internal node of `tree`.
pub(crate) fn extract_bipartitions(
    tree: &Tree,
    taxon_index: &AHashMap<&str, usize>,
) -> Vec<(NodeId, Vec<u64>)> {
    let n_taxa = taxon_index.len();
    let words = n_taxa.div_ceil(64);
    let mut node_bits = vec![0u64; tree.len() * words];
    let mut splits = Vec::new();

    for id in tree.postorder() {
        let base = id.index() * words;
        if tree.is_leaf(id) {
            if let Some(lbl) = tree.label(id)
                && let Some(&tid) = taxon_index.get(lbl)
            {
                node_bits[base + (tid / 64)] |= 1u64 << (tid % 64);
            }
        } else {
            for child in tree.children(id) {
                let c_base = child.index() * words;
                for w in 0..words {
                    node_bits[base + w] |= node_bits[c_base + w];
                }
            }
            if !tree.is_root(id) {
                let mut split = node_bits[base..base + words].to_vec();
                normalize_split(&mut split, n_taxa);
                if is_nontrivial(&split, n_taxa) {
                    splits.push((id, split));
                }
            }
        }
    }

    splits
}

fn extract_split_set(tree: &Tree, taxon_index: &AHashMap<&str, usize>) -> AHashSet<Vec<u64>> {
    extract_bipartitions(tree, taxon_index)
        .into_iter()
        .map(|(_, s)| s)
        .collect()
}

fn normalize_split(bits: &mut [u64], n_taxa: usize) {
    if (bits[0] & 1) != 0 {
        for word in bits.iter_mut() {
            *word = !*word;
        }
        let rem = n_taxa % 64;
        if rem != 0 {
            let last_mask = (1u64 << rem) - 1;
            if let Some(last) = bits.last_mut() {
                *last &= last_mask;
            }
        }
    }
}

fn is_nontrivial(bits: &[u64], n_taxa: usize) -> bool {
    let count: usize = bits.iter().map(|w| w.count_ones() as usize).sum();
    count >= 2 && count <= n_taxa.saturating_sub(2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::tree::estimate::parsimony::{
        BranchSwap, GapMode, ParsimonyBranchLengths, ParsimonySearch,
    };
    use crate::core::tree::estimate::rooting::RootingMethod;
    use crate::stats::distance::DistanceModel;

    #[test]
    fn test_bootstrap_determinism_and_strong_signal() {
        let mut matrix = SeqMatrix::new();
        // Strong phylogenetic signal separating (A,B), (C,D), (E,F)
        let ab = "AAAAAAAACCCCCCCCGGGGGGGG".as_bytes();
        let cd = "CCCCCCCCAAAAAAAATTTTTTTT".as_bytes();
        let ef = "GGGGGGGGTTTTTTTTCCCCCCCC".as_bytes();
        let mut a = ab.to_vec();
        a[0] = b'T';
        let mut b = ab.to_vec();
        b[1] = b'T';
        let mut c = cd.to_vec();
        c[0] = b'G';
        let mut d = cd.to_vec();
        d[1] = b'G';
        let mut e = ef.to_vec();
        e[0] = b'A';
        let mut f = ef.to_vec();
        f[1] = b'A';

        matrix.insert("A".into(), &a);
        matrix.insert("B".into(), &b);
        matrix.insert("C".into(), &c);
        matrix.insert("D".into(), &d);
        matrix.insert("E".into(), &e);
        matrix.insert("F".into(), &f);

        let opts1 = TreeEstimateOptions {
            method: TreeMethod::Nj,
            model: DistanceModel::JukesCantor,
            bootstrap: 25,
            seed: 99,
            rooting: RootingMethod::Outgroup(vec!["E".into(), "F".into()]),
            ..Default::default()
        };
        let opts2 = opts1.clone();

        let res1 =
            crate::core::tree::estimate::estimate_matrix(&matrix, &DataType::Dna, &opts1).unwrap();
        let res2 =
            crate::core::tree::estimate::estimate_matrix(&matrix, &DataType::Dna, &opts2).unwrap();

        let nwk1 = crate::writer::tree::to_newick_string(&res1.tree, &Default::default());
        let nwk2 = crate::writer::tree::to_newick_string(&res2.tree, &Default::default());
        assert_eq!(nwk1, nwk2);

        let supports: Vec<f64> = res1
            .tree
            .preorder()
            .filter(|&id| !res1.tree.is_leaf(id) && !res1.tree.is_root(id))
            .filter_map(|id| res1.tree.support_value(id))
            .collect();
        assert!(!supports.is_empty());
        assert!(supports.iter().all(|&s| (s - 100.0).abs() < 1e-9));

        // Also test MP bootstrap determinism
        let mp_opts = TreeEstimateOptions {
            method: TreeMethod::Parsimony,
            parsimony: ParsimonyOptions {
                search: ParsimonySearch::Exact,
                swap: BranchSwap::Tbr,
                replicates: 4,
                gap_mode: GapMode::Missing,
                branch_lengths: ParsimonyBranchLengths::Acctran,
                seed: 7,
            },
            bootstrap: 10,
            seed: 7,
            ..Default::default()
        };
        let mp_res =
            crate::core::tree::estimate::estimate_matrix(&matrix, &DataType::Dna, &mp_opts)
                .unwrap();
        let mp_supports: Vec<f64> = mp_res
            .tree
            .preorder()
            .filter(|&id| !mp_res.tree.is_leaf(id) && !mp_res.tree.is_root(id))
            .filter_map(|id| mp_res.tree.support_value(id))
            .collect();
        assert!(!mp_supports.is_empty());
    }
}
