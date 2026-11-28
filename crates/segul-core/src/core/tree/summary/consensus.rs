//! Consensus tree construction (Strict, MRC, MRE) and MCC/MAP tree selection.

use std::cmp::Ordering;
use std::sync::Arc;

use ahash::AHashSet;
use anyhow::{Result, bail};

use crate::parser::tree::{NodeId, Rootedness, Tree};

use super::bitset::{
    CladeId, are_compatible_rooted, are_compatible_unrooted, cmp_bits, first_set_bit, is_subset,
    popcount,
};
use super::engine::{GlobalCladeTable, TopologyRecord, TreeCladeRecord};
use super::{MccScoring, SummaryMethod};

/// Identifies where a retained input tree lives so MCC/MAP can reload only the single winner.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum TreeLocation {
    InMemory(usize),
    File { file_idx: usize, tree_idx: usize },
}

/// Build an iterative Strict, Majority-Rule (MRC), or Extended Majority-Rule (MRE) tree.
///
/// Returns `(tree, node_clade_ids)` where `node_clade_ids[node.index()]` maps each internal
/// non-root node to its [`CladeId`] in `table`.
pub(crate) fn build_consensus_tree(
    table: &GlobalCladeTable,
    total_trees: u64,
    method: SummaryMethod,
    min_clade_prob: f64,
    rooted: bool,
) -> Result<(Tree, Vec<Option<CladeId>>)> {
    if total_trees == 0 {
        bail!("Cannot build consensus tree from 0 trees");
    }
    let n_taxa = table.taxa.len();
    if n_taxa == 0 {
        bail!("Cannot build consensus tree with 0 taxa");
    }

    let mut accepted: Vec<CladeId> = Vec::new();
    match method {
        SummaryMethod::Strict => {
            for (idx, stats) in table.clades.iter().enumerate() {
                if stats.count == total_trees {
                    accepted.push(CladeId(idx as u32));
                }
            }
        }
        SummaryMethod::Mrc => {
            if !(0.5..=1.0).contains(&min_clade_prob) {
                bail!(
                    "Majority-rule consensus threshold must be in [0.5, 1.0], got {}",
                    min_clade_prob
                );
            }
            let strict_half = (min_clade_prob - 0.5).abs() < 1e-12;
            for (idx, stats) in table.clades.iter().enumerate() {
                let keep = if strict_half {
                    stats.count * 2 > total_trees
                } else {
                    (stats.count as f64 / total_trees as f64) >= min_clade_prob
                };
                if keep {
                    accepted.push(CladeId(idx as u32));
                }
            }
        }
        SummaryMethod::Mre => {
            if !(0.0..=1.0).contains(&min_clade_prob) {
                bail!(
                    "Minimum clade probability must be in [0.0, 1.0], got {}",
                    min_clade_prob
                );
            }
            let mut candidates: Vec<(CladeId, u64)> = Vec::new();
            for (idx, stats) in table.clades.iter().enumerate() {
                let cid = CladeId(idx as u32);
                if stats.count * 2 > total_trees {
                    accepted.push(cid);
                } else if (stats.count as f64 / total_trees as f64) >= min_clade_prob {
                    candidates.push((cid, stats.count));
                }
            }
            candidates.sort_by(|&(a_id, a_cnt), &(b_id, b_cnt)| {
                b_cnt
                    .cmp(&a_cnt)
                    .then_with(|| cmp_bits(table.arena.get(a_id), table.arena.get(b_id)))
            });

            let max_internal = if rooted {
                n_taxa.saturating_sub(2)
            } else {
                n_taxa.saturating_sub(3)
            };
            for (cid, _) in candidates {
                if accepted.len() >= max_internal {
                    break;
                }
                let cand_bits = table.arena.get(cid);
                let compat = accepted.iter().all(|&acc_id| {
                    let acc_bits = table.arena.get(acc_id);
                    if rooted {
                        are_compatible_rooted(cand_bits, acc_bits)
                    } else {
                        are_compatible_unrooted(cand_bits, acc_bits, n_taxa)
                    }
                });
                if compat {
                    accepted.push(cid);
                }
            }
        }
        _ => bail!("build_consensus_tree called with non-consensus method"),
    }

    build_tree_from_compatible_clades(table, &accepted, rooted)
}

/// Construct a [`Tree`] iteratively from a pairwise-compatible set of clades.
fn build_tree_from_compatible_clades(
    table: &GlobalCladeTable,
    accepted: &[CladeId],
    rooted: bool,
) -> Result<(Tree, Vec<Option<CladeId>>)> {
    let n_taxa = table.taxa.len();

    // Order clades from largest to smallest so every clade attaches under its tightest superset.
    let mut ordered: Vec<(CladeId, usize, usize)> = accepted
        .iter()
        .map(|&cid| {
            let bits = table.arena.get(cid);
            (
                cid,
                popcount(bits),
                first_set_bit(bits).unwrap_or(usize::MAX),
            )
        })
        .collect();
    ordered.sort_by(|a, b| {
        b.1.cmp(&a.1)
            .then_with(|| a.2.cmp(&b.2))
            .then_with(|| cmp_bits(table.arena.get(a.0), table.arena.get(b.0)))
    });

    let mut tree = Tree::with_capacity(n_taxa + ordered.len() + 1, 0);
    tree.set_taxa(Arc::clone(&table.taxa));
    tree.metadata.rootedness = if rooted {
        Rootedness::Rooted
    } else {
        Rootedness::Unrooted
    };
    let root = tree.add_node();
    tree.set_root(root);

    // Track internal children of each node along with their clade bitsets for iterative lookup.
    struct NodeEntry {
        node: NodeId,
        clade: Option<CladeId>,
        min_taxon: usize,
        internal_children: Vec<usize>,
    }

    let mut nodes: Vec<NodeEntry> = Vec::with_capacity(ordered.len() + 1);
    nodes.push(NodeEntry {
        node: root,
        clade: None,
        min_taxon: 0,
        internal_children: Vec::new(),
    });

    for (cid, _, first_bit) in ordered {
        let bits = table.arena.get(cid);
        let mut cur_entry = 0usize;
        loop {
            let mut next_entry = None;
            for &child_entry_idx in &nodes[cur_entry].internal_children {
                let child_cid = nodes[child_entry_idx]
                    .clade
                    .expect("Internal entry always has CladeId");
                if is_subset(bits, table.arena.get(child_cid)) {
                    next_entry = Some(child_entry_idx);
                    break;
                }
            }
            if let Some(next) = next_entry {
                cur_entry = next;
            } else {
                break;
            }
        }
        let new_node = tree.add_node();
        let new_entry_idx = nodes.len();
        nodes.push(NodeEntry {
            node: new_node,
            clade: Some(cid),
            min_taxon: first_bit,
            internal_children: Vec::new(),
        });
        nodes[cur_entry].internal_children.push(new_entry_idx);
    }

    // Determine the deepest containing internal node for each tip.
    let mut tip_assignments: Vec<Vec<(usize, NodeId)>> = vec![Vec::new(); nodes.len()];
    for tid in 0..n_taxa {
        let mut cur_entry = 0usize;
        if rooted || tid != 0 {
            loop {
                let mut next_entry = None;
                for &child_entry_idx in &nodes[cur_entry].internal_children {
                    let child_cid = nodes[child_entry_idx]
                        .clade
                        .expect("Internal entry always has CladeId");
                    let bits = table.arena.get(child_cid);
                    if (bits[tid / 64] & (1u64 << (tid % 64))) != 0 {
                        next_entry = Some(child_entry_idx);
                        break;
                    }
                }
                if let Some(next) = next_entry {
                    cur_entry = next;
                } else {
                    break;
                }
            }
        }
        let leaf = tree.add_node();
        let taxon_id = crate::parser::tree::TaxonId::from_raw(tid as u32);
        tree.set_interned_label(leaf, taxon_id, false);
        tip_assignments[cur_entry].push((tid, leaf));
    }

    // Attach children in deterministic smallest-taxon order.
    enum ChildKind {
        Tip(usize, NodeId),
        Internal(usize, NodeId),
    }
    impl ChildKind {
        fn min_taxon(&self) -> usize {
            match *self {
                Self::Tip(t, _) | Self::Internal(t, _) => t,
            }
        }
        fn node_id(&self) -> NodeId {
            match *self {
                Self::Tip(_, n) | Self::Internal(_, n) => n,
            }
        }
    }

    let mut node_clade_ids = vec![None; tree.len()];
    for (entry_idx, entry) in nodes.iter().enumerate() {
        if let Some(cid) = entry.clade {
            node_clade_ids[entry.node.index()] = Some(cid);
        }
        let mut children: Vec<ChildKind> =
            Vec::with_capacity(tip_assignments[entry_idx].len() + entry.internal_children.len());
        for &(tid, leaf) in &tip_assignments[entry_idx] {
            children.push(ChildKind::Tip(tid, leaf));
        }
        for &child_entry_idx in &entry.internal_children {
            children.push(ChildKind::Internal(
                nodes[child_entry_idx].min_taxon,
                nodes[child_entry_idx].node,
            ));
        }
        children.sort_by_key(ChildKind::min_taxon);
        let mut prev = None;
        for ch in children {
            let cid_node = ch.node_id();
            tree.attach_child(entry.node, cid_node, prev);
            prev = Some(cid_node);
        }
    }

    Ok((tree, node_clade_ids))
}

/// Select the Maximum Clade Credibility (MCC) tree location and score.
pub(crate) fn select_mcc_tree(
    records: &[TreeCladeRecord],
    table: &GlobalCladeTable,
    topologies: &[TopologyRecord],
    total_trees: u64,
    scoring: MccScoring,
) -> Result<(TreeLocation, f64)> {
    if records.is_empty() || total_trees == 0 {
        bail!("Cannot select MCC tree from 0 trees");
    }

    let mut topo_counts = ahash::AHashMap::with_capacity(topologies.len());
    for topo in topologies {
        topo_counts.insert(topo.hash, topo.count);
    }

    let total_f = total_trees as f64;
    let mut best_loc = records[0].location;
    let mut best_score = f64::NEG_INFINITY;
    let mut best_topo_count = 0u64;

    for rec in records {
        let score = score_clade_list(&rec.clade_ids, table, total_f, scoring);
        let topo_cnt = topo_counts
            .get(&hash_clade_list(&rec.clade_ids))
            .copied()
            .unwrap_or(1);

        let better = match cmp_score(score, best_score) {
            Ordering::Greater => true,
            Ordering::Less => false,
            Ordering::Equal => {
                topo_cnt > best_topo_count
                    || (topo_cnt == best_topo_count && rec.location < best_loc)
            }
        };
        if better {
            best_score = score;
            best_loc = rec.location;
            best_topo_count = topo_cnt;
        }
    }

    Ok((best_loc, best_score))
}

/// Select the Maximum A Posteriori (MAP) topology tree location, frequency, and MCC score.
pub(crate) fn select_map_tree(
    topologies: &[TopologyRecord],
    table: &GlobalCladeTable,
    total_trees: u64,
    scoring: MccScoring,
) -> Result<(TreeLocation, f64, f64)> {
    if topologies.is_empty() || total_trees == 0 {
        bail!("Cannot select MAP tree from 0 trees");
    }
    let total_f = total_trees as f64;
    let mut best_topo = &topologies[0];
    let mut best_score = score_clade_list(&best_topo.clade_ids, table, total_f, scoring);

    for topo in topologies.iter().skip(1) {
        let score = score_clade_list(&topo.clade_ids, table, total_f, scoring);
        let better = match topo.count.cmp(&best_topo.count) {
            Ordering::Greater => true,
            Ordering::Less => false,
            Ordering::Equal => match cmp_score(score, best_score) {
                Ordering::Greater => true,
                Ordering::Less => false,
                Ordering::Equal => topo.first_location < best_topo.first_location,
            },
        };
        if better {
            best_topo = topo;
            best_score = score;
        }
    }

    let freq = best_topo.count as f64 / total_f;
    Ok((best_topo.first_location, freq, best_score))
}

pub(crate) fn score_clade_list(
    clades: &[CladeId],
    table: &GlobalCladeTable,
    total_f: f64,
    scoring: MccScoring,
) -> f64 {
    let mut acc = 0.0f64;
    for &cid in clades {
        let freq = (table.clades[cid.index()].count as f64 / total_f).max(f64::MIN_POSITIVE);
        match scoring {
            MccScoring::LogProduct => acc += freq.ln(),
            MccScoring::Sum => acc += freq,
        }
    }
    acc
}

pub(crate) fn hash_clade_list(sorted_clades: &[CladeId]) -> (u64, u64) {
    let mut h1 = 0x517c_c1b7_2722_0a95u64;
    let mut h2 = 0x9e37_79b9_7f4a_7c15u64;
    for &cid in sorted_clades {
        let x = (cid.0 as u64).wrapping_add(1);
        h1 = h1.wrapping_mul(0x9e37_79b1_85eb_ca87) ^ x.rotate_left(17);
        h2 = h2.wrapping_mul(0xc2b2_ae3d_27d4_eb4f) ^ x.rotate_right(13);
    }
    (h1, h2)
}

fn cmp_score(a: f64, b: f64) -> Ordering {
    if (a - b).abs() <= 1e-12 {
        Ordering::Equal
    } else if a > b {
        Ordering::Greater
    } else {
        Ordering::Less
    }
}

pub(crate) type TreeNodeMapping = (Vec<Option<CladeId>>, Vec<Option<usize>>);

/// Map every internal node of an existing `tree` to its `CladeId` in `table` (inserting missing
/// clades with `count = 0` when mapping onto a user-provided `--target` tree).
pub(crate) fn map_tree_nodes_to_clades(
    tree: &Tree,
    table: &mut GlobalCladeTable,
    rooted: bool,
    insert_missing: bool,
) -> Result<TreeNodeMapping> {
    let n_taxa = table.taxa.len();
    let words = table.zobrist.words();
    let n_nodes = tree.len();
    let mut node_bits = vec![0u64; n_nodes * words];
    let mut node_hashes = vec![(0u64, 0u64); n_nodes];
    let mut node_clades = vec![None; n_nodes];
    let mut node_taxa = vec![None; n_nodes];
    let mut seen_taxa = AHashSet::with_capacity(n_taxa);

    let root = tree.root();
    let root_is_binary = !rooted && tree.child_count(root) == 2;
    let first_rc = if root_is_binary {
        tree.first_child(root)
    } else {
        None
    };
    let second_rc = first_rc.and_then(|c| tree.next_sibling(c));

    for id in tree.postorder() {
        let idx = id.index();
        if tree.is_leaf(id) {
            let lbl = tree
                .label(id)
                .ok_or_else(|| anyhow::anyhow!("Tree contains an unlabeled tip"))?;
            let tid = table
                .taxa
                .get_index_of(lbl)
                .map(|t| t.index())
                .ok_or_else(|| anyhow::anyhow!("Taxon '{}' not found in taxon table", lbl))?;
            if !seen_taxa.insert(tid) {
                bail!("Duplicate taxon '{}' in tree", lbl);
            }
            node_taxa[idx] = Some(tid);
            let base = idx * words;
            super::bitset::set_bit(&mut node_bits[base..base + words], tid);
            node_hashes[idx] = table.zobrist.taxon_hash(tid);
        } else {
            let base = idx * words;
            let mut h = (0u64, 0u64);
            for child in tree.children(id) {
                let c_idx = child.index();
                let c_base = c_idx * words;
                for w in 0..words {
                    node_bits[base + w] |= node_bits[c_base + w];
                }
                h.0 ^= node_hashes[c_idx].0;
                h.1 ^= node_hashes[c_idx].1;
            }
            node_hashes[idx] = h;

            if tree.is_root(id) {
                continue;
            }
            if root_is_binary && Some(id) == second_rc {
                continue;
            }
            let pop = popcount(&node_bits[base..base + words]);
            let nontrivial = if rooted {
                pop >= 2 && pop < n_taxa
            } else {
                pop >= 2 && pop <= n_taxa.saturating_sub(2)
            };
            if !nontrivial {
                continue;
            }

            let mut key_bits = node_bits[base..base + words].to_vec();
            let mut key_hash = h;
            if !rooted {
                super::bitset::canonicalize_unrooted(
                    &mut key_bits,
                    &mut key_hash,
                    n_taxa,
                    &table.zobrist,
                );
            }
            if let Some(cid) = table.lookup(&key_bits, key_hash) {
                node_clades[idx] = Some(cid);
            } else if insert_missing {
                let cid = table.insert_zero_clade(&key_bits, key_hash);
                node_clades[idx] = Some(cid);
            }
        }
    }

    if seen_taxa.len() != n_taxa {
        bail!(
            "Tree taxon count ({}) does not match expected taxon set ({})",
            seen_taxa.len(),
            n_taxa
        );
    }

    Ok((node_clades, node_taxa))
}
