//! Maximum Parsimony (Fitch) tree estimation, branch swapping, and MPR branch lengths.

use ahash::AHashMap;
use anyhow::{Result, bail};
use rayon::prelude::*;

use crate::helper::types::{DataType, SeqMatrix};
use crate::parser::tree::{NodeId, Rootedness, Tree};

use super::rng::SeededRng;

/// Search strategy for Maximum Parsimony tree estimation.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum ParsimonySearch {
    /// Exact Branch-and-Bound for `<= 10` taxa; parallel heuristic search otherwise.
    #[default]
    Auto,
    /// Exact Branch-and-Bound search (errors above 20 taxa).
    Exact,
    /// Parallel random stepwise addition followed by branch swapping.
    Heuristic,
}

impl std::fmt::Display for ParsimonySearch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Auto => write!(f, "auto"),
            Self::Exact => write!(f, "exact"),
            Self::Heuristic => write!(f, "heuristic"),
        }
    }
}

impl std::str::FromStr for ParsimonySearch {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "auto" => Ok(Self::Auto),
            "exact" => Ok(Self::Exact),
            "heuristic" => Ok(Self::Heuristic),
            _ => Err(format!("{} is not a valid parsimony search mode", s)),
        }
    }
}

/// Branch-swapping neighborhood used during heuristic Maximum Parsimony search.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum BranchSwap {
    /// Stepwise addition only (no branch swapping).
    None,
    /// Nearest-Neighbor Interchange.
    Nni,
    /// Subtree Pruning and Regrafting.
    Spr,
    /// Tree Bisection and Reconnection (with SPR warm-start).
    #[default]
    Tbr,
}

impl std::fmt::Display for BranchSwap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::None => write!(f, "none"),
            Self::Nni => write!(f, "nni"),
            Self::Spr => write!(f, "spr"),
            Self::Tbr => write!(f, "tbr"),
        }
    }
}

impl std::str::FromStr for BranchSwap {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "none" => Ok(Self::None),
            "nni" => Ok(Self::Nni),
            "spr" => Ok(Self::Spr),
            "tbr" => Ok(Self::Tbr),
            _ => Err(format!("{} is not a valid branch swapping method", s)),
        }
    }
}

/// How gap characters (`-`) are treated in Maximum Parsimony.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum GapMode {
    /// Treat gaps as missing/ambiguous (all states allowed).
    #[default]
    Missing,
    /// Treat gaps as a 5th DNA state (or 21st amino acid state).
    FifthState,
}

impl std::fmt::Display for GapMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Missing => write!(f, "missing"),
            Self::FifthState => write!(f, "fifth"),
        }
    }
}

impl std::str::FromStr for GapMode {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "missing" => Ok(Self::Missing),
            "fifth" | "fifth-state" | "5th" => Ok(Self::FifthState),
            _ => Err(format!("{} is not a valid gap mode", s)),
        }
    }
}

/// Branch length reconstruction mode for Maximum Parsimony trees.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum ParsimonyBranchLengths {
    /// Accelerated transformation (pushes changes toward root).
    #[default]
    Acctran,
    /// Delayed transformation (delays changes toward tips).
    Deltran,
    /// Do not assign branch lengths on parsimony trees.
    None,
}

impl std::fmt::Display for ParsimonyBranchLengths {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Acctran => write!(f, "acctran"),
            Self::Deltran => write!(f, "deltran"),
            Self::None => write!(f, "none"),
        }
    }
}

impl std::str::FromStr for ParsimonyBranchLengths {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "acctran" => Ok(Self::Acctran),
            "deltran" => Ok(Self::Deltran),
            "none" => Ok(Self::None),
            _ => Err(format!("{} is not a valid MP branch length mode", s)),
        }
    }
}

/// Parameters controlling Maximum Parsimony tree estimation.
#[derive(Debug, Clone)]
pub struct ParsimonyOptions {
    pub search: ParsimonySearch,
    pub swap: BranchSwap,
    pub replicates: usize,
    pub gap_mode: GapMode,
    pub branch_lengths: ParsimonyBranchLengths,
    pub seed: u64,
}

impl Default for ParsimonyOptions {
    fn default() -> Self {
        Self {
            search: ParsimonySearch::Auto,
            swap: BranchSwap::Tbr,
            replicates: 10,
            gap_mode: GapMode::Missing,
            branch_lengths: ParsimonyBranchLengths::Acctran,
            seed: 42,
        }
    }
}

/// Output of Maximum Parsimony tree estimation.
#[derive(Debug, Clone)]
pub struct ParsimonyResult {
    pub tree: Tree,
    pub score: usize,
    pub consistency_index: Option<f64>,
    pub retention_index: Option<f64>,
}

const EXACT_AUTO_MAX_TAXA: usize = 10;
const EXACT_WARN_TAXA: usize = 12;
const EXACT_HARD_MAX_TAXA: usize = 20;

const DNA_MASK_TABLE: [u32; 256] = build_dna_masks(false);
const DNA_MASK_TABLE_FIFTH: [u32; 256] = build_dna_masks(true);
const AA_MASK_TABLE: [u32; 256] = build_aa_masks(false);
const AA_MASK_TABLE_FIFTH: [u32; 256] = build_aa_masks(true);

const fn build_dna_masks(fifth_gap: bool) -> [u32; 256] {
    const A: u32 = 1 << 0;
    const C: u32 = 1 << 1;
    const G: u32 = 1 << 2;
    const T: u32 = 1 << 3;
    const GAP: u32 = 1 << 4;
    const ALL_4: u32 = A | C | G | T;
    let gap_mask = if fifth_gap { GAP } else { ALL_4 };
    let missing_mask = if fifth_gap { ALL_4 | GAP } else { ALL_4 };

    let mut table = [missing_mask; 256];
    let pairs = [
        (b'A', A),
        (b'C', C),
        (b'G', G),
        (b'T', T),
        (b'U', T),
        (b'R', A | G),
        (b'Y', C | T),
        (b'S', G | C),
        (b'W', A | T),
        (b'K', G | T),
        (b'M', A | C),
        (b'B', C | G | T),
        (b'D', A | G | T),
        (b'H', A | C | T),
        (b'V', A | C | G),
        (b'N', ALL_4),
        (b'-', gap_mask),
        (b'?', missing_mask),
        (b'.', missing_mask),
    ];
    let mut i = 0;
    while i < pairs.len() {
        let (ch, mask) = pairs[i];
        table[ch as usize] = mask;
        if ch.is_ascii_alphabetic() {
            table[ch.to_ascii_lowercase() as usize] = mask;
        }
        i += 1;
    }
    table
}

const fn build_aa_masks(fifth_gap: bool) -> [u32; 256] {
    const RESIDUES: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";
    const ALL_20: u32 = (1u32 << 20) - 1;
    const GAP: u32 = 1 << 20;
    let gap_mask = if fifth_gap { GAP } else { ALL_20 };
    let missing_mask = if fifth_gap { ALL_20 | GAP } else { ALL_20 };

    let mut table = [missing_mask; 256];
    let mut i = 0;
    while i < RESIDUES.len() {
        let ch = RESIDUES[i];
        let mask = 1u32 << i;
        table[ch as usize] = mask;
        table[ch.to_ascii_lowercase() as usize] = mask;
        i += 1;
    }
    // B = D(3) | N(2), Z = E(6) | Q(5), J = I(8) | L(9)
    let ambig = [
        (b'B', (1u32 << 3) | (1u32 << 2)),
        (b'Z', (1u32 << 6) | (1u32 << 5)),
        (b'J', (1u32 << 8) | (1u32 << 9)),
        (b'X', ALL_20),
        (b'-', gap_mask),
        (b'?', missing_mask),
        (b'.', missing_mask),
        (b'*', missing_mask),
        (b'~', missing_mask),
    ];
    let mut j = 0;
    while j < ambig.len() {
        let (ch, mask) = ambig[j];
        table[ch as usize] = mask;
        if ch.is_ascii_alphabetic() {
            table[ch.to_ascii_lowercase() as usize] = mask;
        }
        j += 1;
    }
    table
}

/// Site-pattern compressed alignment with bitmask states for Fitch parsimony.
#[derive(Debug, Clone)]
pub(crate) struct CompressedAlignment {
    pub(crate) taxa: Vec<String>,
    pub(crate) ntax: usize,
    pub(crate) nchar: usize,
    /// Deduplicated patterns across the full alignment (`taxon * n_all + p`).
    pub(crate) n_all: usize,
    pub(crate) all_weights: Vec<u32>,
    pub(crate) all_states: Vec<u32>,
    /// Mapping from original site column `0..nchar` to pattern index `0..n_all`.
    pub(crate) col_to_pattern: Vec<usize>,
    /// Subset of patterns needed for topology scoring (parsimony-informative or ambiguous variable).
    pub(crate) n_inf: usize,
    pub(crate) inf_weights: Vec<u32>,
    pub(crate) inf_states: Vec<u32>,
    /// Exact constant + unambiguous-singleton step offset independent of topology.
    pub(crate) const_score: usize,
    /// Minimum possible steps (`m`) summed over all sites for CI/RI.
    pub(crate) min_steps: usize,
    /// Maximum possible steps (`g`) summed over all sites for RI.
    pub(crate) max_steps: usize,
}

impl CompressedAlignment {
    /// Compress an aligned [`SeqMatrix`] into weighted bitmask site patterns.
    pub(crate) fn from_seq_matrix(
        matrix: &SeqMatrix,
        datatype: &DataType,
        gap_mode: GapMode,
    ) -> Result<Self> {
        let ntax = matrix.len();
        if ntax < 3 {
            bail!(
                "At least 3 taxa are required for tree estimation (found {})",
                ntax
            );
        }
        if matches!(datatype, DataType::Ignore) {
            bail!("Tree estimation requires DNA or amino acid datatype, not 'ignore'");
        }
        let table = match (datatype, gap_mode) {
            (DataType::Dna, GapMode::Missing) => &DNA_MASK_TABLE,
            (DataType::Dna, GapMode::FifthState) => &DNA_MASK_TABLE_FIFTH,
            (DataType::Aa, GapMode::Missing) => &AA_MASK_TABLE,
            (DataType::Aa, GapMode::FifthState) => &AA_MASK_TABLE_FIFTH,
            (DataType::Ignore, _) => unreachable!(),
        };

        let mut nchar = None;
        for (id, seq) in matrix.iter() {
            if seq.is_empty() {
                bail!("Sequence '{}' is empty", id);
            }
            match nchar {
                None => nchar = Some(seq.len()),
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
        let nchar = nchar.unwrap_or(0);
        let taxa: Vec<String> = matrix.keys().cloned().collect();
        let rows: Vec<&[u8]> = matrix.values().collect();

        let mut pattern_map: AHashMap<Vec<u32>, usize> = AHashMap::new();
        let mut unique_cols: Vec<Vec<u32>> = Vec::new();
        let mut all_weights: Vec<u32> = Vec::new();
        let mut col_to_pattern = Vec::with_capacity(nchar);
        let mut col_buf = vec![0u32; ntax];

        for site in 0..nchar {
            for (t, row) in rows.iter().enumerate() {
                col_buf[t] = table[row[site] as usize];
            }
            if let Some(&pid) = pattern_map.get(&col_buf) {
                all_weights[pid] += 1;
                col_to_pattern.push(pid);
            } else {
                let pid = unique_cols.len();
                let owned = col_buf.clone();
                pattern_map.insert(owned.clone(), pid);
                unique_cols.push(owned);
                all_weights.push(1);
                col_to_pattern.push(pid);
            }
        }

        let n_all = unique_cols.len();
        let mut all_states = vec![0u32; ntax * n_all];
        for (p, col) in unique_cols.iter().enumerate() {
            for (t, &state) in col.iter().enumerate() {
                all_states[t * n_all + p] = state;
            }
        }

        Ok(Self::build_summary(
            taxa,
            ntax,
            nchar,
            n_all,
            all_weights,
            all_states,
            col_to_pattern,
        ))
    }

    /// Build a resampled alignment using bootstrap weights over `0..n_all`.
    pub(crate) fn with_resampled_weights(&self, all_weights: Vec<u32>) -> Self {
        Self::build_summary(
            self.taxa.clone(),
            self.ntax,
            self.nchar,
            self.n_all,
            all_weights,
            self.all_states.clone(),
            self.col_to_pattern.clone(),
        )
    }

    fn build_summary(
        taxa: Vec<String>,
        ntax: usize,
        nchar: usize,
        n_all: usize,
        all_weights: Vec<u32>,
        all_states: Vec<u32>,
        col_to_pattern: Vec<usize>,
    ) -> Self {
        let mut const_score = 0usize;
        let mut min_steps = 0usize;
        let mut max_steps = 0usize;
        let mut inf_indices = Vec::new();

        for (p, &w) in all_weights.iter().enumerate() {
            if w == 0 {
                continue;
            }
            let w_usize = w as usize;
            let mut all_unambiguous = true;
            let mut singleton_counts = [0u32; 32];
            let mut any_union = 0u32;
            let mut exact_union = 0u32;
            let mut exact_intersection = u32::MAX;

            for t in 0..ntax {
                let s = all_states[t * n_all + p];
                any_union |= s;
                if s.count_ones() == 1 {
                    let bit = s.trailing_zeros() as usize;
                    singleton_counts[bit] += 1;
                    exact_union |= s;
                    exact_intersection &= s;
                } else {
                    all_unambiguous = false;
                }
            }

            if all_unambiguous {
                let distinct = exact_union.count_ones() as usize;
                let site_min = distinct.saturating_sub(1);
                let max_freq = singleton_counts.iter().copied().max().unwrap_or(0) as usize;
                let site_max = ntax.saturating_sub(max_freq);
                min_steps += site_min * w_usize;
                max_steps += site_max * w_usize;

                let mut multi_states = 0usize;
                for &c in &singleton_counts {
                    if c >= 2 {
                        multi_states += 1;
                    }
                }
                if multi_states <= 1 {
                    const_score += site_min * w_usize;
                } else {
                    inf_indices.push(p);
                }
            } else {
                // Estimate site min/max steps with ambiguity present:
                // minimum steps from unambiguous subset if non-empty (exact lower bound).
                let site_min = if exact_intersection != 0 && exact_union != 0 {
                    // Check if all masks can share a single state
                    let mut common = u32::MAX;
                    for t in 0..ntax {
                        common &= all_states[t * n_all + p];
                    }
                    if common != 0 {
                        0
                    } else {
                        (exact_union.count_ones() as usize).saturating_sub(1).max(1)
                    }
                } else if exact_union != 0 {
                    let mut common = u32::MAX;
                    for t in 0..ntax {
                        common &= all_states[t * n_all + p];
                    }
                    if common != 0 {
                        0
                    } else {
                        (exact_union.count_ones() as usize).saturating_sub(1).max(1)
                    }
                } else {
                    let mut common = u32::MAX;
                    for t in 0..ntax {
                        common &= all_states[t * n_all + p];
                    }
                    usize::from(common == 0)
                };
                // Upper bound for RI: total unambiguous taxa minus max single-state count,
                // or conservative star-tree cost.
                let mut unambig_total = 0usize;
                let mut best_count = 0usize;
                for &c in &singleton_counts {
                    unambig_total += c as usize;
                    best_count = best_count.max(c as usize);
                }
                let _ = any_union;
                let site_max = unambig_total.saturating_sub(best_count);
                min_steps += site_min * w_usize;
                max_steps += site_max.max(site_min) * w_usize;

                // Keep variable ambiguous patterns unless all taxa share at least one state.
                let mut common = u32::MAX;
                for t in 0..ntax {
                    common &= all_states[t * n_all + p];
                }
                if common == 0 {
                    inf_indices.push(p);
                }
            }
        }

        let n_inf = inf_indices.len();
        let mut inf_weights = Vec::with_capacity(n_inf);
        let mut inf_states = vec![0u32; ntax * n_inf];
        for (new_p, &old_p) in inf_indices.iter().enumerate() {
            inf_weights.push(all_weights[old_p]);
            for t in 0..ntax {
                inf_states[t * n_inf + new_p] = all_states[t * n_all + old_p];
            }
        }

        Self {
            taxa,
            ntax,
            nchar,
            n_all,
            all_weights,
            all_states,
            col_to_pattern,
            n_inf,
            inf_weights,
            inf_states,
            const_score,
            min_steps,
            max_steps,
        }
    }

    /// Evaluate exact uncompressed-equivalent Fitch score on a [`Tree`] directly.
    #[cfg(test)]
    pub(crate) fn score_tree_direct(&self, tree: &Tree) -> usize {
        let mut taxon_map = AHashMap::with_capacity(self.ntax);
        for (i, name) in self.taxa.iter().enumerate() {
            taxon_map.insert(name.as_str(), i);
        }
        let n_nodes = tree.len();
        let mut states = vec![0u32; n_nodes * self.n_all];
        let mut total = 0usize;

        for id in tree.postorder() {
            if tree.is_leaf(id) {
                let tid = taxon_map[tree.label(id).expect("Leaf missing label")];
                let dst = &mut states[id.index() * self.n_all..(id.index() + 1) * self.n_all];
                dst.copy_from_slice(&self.all_states[tid * self.n_all..(tid + 1) * self.n_all]);
            } else {
                let children: Vec<NodeId> = tree.children(id).collect();
                for p in 0..self.n_all {
                    let w = self.all_weights[p] as usize;
                    if w == 0 {
                        continue;
                    }
                    let mut cur = states[children[0].index() * self.n_all + p];
                    for &ch in &children[1..] {
                        let cs = states[ch.index() * self.n_all + p];
                        let inter = cur & cs;
                        if inter != 0 {
                            cur = inter;
                        } else {
                            cur |= cs;
                            total += w;
                        }
                    }
                    states[id.index() * self.n_all + p] = cur;
                }
            }
        }
        total
    }
}

/// Compact unrooted binary tree over `0..n_taxa` leaves and `n_taxa..2*n_taxa-2` internals.
#[derive(Debug, Clone)]
struct UnrootedTree {
    n_taxa: usize,
    active_internals: usize,
    leaf_nbr: Vec<usize>,
    internal_nbrs: Vec<[usize; 3]>,
}

impl UnrootedTree {
    fn star(n_taxa: usize, t0: usize, t1: usize, t2: usize) -> Self {
        let mut leaf_nbr = vec![usize::MAX; n_taxa];
        let mut internal_nbrs = vec![[usize::MAX; 3]; n_taxa.saturating_sub(2)];
        let root_int = n_taxa;
        leaf_nbr[t0] = root_int;
        leaf_nbr[t1] = root_int;
        leaf_nbr[t2] = root_int;
        internal_nbrs[0] = [t0, t1, t2];
        Self {
            n_taxa,
            active_internals: 1,
            leaf_nbr,
            internal_nbrs,
        }
    }

    #[inline]
    fn is_leaf(&self, node: usize) -> bool {
        node < self.n_taxa
    }

    #[inline]
    fn neighbors(&self, node: usize) -> &[usize] {
        if node < self.n_taxa {
            std::slice::from_ref(&self.leaf_nbr[node])
        } else {
            &self.internal_nbrs[node - self.n_taxa]
        }
    }

    #[inline]
    fn set_neighbor(&mut self, node: usize, old: usize, new: usize) {
        if node < self.n_taxa {
            debug_assert_eq!(self.leaf_nbr[node], old);
            self.leaf_nbr[node] = new;
        } else {
            let slot = &mut self.internal_nbrs[node - self.n_taxa];
            for item in slot.iter_mut() {
                if *item == old {
                    *item = new;
                    return;
                }
            }
            unreachable!("Old neighbor not found in internal node");
        }
    }

    #[inline]
    fn neighbor_slot(&self, node: usize, target: usize) -> usize {
        if node < self.n_taxa {
            debug_assert_eq!(self.leaf_nbr[node], target);
            0
        } else {
            let nbrs = &self.internal_nbrs[node - self.n_taxa];
            if nbrs[0] == target {
                0
            } else if nbrs[1] == target {
                1
            } else {
                debug_assert_eq!(nbrs[2], target);
                2
            }
        }
    }

    /// Insert `taxon` onto undirected edge `(u, v)` using next internal slot.
    fn insert_on_edge(&mut self, u: usize, v: usize, taxon: usize) -> usize {
        let new_int = self.n_taxa + self.active_internals;
        self.active_internals += 1;
        self.set_neighbor(u, v, new_int);
        self.set_neighbor(v, u, new_int);
        self.internal_nbrs[new_int - self.n_taxa] = [u, v, taxon];
        self.leaf_nbr[taxon] = new_int;
        new_int
    }

    /// Undo the most recent `insert_on_edge` on `(u, v)`.
    fn undo_insert(&mut self, u: usize, v: usize, taxon: usize, new_int: usize) {
        self.set_neighbor(u, new_int, v);
        self.set_neighbor(v, new_int, u);
        self.leaf_nbr[taxon] = usize::MAX;
        self.internal_nbrs[new_int - self.n_taxa] = [usize::MAX; 3];
        self.active_internals -= 1;
    }

    /// Collect all current undirected edges `(u, v)` (`u < v`).
    fn collect_edges(&self, out: &mut Vec<(usize, usize)>) {
        out.clear();
        for idx in 0..self.active_internals {
            let u = self.n_taxa + idx;
            for &v in &self.internal_nbrs[idx] {
                if u < v || v < self.n_taxa {
                    out.push((u, v));
                }
            }
        }
    }
}

/// Directed edge Fitch state cache for `UnrootedTree`.
///
/// For any node `u` and neighbor slot `s` (0 for leaf, 0..3 for internal),
/// stores the Fitch state vector and accumulated subtree cost looking away
/// from neighbor `s`.
struct FitchCache {
    n_taxa: usize,
    n_pat: usize,
    states: Vec<u32>,
    costs: Vec<usize>,
    order: Vec<(usize, usize)>,
}

impl FitchCache {
    fn new(n_taxa: usize, n_pat: usize) -> Self {
        let dir_count = n_taxa + n_taxa.saturating_sub(2) * 3;
        Self {
            n_taxa,
            n_pat,
            states: vec![0u32; dir_count * n_pat],
            costs: vec![0usize; dir_count],
            order: Vec::with_capacity(2 * n_taxa),
        }
    }

    #[inline]
    fn dir_index(&self, node: usize, slot: usize) -> usize {
        if node < self.n_taxa {
            node
        } else {
            self.n_taxa + (node - self.n_taxa) * 3 + slot
        }
    }

    #[inline]
    fn state_slice(&self, dir_idx: usize) -> &[u32] {
        let start = dir_idx * self.n_pat;
        &self.states[start..start + self.n_pat]
    }

    /// Populate all directed edge states for `tree` using pattern states & weights,
    /// returning the tree's informative-pattern score (excluding `const_score`).
    fn compute(&mut self, tree: &UnrootedTree, pat_states: &[u32], pat_weights: &[u32]) -> usize {
        let n_pat = self.n_pat;
        let root = tree.n_taxa;
        self.order.clear();
        let mut stack = vec![(root, usize::MAX)];
        while let Some((u, parent)) = stack.pop() {
            self.order.push((u, parent));
            for &v in tree.neighbors(u) {
                if v != parent {
                    stack.push((v, u));
                }
            }
        }

        // Postorder pass: compute `node -> parent` for every non-root node.
        for idx in (1..self.order.len()).rev() {
            let (u, parent) = self.order[idx];
            let p_slot = tree.neighbor_slot(u, parent);
            let out_dir = self.dir_index(u, p_slot);

            if tree.is_leaf(u) {
                let src_start = u * n_pat;
                let dst_start = out_dir * n_pat;
                self.states[dst_start..dst_start + n_pat]
                    .copy_from_slice(&pat_states[src_start..src_start + n_pat]);
                self.costs[out_dir] = 0;
            } else {
                let nbrs = tree.internal_nbrs[u - tree.n_taxa];
                let c1 = nbrs[(p_slot + 1) % 3];
                let c2 = nbrs[(p_slot + 2) % 3];
                let d1 = self.dir_index(c1, tree.neighbor_slot(c1, u));
                let d2 = self.dir_index(c2, tree.neighbor_slot(c2, u));
                self.costs[out_dir] = combine_dir_pair(
                    &mut self.states,
                    out_dir,
                    d1,
                    d2,
                    self.costs[d1],
                    self.costs[d2],
                    n_pat,
                    pat_weights,
                );
            }
        }

        // Preorder pass: compute `u -> child` for every node including root.
        for idx in 0..self.order.len() {
            let (u, _) = self.order[idx];
            if tree.is_leaf(u) {
                continue;
            }
            let nbrs = tree.internal_nbrs[u - tree.n_taxa];
            for slot in 0..3 {
                let child = nbrs[slot];
                if idx > 0 && child == self.order[idx].1 {
                    continue;
                }
                let o1 = nbrs[(slot + 1) % 3];
                let o2 = nbrs[(slot + 2) % 3];
                let d1 = self.dir_index(o1, tree.neighbor_slot(o1, u));
                let d2 = self.dir_index(o2, tree.neighbor_slot(o2, u));
                let out_dir = self.dir_index(u, slot);
                self.costs[out_dir] = combine_dir_pair(
                    &mut self.states,
                    out_dir,
                    d1,
                    d2,
                    self.costs[d1],
                    self.costs[d2],
                    n_pat,
                    pat_weights,
                );
            }
        }

        let first_child = tree.internal_nbrs[0][0];
        let d_root = self.dir_index(root, 0);
        let d_child = self.dir_index(first_child, tree.neighbor_slot(first_child, root));
        self.edge_total_score(d_root, d_child, pat_weights)
    }

    #[inline]
    fn edge_total_score(&self, d_u: usize, d_v: usize, pat_weights: &[u32]) -> usize {
        let su = self.state_slice(d_u);
        let sv = self.state_slice(d_v);
        let mut diff = 0usize;
        for p in 0..self.n_pat {
            if (su[p] & sv[p]) == 0 {
                diff += pat_weights[p] as usize;
            }
        }
        self.costs[d_u] + self.costs[d_v] + diff
    }

    /// Evaluate the extra steps incurred when inserting `leaf` on edge `(u, v)`.
    #[inline]
    fn eval_insert_extra(
        &self,
        tree: &UnrootedTree,
        u: usize,
        v: usize,
        leaf_states: &[u32],
        pat_weights: &[u32],
    ) -> usize {
        let du = self.dir_index(u, tree.neighbor_slot(u, v));
        let dv = self.dir_index(v, tree.neighbor_slot(v, u));
        let su = self.state_slice(du);
        let sv = self.state_slice(dv);
        let mut extra = 0usize;
        for p in 0..self.n_pat {
            let inter = su[p] & sv[p];
            let comb = if inter != 0 { inter } else { su[p] | sv[p] };
            if (comb & leaf_states[p]) == 0 {
                extra += pat_weights[p] as usize;
            }
        }
        extra
    }
}

#[allow(clippy::too_many_arguments)]
#[inline]
fn combine_dir_pair(
    states: &mut [u32],
    out_dir: usize,
    d1: usize,
    d2: usize,
    cost1: usize,
    cost2: usize,
    n_pat: usize,
    pat_weights: &[u32],
) -> usize {
    let mut cost = cost1 + cost2;
    let o_base = out_dir * n_pat;
    let s1_base = d1 * n_pat;
    let s2_base = d2 * n_pat;
    for p in 0..n_pat {
        let a = states[s1_base + p];
        let b = states[s2_base + p];
        let inter = a & b;
        if inter != 0 {
            states[o_base + p] = inter;
        } else {
            states[o_base + p] = a | b;
            cost += pat_weights[p] as usize;
        }
    }
    cost
}

#[inline]
fn fitch_pair_cost(s1: &[u32], s2: &[u32], pat_weights: &[u32], out_comb: &mut [u32]) -> usize {
    let mut cost = 0usize;
    for p in 0..pat_weights.len() {
        let a = s1[p];
        let b = s2[p];
        let inter = a & b;
        if inter != 0 {
            out_comb[p] = inter;
        } else {
            out_comb[p] = a | b;
            cost += pat_weights[p] as usize;
        }
    }
    cost
}

#[inline]
fn fitch_diff_cost(s1: &[u32], s2: &[u32], pat_weights: &[u32]) -> usize {
    let mut cost = 0usize;
    for p in 0..pat_weights.len() {
        if (s1[p] & s2[p]) == 0 {
            cost += pat_weights[p] as usize;
        }
    }
    cost
}

/// Estimate the Maximum Parsimony tree from a compressed alignment.
pub(crate) fn estimate_parsimony_from_compressed(
    aln: &CompressedAlignment,
    opts: &ParsimonyOptions,
    with_branch_lengths: bool,
) -> Result<ParsimonyResult> {
    let n = aln.ntax;
    if n < 3 {
        bail!("At least 3 taxa are required for parsimony (found {})", n);
    }

    let use_exact = match opts.search {
        ParsimonySearch::Auto => n <= EXACT_AUTO_MAX_TAXA,
        ParsimonySearch::Exact => {
            if n > EXACT_HARD_MAX_TAXA {
                bail!(
                    "Exact branch-and-bound search supports up to {} taxa (found {}). \
                    Use --search heuristic instead.",
                    EXACT_HARD_MAX_TAXA,
                    n
                );
            }
            if n > EXACT_WARN_TAXA {
                log::warn!(
                    "Exact branch-and-bound on {} taxa may take significant time; \
                    consider --search heuristic.",
                    n
                );
            }
            true
        }
        ParsimonySearch::Heuristic => false,
    };

    let (best_unrooted, inf_score) = if n == 3 {
        let tree = UnrootedTree::star(3, 0, 1, 2);
        let mut cache = FitchCache::new(3, aln.n_inf);
        let s = cache.compute(&tree, &aln.inf_states, &aln.inf_weights);
        (tree, s)
    } else if use_exact {
        exact_branch_and_bound(aln, opts)
    } else {
        heuristic_search(aln, opts)
    };

    let total_score = aln.const_score + inf_score;
    let ci = if total_score > 0 {
        Some((aln.min_steps as f64 / total_score as f64).min(1.0))
    } else {
        None
    };
    let ri = if aln.max_steps > aln.min_steps {
        Some(
            ((aln.max_steps as f64 - total_score as f64)
                / (aln.max_steps as f64 - aln.min_steps as f64))
                .clamp(0.0, 1.0),
        )
    } else {
        None
    };

    let bl_mode = if with_branch_lengths {
        opts.branch_lengths
    } else {
        ParsimonyBranchLengths::None
    };
    let tree = unrooted_to_tree(&best_unrooted, aln, bl_mode);

    Ok(ParsimonyResult {
        tree,
        score: total_score,
        consistency_index: ci,
        retention_index: ri,
    })
}

/// Estimate Maximum Parsimony tree directly from a [`SeqMatrix`].
pub fn estimate_parsimony(
    matrix: &SeqMatrix,
    datatype: &DataType,
    opts: &ParsimonyOptions,
) -> Result<ParsimonyResult> {
    let aln = CompressedAlignment::from_seq_matrix(matrix, datatype, opts.gap_mode)?;
    estimate_parsimony_from_compressed(&aln, opts, true)
}

fn heuristic_search(aln: &CompressedAlignment, opts: &ParsimonyOptions) -> (UnrootedTree, usize) {
    let reps = opts.replicates.max(1);
    let results: Vec<(usize, usize, UnrootedTree)> = if reps == 1 {
        vec![run_single_replicate(aln, opts, 0)]
    } else {
        (0..reps)
            .into_par_iter()
            .map(|rep| run_single_replicate(aln, opts, rep))
            .collect()
    };

    let mut best = results.into_iter().min_by_key(|(s, rep, _)| (*s, *rep));
    let (score, _, tree) = best.take().expect("At least 1 replicate runs");
    (tree, score)
}

fn run_single_replicate(
    aln: &CompressedAlignment,
    opts: &ParsimonyOptions,
    rep: usize,
) -> (usize, usize, UnrootedTree) {
    let mut rng = SeededRng::derive(opts.seed, rep as u64);
    let mut order: Vec<usize> = (0..aln.ntax).collect();
    rng.shuffle(&mut order);

    let (mut tree, mut score) = stepwise_addition(aln, &order);
    match opts.swap {
        BranchSwap::None => {}
        BranchSwap::Nni => score = nni_search(&mut tree, aln),
        BranchSwap::Spr => score = spr_search(&mut tree, aln),
        BranchSwap::Tbr => {
            let _ = spr_search(&mut tree, aln);
            score = tbr_search(&mut tree, aln);
        }
    }
    (score, rep, tree)
}

fn stepwise_addition(aln: &CompressedAlignment, order: &[usize]) -> (UnrootedTree, usize) {
    let mut tree = UnrootedTree::star(aln.ntax, order[0], order[1], order[2]);
    let mut cache = FitchCache::new(aln.ntax, aln.n_inf);
    let mut edges = Vec::with_capacity(2 * aln.ntax);
    let mut score = cache.compute(&tree, &aln.inf_states, &aln.inf_weights);

    for &taxon in &order[3..] {
        let leaf_states = &aln.inf_states[taxon * aln.n_inf..(taxon + 1) * aln.n_inf];
        tree.collect_edges(&mut edges);
        let mut best_edge = edges[0];
        let mut best_extra = usize::MAX;
        for &(u, v) in &edges {
            let extra = cache.eval_insert_extra(&tree, u, v, leaf_states, &aln.inf_weights);
            if extra < best_extra {
                best_extra = extra;
                best_edge = (u, v);
                if extra == 0 {
                    break;
                }
            }
        }
        tree.insert_on_edge(best_edge.0, best_edge.1, taxon);
        score = cache.compute(&tree, &aln.inf_states, &aln.inf_weights);
    }

    (tree, score)
}

fn nni_search(tree: &mut UnrootedTree, aln: &CompressedAlignment) -> usize {
    let mut cache = FitchCache::new(aln.ntax, aln.n_inf);
    let mut current_score = cache.compute(tree, &aln.inf_states, &aln.inf_weights);
    let mut edges = Vec::with_capacity(2 * aln.ntax);
    let mut comb1 = vec![0u32; aln.n_inf];
    let mut comb2 = vec![0u32; aln.n_inf];

    loop {
        tree.collect_edges(&mut edges);
        let mut best_move: Option<(usize, usize, usize, usize)> = None;
        let mut best_score = current_score;

        for &(u, v) in &edges {
            if tree.is_leaf(u) || tree.is_leaf(v) {
                continue;
            }
            let u_nbrs = tree.internal_nbrs[u - tree.n_taxa];
            let v_nbrs = tree.internal_nbrs[v - tree.n_taxa];
            let u_others: [usize; 2] = other_two(u_nbrs, v);
            let v_others: [usize; 2] = other_two(v_nbrs, u);

            let a = u_others[0];
            let b = u_others[1];
            let c = v_others[0];
            let d = v_others[1];

            let da = cache.dir_index(a, tree.neighbor_slot(a, u));
            let db = cache.dir_index(b, tree.neighbor_slot(b, u));
            let dc = cache.dir_index(c, tree.neighbor_slot(c, v));
            let dd = cache.dir_index(d, tree.neighbor_slot(d, v));

            let base_cost = cache.costs[da] + cache.costs[db] + cache.costs[dc] + cache.costs[dd];

            // Swap 1: (a, c) at u, (b, d) at v
            let c1 = fitch_pair_cost(
                cache.state_slice(da),
                cache.state_slice(dc),
                &aln.inf_weights,
                &mut comb1,
            );
            let c2 = fitch_pair_cost(
                cache.state_slice(db),
                cache.state_slice(dd),
                &aln.inf_weights,
                &mut comb2,
            );
            let s1 = base_cost + c1 + c2 + fitch_diff_cost(&comb1, &comb2, &aln.inf_weights);
            if s1 < best_score {
                best_score = s1;
                best_move = Some((u, v, b, c));
            }

            // Swap 2: (a, d) at u, (c, b) at v
            let c1 = fitch_pair_cost(
                cache.state_slice(da),
                cache.state_slice(dd),
                &aln.inf_weights,
                &mut comb1,
            );
            let c2 = fitch_pair_cost(
                cache.state_slice(dc),
                cache.state_slice(db),
                &aln.inf_weights,
                &mut comb2,
            );
            let s2 = base_cost + c1 + c2 + fitch_diff_cost(&comb1, &comb2, &aln.inf_weights);
            if s2 < best_score {
                best_score = s2;
                best_move = Some((u, v, b, d));
            }
        }

        if let Some((u, v, u_child, v_child)) = best_move {
            tree.set_neighbor(u, u_child, v_child);
            tree.set_neighbor(v, v_child, u_child);
            tree.set_neighbor(u_child, u, v);
            tree.set_neighbor(v_child, v, u);
            current_score = cache.compute(tree, &aln.inf_states, &aln.inf_weights);
        } else {
            break;
        }
    }

    current_score
}

fn spr_search(tree: &mut UnrootedTree, aln: &CompressedAlignment) -> usize {
    let mut cache = FitchCache::new(aln.ntax, aln.n_inf);
    let mut current_score = cache.compute(tree, &aln.inf_states, &aln.inf_weights);
    let mut edges = Vec::with_capacity(2 * aln.ntax);

    loop {
        let mut best_move: Option<(usize, usize, usize, usize)> = None;
        let mut best_score = current_score;

        let n_int = tree.active_internals;
        for idx in 0..n_int {
            let p = tree.n_taxa + idx;
            let p_nbrs = tree.internal_nbrs[idx];
            for slot in 0..3 {
                let sub = p_nbrs[slot];
                let a = p_nbrs[(slot + 1) % 3];
                let b = p_nbrs[(slot + 2) % 3];
                if tree.is_leaf(a) && tree.is_leaf(b) {
                    continue;
                }

                let sub_dir = cache.dir_index(sub, tree.neighbor_slot(sub, p));
                let sub_cost = cache.costs[sub_dir];
                let sub_states = cache.state_slice(sub_dir).to_vec();

                // Prune `sub` and suppress degree-2 node `p` by linking `(a, b)`.
                let mut pruned = tree.clone();
                pruned.set_neighbor(a, p, b);
                pruned.set_neighbor(b, p, a);
                // Root the remaining component at an internal node on the `(a, b)` side.
                let pruned_root = if !pruned.is_leaf(a) { a } else { b };
                let (pruned_score, pruned_cache) = compute_subcomponent_cache(
                    &pruned,
                    pruned_root,
                    &aln.inf_states,
                    &aln.inf_weights,
                );

                collect_reachable_edges(&pruned, pruned_root, &mut edges);
                for &(x, y) in &edges {
                    if (x == a && y == b) || (x == b && y == a) {
                        continue;
                    }
                    let extra = pruned_cache.eval_insert_extra(
                        &pruned,
                        x,
                        y,
                        &sub_states,
                        &aln.inf_weights,
                    );
                    let total = pruned_score + sub_cost + extra;
                    if total < best_score {
                        best_score = total;
                        best_move = Some((p, sub, x, y));
                    }
                }
            }
        }

        if let Some((p, sub, x, y)) = best_move {
            let p_nbrs = tree.internal_nbrs[p - tree.n_taxa];
            let [a, b] = other_two(p_nbrs, sub);
            // Detach p from a, b
            tree.set_neighbor(a, p, b);
            tree.set_neighbor(b, p, a);
            // Re-insert p on (x, y)
            tree.set_neighbor(x, y, p);
            tree.set_neighbor(y, x, p);
            tree.internal_nbrs[p - tree.n_taxa] = [x, y, sub];
            current_score = cache.compute(tree, &aln.inf_states, &aln.inf_weights);
        } else {
            break;
        }
    }

    current_score
}

fn tbr_search(tree: &mut UnrootedTree, aln: &CompressedAlignment) -> usize {
    let mut cache = FitchCache::new(aln.ntax, aln.n_inf);
    let mut current_score = cache.compute(tree, &aln.inf_states, &aln.inf_weights);
    if aln.ntax <= 4 {
        return current_score;
    }

    let mut edges = Vec::with_capacity(2 * aln.ntax);
    let mut comp1_edges = Vec::with_capacity(2 * aln.ntax);
    let mut comp2_edges = Vec::with_capacity(2 * aln.ntax);

    loop {
        tree.collect_edges(&mut edges);
        let mut best_move: Option<(usize, usize, usize, usize, usize, usize)> = None;
        let mut best_score = current_score;

        for &(u, v) in &edges {
            if tree.is_leaf(u) || tree.is_leaf(v) {
                continue;
            }
            let [a, b] = other_two(tree.internal_nbrs[u - tree.n_taxa], v);
            let [c, d] = other_two(tree.internal_nbrs[v - tree.n_taxa], u);

            let mut split_tree = tree.clone();
            split_tree.set_neighbor(a, u, b);
            split_tree.set_neighbor(b, u, a);
            split_tree.set_neighbor(c, v, d);
            split_tree.set_neighbor(d, v, c);

            let (s1, e1_roots) = component_edge_roots(
                &split_tree,
                a,
                b,
                &aln.inf_states,
                &aln.inf_weights,
                &mut comp1_edges,
            );
            let (s2, e2_roots) = component_edge_roots(
                &split_tree,
                c,
                d,
                &aln.inf_states,
                &aln.inf_weights,
                &mut comp2_edges,
            );

            for (i1, &(x1, y1)) in comp1_edges.iter().enumerate() {
                let r1 = &e1_roots[i1 * aln.n_inf..(i1 + 1) * aln.n_inf];
                for (i2, &(x2, y2)) in comp2_edges.iter().enumerate() {
                    if ((x1 == a && y1 == b) || (x1 == b && y1 == a))
                        && ((x2 == c && y2 == d) || (x2 == d && y2 == c))
                    {
                        continue;
                    }
                    let r2 = &e2_roots[i2 * aln.n_inf..(i2 + 1) * aln.n_inf];
                    let bridge_cost = fitch_diff_cost(r1, r2, &aln.inf_weights);
                    let total = s1 + s2 + bridge_cost;
                    if total < best_score {
                        best_score = total;
                        best_move = Some((u, v, x1, y1, x2, y2));
                    }
                }
            }
        }

        if let Some((u, v, x1, y1, x2, y2)) = best_move {
            let [a, b] = other_two(tree.internal_nbrs[u - tree.n_taxa], v);
            let [c, d] = other_two(tree.internal_nbrs[v - tree.n_taxa], u);

            // Suppress old u and v
            tree.set_neighbor(a, u, b);
            tree.set_neighbor(b, u, a);
            tree.set_neighbor(c, v, d);
            tree.set_neighbor(d, v, c);

            // Reattach u on (x1, y1) and v on (x2, y2) with edge (u, v)
            tree.set_neighbor(x1, y1, u);
            tree.set_neighbor(y1, x1, u);
            tree.internal_nbrs[u - tree.n_taxa] = [x1, y1, v];

            tree.set_neighbor(x2, y2, v);
            tree.set_neighbor(y2, x2, v);
            tree.internal_nbrs[v - tree.n_taxa] = [x2, y2, u];

            current_score = cache.compute(tree, &aln.inf_states, &aln.inf_weights);
        } else {
            break;
        }
    }

    current_score
}

/// Compute component score and combined Fitch root states for every edge in a bisected component.
fn component_edge_roots(
    tree: &UnrootedTree,
    seed_a: usize,
    seed_b: usize,
    pat_states: &[u32],
    pat_weights: &[u32],
    comp_edges: &mut Vec<(usize, usize)>,
) -> (usize, Vec<u32>) {
    let n_pat = pat_weights.len();
    if tree.is_leaf(seed_a) && tree.is_leaf(seed_b) {
        comp_edges.clear();
        comp_edges.push((seed_a, seed_b));
        let mut root_state = vec![0u32; n_pat];
        let sa = &pat_states[seed_a * n_pat..(seed_a + 1) * n_pat];
        let sb = &pat_states[seed_b * n_pat..(seed_b + 1) * n_pat];
        let score = fitch_pair_cost(sa, sb, pat_weights, &mut root_state);
        return (score, root_state);
    }

    let root = if !tree.is_leaf(seed_a) {
        seed_a
    } else {
        seed_b
    };
    let (score, cache) = compute_subcomponent_cache(tree, root, pat_states, pat_weights);
    collect_reachable_edges(tree, root, comp_edges);
    let mut roots = vec![0u32; comp_edges.len() * n_pat];
    for (ei, &(x, y)) in comp_edges.iter().enumerate() {
        let dx = cache.dir_index(x, tree.neighbor_slot(x, y));
        let dy = cache.dir_index(y, tree.neighbor_slot(y, x));
        let sx = cache.state_slice(dx);
        let sy = cache.state_slice(dy);
        let dst = &mut roots[ei * n_pat..(ei + 1) * n_pat];
        for p in 0..n_pat {
            let inter = sx[p] & sy[p];
            dst[p] = if inter != 0 { inter } else { sx[p] | sy[p] };
        }
    }
    (score, roots)
}

fn collect_reachable_edges(tree: &UnrootedTree, start: usize, out: &mut Vec<(usize, usize)>) {
    out.clear();
    let mut stack = vec![(start, usize::MAX)];
    while let Some((u, parent)) = stack.pop() {
        for &v in tree.neighbors(u) {
            if v != parent {
                out.push((u, v));
                stack.push((v, u));
            }
        }
    }
}

fn compute_subcomponent_cache(
    tree: &UnrootedTree,
    root: usize,
    pat_states: &[u32],
    pat_weights: &[u32],
) -> (usize, FitchCache) {
    let n_pat = pat_weights.len();
    let mut cache = FitchCache::new(tree.n_taxa, n_pat);
    cache.order.clear();
    let mut stack = vec![(root, usize::MAX)];
    while let Some((u, parent)) = stack.pop() {
        cache.order.push((u, parent));
        for &v in tree.neighbors(u) {
            if v != parent {
                stack.push((v, u));
            }
        }
    }

    for idx in (1..cache.order.len()).rev() {
        let (u, parent) = cache.order[idx];
        let p_slot = tree.neighbor_slot(u, parent);
        let out_dir = cache.dir_index(u, p_slot);
        if tree.is_leaf(u) {
            let src_start = u * n_pat;
            let dst_start = out_dir * n_pat;
            cache.states[dst_start..dst_start + n_pat]
                .copy_from_slice(&pat_states[src_start..src_start + n_pat]);
            cache.costs[out_dir] = 0;
        } else {
            let nbrs = tree.internal_nbrs[u - tree.n_taxa];
            let c1 = nbrs[(p_slot + 1) % 3];
            let c2 = nbrs[(p_slot + 2) % 3];
            let d1 = cache.dir_index(c1, tree.neighbor_slot(c1, u));
            let d2 = cache.dir_index(c2, tree.neighbor_slot(c2, u));
            cache.costs[out_dir] = combine_dir_pair(
                &mut cache.states,
                out_dir,
                d1,
                d2,
                cache.costs[d1],
                cache.costs[d2],
                n_pat,
                pat_weights,
            );
        }
    }

    for idx in 0..cache.order.len() {
        let (u, _) = cache.order[idx];
        if tree.is_leaf(u) {
            continue;
        }
        let nbrs = tree.internal_nbrs[u - tree.n_taxa];
        for slot in 0..3 {
            let child = nbrs[slot];
            if idx > 0 && child == cache.order[idx].1 {
                continue;
            }
            let o1 = nbrs[(slot + 1) % 3];
            let o2 = nbrs[(slot + 2) % 3];
            let d1 = cache.dir_index(o1, tree.neighbor_slot(o1, u));
            let d2 = cache.dir_index(o2, tree.neighbor_slot(o2, u));
            let out_dir = cache.dir_index(u, slot);
            cache.costs[out_dir] = combine_dir_pair(
                &mut cache.states,
                out_dir,
                d1,
                d2,
                cache.costs[d1],
                cache.costs[d2],
                n_pat,
                pat_weights,
            );
        }
    }

    let first_child = tree.internal_nbrs[root - tree.n_taxa][0];
    let d_root = cache.dir_index(root, 0);
    let d_child = cache.dir_index(first_child, tree.neighbor_slot(first_child, root));
    let score = cache.edge_total_score(d_root, d_child, pat_weights);
    (score, cache)
}

#[inline]
fn other_two(nbrs: [usize; 3], excluded: usize) -> [usize; 2] {
    if nbrs[0] == excluded {
        [nbrs[1], nbrs[2]]
    } else if nbrs[1] == excluded {
        [nbrs[0], nbrs[2]]
    } else {
        debug_assert_eq!(nbrs[2], excluded);
        [nbrs[0], nbrs[1]]
    }
}

/// Exact Branch-and-Bound search initialized with a heuristic upper bound.
fn exact_branch_and_bound(
    aln: &CompressedAlignment,
    opts: &ParsimonyOptions,
) -> (UnrootedTree, usize) {
    let heuristic_opts = ParsimonyOptions {
        search: ParsimonySearch::Heuristic,
        swap: BranchSwap::Tbr,
        replicates: opts.replicates.max(4),
        gap_mode: opts.gap_mode,
        branch_lengths: ParsimonyBranchLengths::None,
        seed: opts.seed,
    };
    let (mut best_tree, mut best_score) = heuristic_search(aln, &heuristic_opts);

    if aln.n_inf == 0 {
        return (best_tree, 0);
    }

    // Order taxa by descending pairwise incompatibility / informativeness so
    // the partial score grows rapidly near the root of the search tree.
    let order = branch_and_bound_taxon_order(aln);

    // Precompute admissible suffix lower bounds for each depth k in 0..=n_taxa.
    let suffix_lb = compute_suffix_lower_bounds(aln, &order);

    let mut work_tree = UnrootedTree::star(aln.ntax, order[0], order[1], order[2]);
    let mut cache = FitchCache::new(aln.ntax, aln.n_inf);
    let initial_score = cache.compute(&work_tree, &aln.inf_states, &aln.inf_weights);

    bnb_recurse(
        aln,
        &order,
        &suffix_lb,
        3,
        initial_score,
        &mut work_tree,
        &mut best_score,
        &mut best_tree,
    );

    (best_tree, best_score)
}

fn branch_and_bound_taxon_order(aln: &CompressedAlignment) -> Vec<usize> {
    let n = aln.ntax;
    // Pick initial pair with maximum weighted mismatch, then greedily add max-distance taxon.
    let mut best_pair = (0, 1);
    let mut max_diff = 0usize;
    for i in 0..n {
        let si = &aln.inf_states[i * aln.n_inf..(i + 1) * aln.n_inf];
        for j in (i + 1)..n {
            let sj = &aln.inf_states[j * aln.n_inf..(j + 1) * aln.n_inf];
            let d = fitch_diff_cost(si, sj, &aln.inf_weights);
            if d > max_diff {
                max_diff = d;
                best_pair = (i, j);
            }
        }
    }

    let mut chosen = vec![false; n];
    let mut order = Vec::with_capacity(n);
    order.push(best_pair.0);
    order.push(best_pair.1);
    chosen[best_pair.0] = true;
    chosen[best_pair.1] = true;

    while order.len() < n {
        let mut next_t = 0usize;
        let mut best_min_d = 0usize;
        let mut found = false;
        for (t, &is_chosen) in chosen.iter().enumerate() {
            if is_chosen {
                continue;
            }
            let st = &aln.inf_states[t * aln.n_inf..(t + 1) * aln.n_inf];
            let min_d = order
                .iter()
                .map(|&prev| {
                    let sp = &aln.inf_states[prev * aln.n_inf..(prev + 1) * aln.n_inf];
                    fitch_diff_cost(st, sp, &aln.inf_weights)
                })
                .min()
                .unwrap_or(0);
            if !found || min_d > best_min_d {
                best_min_d = min_d;
                next_t = t;
                found = true;
            }
        }
        chosen[next_t] = true;
        order.push(next_t);
    }
    order
}

fn compute_suffix_lower_bounds(aln: &CompressedAlignment, order: &[usize]) -> Vec<usize> {
    let n = order.len();
    let mut lb = vec![0usize; n + 1];
    for depth in 3..n {
        let mut bound = 0usize;
        for p in 0..aln.n_inf {
            let mut placed_union = 0u32;
            for &t in &order[..depth] {
                placed_union |= aln.inf_states[t * aln.n_inf + p];
            }
            let mut extra_single_bits = 0u32;
            for &t in &order[depth..] {
                let s = aln.inf_states[t * aln.n_inf + p];
                if s.count_ones() == 1 && (s & placed_union) == 0 {
                    extra_single_bits |= s;
                }
            }
            bound += (extra_single_bits.count_ones() as usize) * aln.inf_weights[p] as usize;
        }
        lb[depth] = bound;
    }
    lb
}

#[allow(clippy::too_many_arguments)]
fn bnb_recurse(
    aln: &CompressedAlignment,
    order: &[usize],
    suffix_lb: &[usize],
    depth: usize,
    current_score: usize,
    tree: &mut UnrootedTree,
    best_score: &mut usize,
    best_tree: &mut UnrootedTree,
) {
    if current_score + suffix_lb[depth] >= *best_score {
        return;
    }
    if depth == order.len() {
        *best_score = current_score;
        *best_tree = tree.clone();
        return;
    }

    let taxon = order[depth];
    let leaf_states = &aln.inf_states[taxon * aln.n_inf..(taxon + 1) * aln.n_inf];

    let mut cache = FitchCache::new(aln.ntax, aln.n_inf);
    let _ = cache.compute(tree, &aln.inf_states, &aln.inf_weights);

    let mut edges = Vec::with_capacity(2 * depth);
    tree.collect_edges(&mut edges);

    let mut candidates: Vec<(usize, (usize, usize))> = edges
        .into_iter()
        .map(|(u, v)| {
            let extra = cache.eval_insert_extra(tree, u, v, leaf_states, &aln.inf_weights);
            (extra, (u, v))
        })
        .collect();
    candidates.sort_by_key(|(extra, _)| *extra);

    for (extra, (u, v)) in candidates {
        let next_score = current_score + extra;
        if next_score + suffix_lb[depth + 1] >= *best_score {
            break;
        }
        let new_int = tree.insert_on_edge(u, v, taxon);
        bnb_recurse(
            aln,
            order,
            suffix_lb,
            depth + 1,
            next_score,
            tree,
            best_score,
            best_tree,
        );
        tree.undo_insert(u, v, taxon, new_int);
    }
}

/// Convert `UnrootedTree` into `Tree`, optionally computing ACCTRAN/DELTRAN edge lengths.
fn unrooted_to_tree(
    unrooted: &UnrootedTree,
    aln: &CompressedAlignment,
    bl_mode: ParsimonyBranchLengths,
) -> Tree {
    let n = unrooted.n_taxa;
    let total_nodes = 2 * n - 2;
    let label_cap = aln.taxa.iter().map(String::len).sum::<usize>();
    let mut tree = Tree::with_capacity(total_nodes, label_cap);

    let mut mapped_nodes: Vec<Option<NodeId>> = vec![None; total_nodes];

    let root_u = n;
    let root_id = tree.add_node();
    tree.set_root(root_id);
    mapped_nodes[root_u] = Some(root_id);

    let edge_lengths = match bl_mode {
        ParsimonyBranchLengths::None => None,
        ParsimonyBranchLengths::Acctran | ParsimonyBranchLengths::Deltran => {
            Some(compute_mpr_edge_lengths(unrooted, aln, bl_mode))
        }
    };

    // Build rooted LCRS tree from `root_u`.
    let mut stack = vec![(root_u, usize::MAX)];
    while let Some((u, parent)) = stack.pop() {
        let parent_id = mapped_nodes[u].expect("Mapped parent node");
        let mut prev_child: Option<NodeId> = None;
        for &v in unrooted.neighbors(u) {
            if v == parent {
                continue;
            }
            let child_id = tree.add_node();
            mapped_nodes[v] = Some(child_id);
            if unrooted.is_leaf(v) {
                tree.set_label(child_id, &aln.taxa[v], false);
            }
            if let Some(ref lens) = edge_lengths {
                let len = lens.get(&(u.min(v), u.max(v))).copied().unwrap_or(0.0);
                tree.set_branch_length(child_id, len);
            }
            tree.attach_child(parent_id, child_id, prev_child);
            prev_child = Some(child_id);
            if !unrooted.is_leaf(v) {
                stack.push((v, u));
            }
        }
    }

    tree.metadata.rootedness = Rootedness::Unrooted;
    tree
}

/// Swofford & Maddison (1987) MPR two-pass ACCTRAN / DELTRAN branch length calculation.
fn compute_mpr_edge_lengths(
    unrooted: &UnrootedTree,
    aln: &CompressedAlignment,
    mode: ParsimonyBranchLengths,
) -> AHashMap<(usize, usize), f64> {
    let n_nodes = 2 * unrooted.n_taxa - 2;
    let n_pat = aln.n_all;

    // Pick an arbitrary edge `(r_a, r_b)` and temporarily root at virtual node `r_a` -> `r_b`
    // using FitchCache's already-defined directed semantics on a virtual binary root!
    let r_a = unrooted.n_taxa;
    let r_b = unrooted.internal_nbrs[0][0];

    // Collect children in rooted view where virtual root has children `(r_a, r_b)`.
    let mut parent_of = vec![usize::MAX; n_nodes];
    let mut children_of = vec![Vec::<usize>::new(); n_nodes];
    let mut postorder = Vec::with_capacity(n_nodes);

    let mut stack = vec![(r_a, r_b, false), (r_b, r_a, false)];
    parent_of[r_a] = usize::MAX;
    parent_of[r_b] = usize::MAX;

    while let Some((u, avoid, visited)) = stack.pop() {
        if visited {
            postorder.push(u);
            continue;
        }
        stack.push((u, avoid, true));
        for &v in unrooted.neighbors(u) {
            if v != avoid {
                parent_of[v] = u;
                children_of[u].push(v);
                stack.push((v, u, false));
            }
        }
    }

    // 1. Downpass sets D[u][p]
    let mut down = vec![0u32; n_nodes * n_pat];
    for &u in &postorder {
        if unrooted.is_leaf(u) {
            let src = &aln.all_states[u * n_pat..(u + 1) * n_pat];
            down[u * n_pat..(u + 1) * n_pat].copy_from_slice(src);
        } else {
            let c1 = children_of[u][0];
            let c2 = children_of[u][1];
            for p in 0..n_pat {
                let s1 = down[c1 * n_pat + p];
                let s2 = down[c2 * n_pat + p];
                let inter = s1 & s2;
                down[u * n_pat + p] = if inter != 0 { inter } else { s1 | s2 };
            }
        }
    }

    // Virtual root downpass set
    let mut vroot_down = vec![0u32; n_pat];
    for p in 0..n_pat {
        let sa = down[r_a * n_pat + p];
        let sb = down[r_b * n_pat + p];
        let inter = sa & sb;
        vroot_down[p] = if inter != 0 { inter } else { sa | sb };
    }

    // 2. Uppass & final MPR sets F[u][p] (only needed for DELTRAN, cheap to always compute)
    let mut up = vec![0u32; n_nodes * n_pat];
    let mut final_mpr = vec![0u32; n_nodes * n_pat];

    // Seed r_a (outgroup side is r_b) and r_b (outgroup side is r_a)
    up[r_a * n_pat..(r_a + 1) * n_pat].copy_from_slice(&down[r_b * n_pat..(r_b + 1) * n_pat]);
    up[r_b * n_pat..(r_b + 1) * n_pat].copy_from_slice(&down[r_a * n_pat..(r_a + 1) * n_pat]);

    for &u in postorder.iter().rev() {
        if unrooted.is_leaf(u) {
            for p in 0..n_pat {
                let d_u = down[u * n_pat + p];
                let u_up = up[u * n_pat + p];
                let inter = d_u & u_up;
                final_mpr[u * n_pat + p] = if inter != 0 { inter } else { d_u };
            }
        } else {
            let c1 = children_of[u][0];
            let c2 = children_of[u][1];
            for p in 0..n_pat {
                let u_up = up[u * n_pat + p];
                let d1 = down[c1 * n_pat + p];
                let d2 = down[c2 * n_pat + p];

                // Swofford & Maddison (1987) MPR final set from 3 neighbors (d1, d2, u_up)
                let maj = (d1 & d2) | (d1 & u_up) | (d2 & u_up);
                final_mpr[u * n_pat + p] = if maj != 0 { maj } else { d1 | d2 | u_up };

                let inter1 = u_up & d2;
                up[c1 * n_pat + p] = if inter1 != 0 { inter1 } else { u_up | d2 };
                let inter2 = u_up & d1;
                up[c2 * n_pat + p] = if inter2 != 0 { inter2 } else { u_up | d1 };
            }
        }
    }

    // 3. Top-down single-state assignment
    let candidate_sets = match mode {
        ParsimonyBranchLengths::Acctran => &down,
        ParsimonyBranchLengths::Deltran => &final_mpr,
        ParsimonyBranchLengths::None => unreachable!(),
    };

    let mut assigned = vec![0u32; n_nodes * n_pat];
    let mut root_state = vec![0u32; n_pat];
    for p in 0..n_pat {
        let mask = vroot_down[p];
        root_state[p] = 1u32 << mask.trailing_zeros();
    }

    let mut bridge_changes = 0usize;
    for p in 0..n_pat {
        let w = aln.all_weights[p] as usize;
        let r_s = root_state[p];
        let cand_a = candidate_sets[r_a * n_pat + p];
        let s_a = if (cand_a & r_s) != 0 {
            r_s
        } else {
            1u32 << cand_a.trailing_zeros()
        };
        assigned[r_a * n_pat + p] = s_a;

        let cand_b = candidate_sets[r_b * n_pat + p];
        let s_b = if (cand_b & r_s) != 0 {
            r_s
        } else {
            1u32 << cand_b.trailing_zeros()
        };
        assigned[r_b * n_pat + p] = s_b;

        if s_a != s_b {
            bridge_changes += w;
        }
    }

    let mut edge_changes: AHashMap<(usize, usize), usize> = AHashMap::new();
    edge_changes.insert((r_a.min(r_b), r_a.max(r_b)), bridge_changes);

    for &u in postorder.iter().rev() {
        if u == r_a || u == r_b {
            continue;
        }
        let parent = parent_of[u];
        let mut changes = 0usize;
        for p in 0..n_pat {
            let ps = assigned[parent * n_pat + p];
            let cand = candidate_sets[u * n_pat + p];
            let s = if (cand & ps) != 0 {
                ps
            } else {
                1u32 << cand.trailing_zeros()
            };
            assigned[u * n_pat + p] = s;
            if s != ps {
                changes += aln.all_weights[p] as usize;
            }
        }
        edge_changes.insert((u.min(parent), u.max(parent)), changes);
    }

    let denom = aln.nchar.max(1) as f64;
    edge_changes
        .into_iter()
        .map(|(k, count)| (k, count as f64 / denom))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn enumerate_all_unrooted_trees(n: usize) -> Vec<UnrootedTree> {
        let mut trees = vec![UnrootedTree::star(n, 0, 1, 2)];
        let mut edges = Vec::new();
        for taxon in 3..n {
            let mut next = Vec::new();
            for t in trees {
                t.collect_edges(&mut edges);
                for &(u, v) in &edges {
                    let mut clone = t.clone();
                    clone.insert_on_edge(u, v, taxon);
                    next.push(clone);
                }
            }
            trees = next;
        }
        trees
    }

    #[test]
    fn test_exact_bnb_matches_exhaustive_5_and_6_taxa() {
        // 5-taxon alignment with ambiguity and gaps
        let mut m5 = SeqMatrix::new();
        m5.insert("T1".into(), b"AAAAACCCC-GGGG");
        m5.insert("T2".into(), b"AAAATCCCC-GGGG");
        m5.insert("T3".into(), b"CCGGATTTTRAAAA");
        m5.insert("T4".into(), b"CCGGATTTTTAAAT");
        m5.insert("T5".into(), b"TTTTTTTTTCAAAA");

        for gap_mode in [GapMode::Missing, GapMode::FifthState] {
            let aln = CompressedAlignment::from_seq_matrix(&m5, &DataType::Dna, gap_mode).unwrap();
            let all_15 = enumerate_all_unrooted_trees(5);
            assert_eq!(all_15.len(), 15);

            let mut best_brute = usize::MAX;
            for ut in &all_15 {
                let mut cache = FitchCache::new(5, aln.n_inf);
                let s = aln.const_score + cache.compute(ut, &aln.inf_states, &aln.inf_weights);
                // Also confirm compressed score == direct uncompressed Tree score!
                let tree_obj = unrooted_to_tree(ut, &aln, ParsimonyBranchLengths::None);
                assert_eq!(s, aln.score_tree_direct(&tree_obj));
                best_brute = best_brute.min(s);
            }

            let opts = ParsimonyOptions {
                search: ParsimonySearch::Exact,
                gap_mode,
                ..Default::default()
            };
            let res = estimate_parsimony_from_compressed(&aln, &opts, true).unwrap();
            assert_eq!(res.score, best_brute);
            assert_eq!(res.score, aln.score_tree_direct(&res.tree));
        }

        // 6-taxon alignment (105 unrooted trees)
        let mut m6 = SeqMatrix::new();
        m6.insert("T1".into(), b"AAACCCGGGTTTAAA");
        m6.insert("T2".into(), b"AAACCCGGGTTTCCC");
        m6.insert("T3".into(), b"TTTAAAGGGAAACCC");
        m6.insert("T4".into(), b"TTTAAAGGGAAAGGG");
        m6.insert("T5".into(), b"CCCGGGTTTCCCAAA");
        m6.insert("T6".into(), b"CCCGGGTTTCCCGGG");

        let aln6 =
            CompressedAlignment::from_seq_matrix(&m6, &DataType::Dna, GapMode::Missing).unwrap();
        let all_105 = enumerate_all_unrooted_trees(6);
        assert_eq!(all_105.len(), 105);
        let best_brute6 = all_105
            .iter()
            .map(|ut| {
                let mut cache = FitchCache::new(6, aln6.n_inf);
                aln6.const_score + cache.compute(ut, &aln6.inf_states, &aln6.inf_weights)
            })
            .min()
            .unwrap();

        let res6 = estimate_parsimony(
            &m6,
            &DataType::Dna,
            &ParsimonyOptions {
                search: ParsimonySearch::Exact,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(res6.score, best_brute6);
    }

    #[test]
    fn test_acctran_and_deltran_preserve_total_steps_and_differ_on_ambiguous_case() {
        // Canonical 4-taxon case with ambiguous internal reconstruction:
        // ((A:A, B:C), (C:A, D:T)) on site 0: ACAT -> min steps = 2.
        // ACCTRAN pushes change toward root (internal edges), DELTRAN delays toward tips.
        let mut m = SeqMatrix::new();
        m.insert("A".into(), b"AA");
        m.insert("B".into(), b"CA");
        m.insert("C".into(), b"AT");
        m.insert("D".into(), b"TT");

        let acctran = estimate_parsimony(
            &m,
            &DataType::Dna,
            &ParsimonyOptions {
                search: ParsimonySearch::Exact,
                branch_lengths: ParsimonyBranchLengths::Acctran,
                ..Default::default()
            },
        )
        .unwrap();
        let deltran = estimate_parsimony(
            &m,
            &DataType::Dna,
            &ParsimonyOptions {
                search: ParsimonySearch::Exact,
                branch_lengths: ParsimonyBranchLengths::Deltran,
                ..Default::default()
            },
        )
        .unwrap();
        let none = estimate_parsimony(
            &m,
            &DataType::Dna,
            &ParsimonyOptions {
                search: ParsimonySearch::Exact,
                branch_lengths: ParsimonyBranchLengths::None,
                ..Default::default()
            },
        )
        .unwrap();

        let sum_acctran: f64 = acctran
            .tree
            .preorder()
            .filter_map(|id| acctran.tree.branch_length(id))
            .sum();
        let sum_deltran: f64 = deltran
            .tree
            .preorder()
            .filter_map(|id| deltran.tree.branch_length(id))
            .sum();

        // Total steps * nchar == parsimony score
        assert!(((sum_acctran * 2.0) - acctran.score as f64).abs() < 1e-9);
        assert!(((sum_deltran * 2.0) - deltran.score as f64).abs() < 1e-9);
        assert!(
            none.tree
                .preorder()
                .all(|id| none.tree.branch_length(id).is_none())
        );

        // Exhaustively verify ACCTRAN and DELTRAN total-step invariant over all 15 5-taxon
        // unrooted trees including A/C/G/T, ambiguities, and gaps.
        let mut m5 = SeqMatrix::new();
        m5.insert("T1".into(), b"AAAAACCCC-GGGG");
        m5.insert("T2".into(), b"TTTTTCCCC-GGGG");
        m5.insert("T3".into(), b"CCGGATTTTRAAAA");
        m5.insert("T4".into(), b"CCGGATTTTTAAAT");
        m5.insert("T5".into(), b"GGGGTTTTTCAAAA");

        let aln5 =
            CompressedAlignment::from_seq_matrix(&m5, &DataType::Dna, GapMode::Missing).unwrap();
        let all_15 = enumerate_all_unrooted_trees(5);
        for ut in &all_15 {
            let mut cache = FitchCache::new(5, aln5.n_inf);
            let exact_score =
                aln5.const_score + cache.compute(ut, &aln5.inf_states, &aln5.inf_weights);
            for mode in [
                ParsimonyBranchLengths::Acctran,
                ParsimonyBranchLengths::Deltran,
            ] {
                let tree = unrooted_to_tree(ut, &aln5, mode);
                let total_bl: f64 = tree
                    .preorder()
                    .filter_map(|id| tree.branch_length(id))
                    .sum();
                assert!(
                    ((total_bl * aln5.nchar as f64) - exact_score as f64).abs() < 1e-9,
                    "Mode {:?} branch sum {} * {} != exact score {}",
                    mode,
                    total_bl,
                    aln5.nchar,
                    exact_score
                );
            }
        }
    }

    #[test]
    fn test_swapping_monotonicity_and_aa_parsimony() {
        let mut m = SeqMatrix::new();
        m.insert("S1".into(), b"ARNDCQEGHI");
        m.insert("S2".into(), b"ARNDCQEGHL");
        m.insert("S3".into(), b"WWWWWWGGGG");
        m.insert("S4".into(), b"WWWWWWGGGA");
        m.insert("S5".into(), b"YYYYYYAAAA");

        let aln =
            CompressedAlignment::from_seq_matrix(&m, &DataType::Aa, GapMode::Missing).unwrap();
        let order = [0, 1, 2, 3, 4];
        let (mut tree_nni, s_add) = stepwise_addition(&aln, &order);
        let mut tree_spr = tree_nni.clone();
        let mut tree_tbr = tree_nni.clone();

        let s_nni = nni_search(&mut tree_nni, &aln);
        let s_spr = spr_search(&mut tree_spr, &aln);
        let s_tbr = tbr_search(&mut tree_tbr, &aln);

        assert!(s_nni <= s_add);
        assert!(s_spr <= s_add);
        assert!(s_tbr <= s_spr);
    }
}
