//! Cache-compact bitsets and 128-bit SplitMix64 Zobrist hashing for clades and splits.

use std::cmp::Ordering;

const ZOBRIST_SEED_1: u64 = 0x5E60_17EE_50AA_A811;
const ZOBRIST_SEED_2: u64 = 0x9E37_79B9_7F4A_7C15;

/// Compact 32-bit identifier for a unique clade or bipartition in [`CladeBitsArena`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) struct CladeId(pub(crate) u32);

impl CladeId {
    #[inline]
    pub(crate) fn index(self) -> usize {
        self.0 as usize
    }
}

/// Deterministic 128-bit Zobrist hash table for $O(1)$ internal node clade hashing.
#[derive(Debug, Clone)]
pub(crate) struct ZobristTable {
    h1: Vec<u64>,
    h2: Vec<u64>,
    all_h1: u64,
    all_h2: u64,
    n_taxa: usize,
    words: usize,
}

impl ZobristTable {
    /// Build deterministic 128-bit SplitMix64 keys for `n_taxa` taxa.
    pub(crate) fn new(n_taxa: usize) -> Self {
        let words = n_taxa.div_ceil(64).max(1);
        let mut h1 = Vec::with_capacity(n_taxa);
        let mut h2 = Vec::with_capacity(n_taxa);
        let mut all_h1 = 0u64;
        let mut all_h2 = 0u64;

        let mut s1 = ZOBRIST_SEED_1;
        let mut s2 = ZOBRIST_SEED_2;
        for _ in 0..n_taxa {
            let k1 = splitmix64(&mut s1);
            let k2 = splitmix64(&mut s2);
            h1.push(k1);
            h2.push(k2);
            all_h1 ^= k1;
            all_h2 ^= k2;
        }

        Self {
            h1,
            h2,
            all_h1,
            all_h2,
            n_taxa,
            words,
        }
    }

    #[inline]
    pub(crate) fn words(&self) -> usize {
        self.words
    }

    #[inline]
    pub(crate) fn taxon_hash(&self, taxon_idx: usize) -> (u64, u64) {
        (self.h1[taxon_idx], self.h2[taxon_idx])
    }

    #[inline]
    pub(crate) fn all_hash(&self) -> (u64, u64) {
        (self.all_h1, self.all_h2)
    }

    /// Compute the 128-bit Zobrist hash directly from a bitset slice.
    pub(crate) fn hash_bits(&self, bits: &[u64]) -> (u64, u64) {
        let mut x1 = 0u64;
        let mut x2 = 0u64;
        for (word_idx, &word) in bits.iter().enumerate() {
            let mut w = word;
            let base = word_idx * 64;
            while w != 0 {
                let tz = w.trailing_zeros() as usize;
                let idx = base + tz;
                if idx < self.n_taxa {
                    x1 ^= self.h1[idx];
                    x2 ^= self.h2[idx];
                }
                w &= w - 1;
            }
        }
        (x1, x2)
    }
}

#[inline]
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Flat contiguous arena storing fixed-width clade bitsets with zero per-clade allocations.
#[derive(Debug, Clone)]
pub(crate) struct CladeBitsArena {
    data: Vec<u64>,
    words: usize,
    count: usize,
}

impl CladeBitsArena {
    pub(crate) fn new(words: usize) -> Self {
        Self {
            data: Vec::new(),
            words: words.max(1),
            count: 0,
        }
    }

    /// Append `bits` to the arena and return its [`CladeId`].
    #[inline]
    pub(crate) fn push(&mut self, bits: &[u64]) -> CladeId {
        debug_assert_eq!(bits.len(), self.words);
        let id = u32::try_from(self.count).expect("Clade count exceeds u32::MAX");
        self.data.extend_from_slice(bits);
        self.count += 1;
        CladeId(id)
    }

    #[inline]
    pub(crate) fn get(&self, id: CladeId) -> &[u64] {
        let start = id.index() * self.words;
        &self.data[start..start + self.words]
    }
}

#[inline]
pub(crate) fn last_word_mask(n_taxa: usize) -> u64 {
    let rem = n_taxa % 64;
    if rem == 0 && n_taxa > 0 {
        u64::MAX
    } else if n_taxa == 0 {
        0
    } else {
        (1u64 << rem) - 1
    }
}

#[inline]
pub(crate) fn set_bit(bits: &mut [u64], idx: usize) {
    bits[idx / 64] |= 1u64 << (idx % 64);
}

#[inline]
pub(crate) fn bit_is_set(bits: &[u64], idx: usize) -> bool {
    (bits[idx / 64] & (1u64 << (idx % 64))) != 0
}

#[inline]
pub(crate) fn popcount(bits: &[u64]) -> usize {
    bits.iter().map(|w| w.count_ones() as usize).sum()
}

#[inline]
pub(crate) fn first_set_bit(bits: &[u64]) -> Option<usize> {
    for (word_idx, &w) in bits.iter().enumerate() {
        if w != 0 {
            return Some(word_idx * 64 + w.trailing_zeros() as usize);
        }
    }
    None
}

/// Canonicalize an unrooted bipartition so taxon `0` is always in the complement (`0` bit).
#[inline]
pub(crate) fn canonicalize_unrooted(
    bits: &mut [u64],
    hash: &mut (u64, u64),
    n_taxa: usize,
    zobrist: &ZobristTable,
) {
    if (bits[0] & 1) != 0 {
        for word in bits.iter_mut() {
            *word = !*word;
        }
        if let Some(last) = bits.last_mut() {
            *last &= last_word_mask(n_taxa);
        }
        let all = zobrist.all_hash();
        hash.0 ^= all.0;
        hash.1 ^= all.1;
    }
}

/// Canonicalize `bits` relative to `mask` for partial-taxa split comparison.
///
/// Intersects `bits` with `mask`, and when `!rooted`, flips relative to `mask` so the
/// lowest set bit of `mask` is `0`.
pub(crate) fn canonicalize_for_mask(bits: &mut [u64], mask: &[u64], rooted: bool) {
    for (b, &m) in bits.iter_mut().zip(mask.iter()) {
        *b &= m;
    }
    if !rooted
        && let Some(anchor) = first_set_bit(mask)
        && bit_is_set(bits, anchor)
    {
        for (b, &m) in bits.iter_mut().zip(mask.iter()) {
            *b = (!*b) & m;
        }
    }
}

/// Return true if two rooted clades are compatible (disjoint or one contains the other).
#[inline]
pub(crate) fn are_compatible_rooted(a: &[u64], b: &[u64]) -> bool {
    let mut has_ab = false;
    let mut has_a_not_b = false;
    let mut has_b_not_a = false;

    for (&wa, &wb) in a.iter().zip(b.iter()) {
        has_ab |= (wa & wb) != 0;
        has_a_not_b |= (wa & !wb) != 0;
        has_b_not_a |= (!wa & wb) != 0;
        if has_ab && has_a_not_b && has_b_not_a {
            return false;
        }
    }
    true
}

/// Return true if two unrooted bipartitions are compatible.
#[inline]
pub(crate) fn are_compatible_unrooted(a: &[u64], b: &[u64], n_taxa: usize) -> bool {
    let mut has_ab = false;
    let mut has_a_not_b = false;
    let mut has_b_not_a = false;
    let mut has_neither = false;
    let last_idx = a.len().saturating_sub(1);
    let tail_mask = last_word_mask(n_taxa);

    for (idx, (&wa, &wb)) in a.iter().zip(b.iter()).enumerate() {
        let valid = if idx == last_idx { tail_mask } else { u64::MAX };
        has_ab |= (wa & wb) != 0;
        has_a_not_b |= (wa & !wb) != 0;
        has_b_not_a |= (!wa & wb) != 0;
        has_neither |= ((!wa & !wb) & valid) != 0;
        if has_ab && has_a_not_b && has_b_not_a && has_neither {
            return false;
        }
    }
    true
}

#[inline]
pub(crate) fn is_subset(sub: &[u64], sup: &[u64]) -> bool {
    sub.iter().zip(sup.iter()).all(|(&s, &p)| (s & !p) == 0)
}

#[inline]
pub(crate) fn cmp_bits(a: &[u64], b: &[u64]) -> Ordering {
    for (&x, &y) in a.iter().zip(b.iter()) {
        match x.cmp(&y) {
            Ordering::Equal => continue,
            non_eq => return non_eq,
        }
    }
    Ordering::Equal
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zobrist_invariance_and_large_taxa_bitset() {
        let n_taxa = 300;
        let zobrist = ZobristTable::new(n_taxa);
        assert_eq!(zobrist.words(), 5);

        let mut bits = vec![0u64; zobrist.words()];
        let mut h = (0u64, 0u64);
        for tid in [0usize, 63, 64, 150, 299] {
            set_bit(&mut bits, tid);
            let th = zobrist.taxon_hash(tid);
            h.0 ^= th.0;
            h.1 ^= th.1;
        }
        assert_eq!(h, zobrist.hash_bits(&bits));

        // Canonicalize unrooted split (bit 0 was set -> flips to complement)
        let mut comp = vec![u64::MAX; zobrist.words()];
        *comp.last_mut().unwrap() &= last_word_mask(n_taxa);
        for tid in [0usize, 63, 64, 150, 299] {
            comp[tid / 64] &= !(1u64 << (tid % 64));
        }

        canonicalize_unrooted(&mut bits, &mut h, n_taxa, &zobrist);
        assert_eq!(bits, comp);
        assert_eq!(h, zobrist.hash_bits(&bits));
    }

    #[test]
    fn test_rooted_and_unrooted_compatibility() {
        // Rooted: {1,2} and {1,2,3} compatible; {1,2} and {2,3} incompatible
        let a = [0b0110u64];
        let b = [0b1110u64];
        let c = [0b1100u64];
        assert!(are_compatible_rooted(&a, &b));
        assert!(!are_compatible_rooted(&a, &c));

        // Unrooted over 4 taxa: {2,3} vs {0,1} and {1,3} vs {0,2} are incompatible
        let s1 = [0b1100u64];
        let s2 = [0b1010u64];
        assert!(!are_compatible_unrooted(&s1, &s2, 4));
        // {2,3} and {3,4} on 5 taxa share 3, differ on 2 vs 4, and {0,1} are in neither -> incompatible.
        let s3 = [0b10000u64 | 0b1000u64]; // {3,4}
        assert!(!are_compatible_unrooted(&[0b01100u64], &s3, 5));
        let s4 = [0b11100u64]; // {2,3,4}
        assert!(are_compatible_unrooted(&[0b01100u64], &s4, 5));
    }
}
