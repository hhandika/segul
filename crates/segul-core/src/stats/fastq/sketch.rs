//! Bounded-memory sketches for duplication and overrepresented sequences.
//!
//! Both sketches see every read of a file and merge across threads, so
//! large files are summarized from all their reads, not only the first ones.
use std::collections::BinaryHeap;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Mutex, MutexGuard};

use ahash::AHashMap;
use xxhash_rust::xxh3::xxh3_64;

/// Distinct sequences kept to estimate duplication, unless set otherwise.
pub const DEFAULT_DUP_SKETCH_SIZE: usize = 100_000;

/// Counters for overrepresented sequences: four per 1/0.001, so every
/// sequence above 0.1% of reads is found with counts within 0.025%.
pub(crate) const HEAVY_HITTERS: usize = 4_000;
// Long reads are keyed on their start, following FastQC.
const KEY_FULL_MAX: usize = 75;
const KEY_PREFIX: usize = 50;

// Kept sequences per shard when the sketch is split. Small sketches use
// one shard, so a sketch of k keeps exactly k sequences.
const MIN_SHARD_CAPACITY: usize = 1_024;
const MAX_SHARDS: usize = 64;

/// Bottom-k sketch shared by every worker of a file: the `k` distinct
/// sequences with the smallest hashes, each with its count over every read.
///
/// The kept sequences are a uniform sample of the distinct sequences, so
/// their counts estimate the duplication levels of the whole file. Hashes
/// are split into shards by range, each keeping the smallest hashes in its
/// range. A read above its shard's threshold is skipped without a lock,
/// which is almost every read once the sketch is full. The kept set and its
/// counts are the same in whatever order reads arrive, because a kept hash
/// never has enough smaller hashes in its shard to be evicted.
#[derive(Debug)]
pub(crate) struct DuplicationSketch {
    /// Distinct sequences kept. Zero keeps every one.
    capacity: usize,
    /// Sequences kept per shard, or zero for no limit.
    shard_capacity: usize,
    shards: Vec<Mutex<Shard>>,
    gates: Vec<ShardGate>,
}

/// Lock-free state of a shard, read for every read. Each sits on its own
/// cache line, so updates during warm-up do not slow reads of the others.
#[derive(Debug)]
#[repr(align(64))]
struct ShardGate {
    /// Largest kept hash once the shard is full. Larger hashes are skipped.
    threshold: AtomicU64,
    /// Whether a distinct sequence was left out, so the counts are a sample.
    saturated: AtomicBool,
}

impl ShardGate {
    fn new() -> Self {
        Self {
            threshold: AtomicU64::new(u64::MAX),
            saturated: AtomicBool::new(false),
        }
    }

    fn saturate(&self) {
        // Reading first avoids writing the shared line for every read.
        if !self.saturated.load(Ordering::Relaxed) {
            self.saturated.store(true, Ordering::Relaxed);
        }
    }
}

#[derive(Debug, Default)]
struct Shard {
    counts: AHashMap<u64, u32>,
    /// Kept hashes, largest on top, to find the one to evict.
    kept: BinaryHeap<u64>,
}

impl Shard {
    /// A shard for up to `capacity` hashes. Zero grows as needed.
    fn with_capacity(capacity: usize) -> Self {
        Self {
            counts: AHashMap::with_capacity(capacity),
            kept: BinaryHeap::with_capacity(capacity),
        }
    }
}

impl DuplicationSketch {
    pub(crate) fn new(capacity: usize) -> Self {
        let shard_count = if capacity == 0 {
            MAX_SHARDS
        } else {
            (capacity / MIN_SHARD_CAPACITY).clamp(1, MAX_SHARDS)
        };
        let shard_capacity = capacity.div_ceil(shard_count);
        Self {
            capacity,
            shard_capacity,
            // Sized once, so shards do not reallocate as they fill.
            shards: (0..shard_count)
                .map(|_| Mutex::new(Shard::with_capacity(shard_capacity)))
                .collect(),
            gates: (0..shard_count).map(|_| ShardGate::new()).collect(),
        }
    }

    pub(crate) fn add(&self, sequence: &[u8]) {
        self.add_hash(xxh3_64(sequence));
    }

    fn add_hash(&self, hash: u64) {
        let index = self.shard_of(hash);
        let gate = &self.gates[index];
        // A stale threshold is only ever larger, so this never skips a
        // hash that belongs in the sketch.
        if hash > gate.threshold.load(Ordering::Relaxed) {
            gate.saturate();
            return;
        }
        let mut shard = self.lock(index);
        if let Some(count) = shard.counts.get_mut(&hash) {
            *count += 1;
            return;
        }
        if self.shard_capacity == 0 {
            shard.counts.insert(hash, 1);
            return;
        }
        if shard.counts.len() == self.shard_capacity {
            let largest = *shard.kept.peek().expect("Full shard has hashes");
            gate.saturate();
            if hash > largest {
                return;
            }
            shard.kept.pop();
            shard.counts.remove(&largest);
        }
        shard.counts.insert(hash, 1);
        shard.kept.push(hash);
        if shard.counts.len() == self.shard_capacity {
            let largest = *shard.kept.peek().expect("Full shard has hashes");
            gate.threshold.store(largest, Ordering::Relaxed);
        }
    }

    pub(crate) fn capacity(&self) -> usize {
        self.capacity
    }

    /// Whether every distinct sequence was counted.
    pub(crate) fn is_exact(&self) -> bool {
        self.gates
            .iter()
            .all(|g| !g.saturated.load(Ordering::Relaxed))
    }

    /// Distinct sequences in the file: exact unless the sketch was full,
    /// then estimated from the largest kept hash of each shard, with
    /// about 1/√k error.
    pub(crate) fn distinct(&self) -> f64 {
        let width = 2f64.powi(64) / self.shards.len() as f64;
        (0..self.shards.len())
            .map(|i| {
                let shard = self.lock(i);
                match shard.kept.peek() {
                    Some(&largest) if self.gates[i].saturated.load(Ordering::Relaxed) => {
                        let start = i as f64 * width;
                        let fraction = (largest as f64 - start + 1.0) / width;
                        (shard.counts.len() as f64 - 1.0) / fraction
                    }
                    _ => shard.counts.len() as f64,
                }
            })
            .sum()
    }

    /// Copies of each kept sequence.
    pub(crate) fn copies(&self) -> Vec<u64> {
        (0..self.shards.len())
            .flat_map(|i| {
                self.lock(i)
                    .counts
                    .values()
                    .map(|c| u64::from(*c))
                    .collect::<Vec<u64>>()
            })
            .collect()
    }

    fn shard_of(&self, hash: u64) -> usize {
        ((u128::from(hash) * self.shards.len() as u128) >> 64) as usize
    }

    fn lock(&self, index: usize) -> MutexGuard<'_, Shard> {
        self.shards[index]
            .lock()
            .expect("Failed locking the duplication sketch")
    }
}

/// Sequence start counted by [`HeavyHitters`], stored inline.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Counter {
    pub(crate) count: u64,
    key: [u8; KEY_FULL_MAX],
    len: u8,
}

impl Counter {
    pub(crate) fn sequence(&self) -> &[u8] {
        &self.key[..self.len as usize]
    }
}

/// Misra-Gries summary of the most common read starts.
///
/// Every sequence above 1/(k + 1) of the reads keeps a counter. Counts are
/// lower bounds, short by at most reads/(k + 1), and exact while the
/// summary never fills.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct HeavyHitters {
    capacity: usize,
    counters: AHashMap<u64, Counter>,
    reads: u64,
    saturated: bool,
}

impl HeavyHitters {
    pub(crate) fn new(capacity: usize) -> Self {
        Self {
            capacity,
            counters: AHashMap::new(),
            reads: 0,
            saturated: false,
        }
    }

    pub(crate) fn add(&mut self, sequence: &[u8]) {
        let key = if sequence.len() > KEY_FULL_MAX {
            &sequence[..KEY_PREFIX]
        } else {
            sequence
        };
        self.reads += 1;
        let hash = xxh3_64(key);
        if let Some(counter) = self.counters.get_mut(&hash) {
            counter.count += 1;
        } else if self.counters.len() < self.capacity {
            let mut stored = [0u8; KEY_FULL_MAX];
            stored[..key.len()].copy_from_slice(key);
            self.counters.insert(
                hash,
                Counter {
                    count: 1,
                    key: stored,
                    len: key.len() as u8,
                },
            );
        } else {
            // Each decrement removes k + 1 counts, so it runs at most
            // reads/(k + 1) times.
            self.counters.retain(|_, c| {
                c.count -= 1;
                c.count > 0
            });
            self.saturated = true;
        }
    }

    pub(crate) fn merge(&mut self, mut other: Self) {
        if other.counters.len() > self.counters.len() {
            std::mem::swap(&mut self.counters, &mut other.counters);
        }
        other.counters.into_iter().for_each(|(hash, counter)| {
            self.counters
                .entry(hash)
                .and_modify(|c| c.count += counter.count)
                .or_insert(counter);
        });
        self.reads += other.reads;
        self.saturated |= other.saturated;
        if self.counters.len() > self.capacity {
            let mut counts: Vec<u64> = self.counters.values().map(|c| c.count).collect();
            counts.sort_unstable_by(|a, b| b.cmp(a));
            let cut = counts[self.capacity];
            self.counters.retain(|_, c| {
                c.count = c.count.saturating_sub(cut);
                c.count > 0
            });
            self.saturated = true;
        }
    }

    pub(crate) fn reads(&self) -> u64 {
        self.reads
    }

    /// Largest amount a count can be short by.
    pub(crate) fn max_error(&self) -> u64 {
        if self.saturated {
            self.reads / (self.capacity as u64 + 1)
        } else {
            0
        }
    }

    pub(crate) fn counters(&self) -> impl Iterator<Item = &Counter> {
        self.counters.values()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    /// Deterministic pseudo-random hashes.
    fn hashes(n: u64, seed: u64) -> impl Iterator<Item = u64> {
        (0..n).map(move |i| xxh3_64(&(i ^ seed).to_le_bytes()))
    }

    /// Kept hashes and counts, sorted, to compare sketches.
    fn contents(sketch: &DuplicationSketch) -> Vec<(u64, u32)> {
        let mut all: Vec<(u64, u32)> = (0..sketch.shards.len())
            .flat_map(|i| {
                let shard = sketch.lock(i);
                shard
                    .counts
                    .iter()
                    .map(|(h, c)| (*h, *c))
                    .collect::<Vec<_>>()
            })
            .collect();
        all.sort();
        all
    }

    #[test]
    fn test_sketch_exact_below_capacity() {
        let sketch = DuplicationSketch::new(10);
        for seq in [&b"AC"[..], b"AC", b"GT"] {
            sketch.add(seq);
        }
        assert!(sketch.is_exact());
        assert_eq!(sketch.distinct(), 2.0);
        let mut copies = sketch.copies();
        copies.sort();
        assert_eq!(copies, vec![1, 2]);
        let unbounded = DuplicationSketch::new(0);
        hashes(1_000, 1).for_each(|h| unbounded.add_hash(h));
        assert!(unbounded.is_exact());
        assert_eq!(unbounded.distinct(), 1_000.0);
    }

    #[test]
    fn test_sketch_estimates_distinct() {
        let sketch = DuplicationSketch::new(10_000);
        // 500,000 distinct sequences, each seen twice.
        hashes(500_000, 7)
            .chain(hashes(500_000, 7))
            .for_each(|h| sketch.add_hash(h));
        assert!(!sketch.is_exact());
        // Nine shards of 1,112 round the capacity up slightly.
        assert_eq!(sketch.copies().len(), 9 * 1_112);
        let distinct = sketch.distinct();
        assert!((distinct / 500_000.0 - 1.0).abs() < 0.03, "{}", distinct);
        // Kept sequences have exact counts.
        assert!(sketch.copies().iter().all(|c| *c == 2));
    }

    #[test]
    fn test_sketch_same_in_any_order() {
        let all: Vec<u64> = hashes(40_000, 3).chain(hashes(10_000, 3)).collect();
        let sequential = DuplicationSketch::new(2_048);
        all.iter().for_each(|h| sequential.add_hash(*h));
        // Reversed and split across threads, as blocks reach workers.
        let parallel = DuplicationSketch::new(2_048);
        std::thread::scope(|scope| {
            for part in all.rchunks(7_000) {
                let sketch = &parallel;
                scope.spawn(move || part.iter().rev().for_each(|h| sketch.add_hash(*h)));
            }
        });
        assert_eq!(contents(&parallel), contents(&sequential));
        assert_eq!(parallel.distinct(), sequential.distinct());
        assert!(!parallel.is_exact());
    }

    #[test]
    fn test_heavy_hitters_find_frequent() {
        let mut hitters = HeavyHitters::new(10);
        // One start in 20% of reads among 1,000 unique ones.
        for i in 0..1_000u32 {
            hitters.add(format!("ACGT{}", i).as_bytes());
            if i % 4 == 0 {
                hitters.add(b"GGGGCCCC");
            }
        }
        let frequent = hitters
            .counters()
            .find(|c| c.sequence() == b"GGGGCCCC")
            .unwrap();
        let error = hitters.max_error();
        assert!(frequent.count <= 250 && frequent.count + error >= 250);
        assert_eq!(error, 1_250 / 11);
    }

    #[test]
    fn test_heavy_hitters_merge() {
        let mut first = HeavyHitters::new(3);
        let mut second = HeavyHitters::new(3);
        ["AA", "AA", "CC", "GG"]
            .iter()
            .for_each(|s| first.add(s.as_bytes()));
        ["AA", "TT", "TT", "CA"]
            .iter()
            .for_each(|s| second.add(s.as_bytes()));
        first.merge(second);
        assert_eq!(first.reads(), 8);
        let aa = first.counters().find(|c| c.sequence() == b"AA").unwrap();
        // Six distinct starts for three counters: counts drop by at most 8/4.
        assert!(aa.count + first.max_error() >= 3);
        assert!(first.counters().count() <= 3);
        // Reads longer than 75 bp are keyed on their first 50 bp.
        let mut long = HeavyHitters::new(3);
        long.add(&[b'A'; 80]);
        long.add(&[b'A'; 90]);
        assert_eq!(long.counters().next().unwrap().count, 2);
        assert_eq!(long.counters().next().unwrap().sequence().len(), 50);
    }
}
