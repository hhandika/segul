//! Single-pass accumulator for FASTQ statistics.
//!
//! Every field is a count or a histogram, so accumulators built from
//! different blocks of a file merge by addition. The result does not depend
//! on how the file was split between threads.
use std::collections::BTreeMap;
use std::sync::LazyLock;

use anyhow::{Result, bail};
use memchr::memmem::Finder;

use super::parse::RecordParser;
use super::sketch::{DuplicationSketch, HEAVY_HITTERS, HeavyHitters};

/// Highest quality score tracked. Phred+33 ends at '~' (Q93).
pub(crate) const MAX_QSCORE: usize = 93;
pub(crate) const QSCORE_BINS: usize = MAX_QSCORE + 1;
/// Quality below this is low quality.
pub(crate) const LOW_QSCORE: usize = 20;
/// Positions tracked one by one. Later positions are grouped in bins.
pub(crate) const EXACT_POSITIONS: usize = 1_000;
pub(crate) const POSITION_BIN_WIDTH: usize = 1_000;
// Read lengths up to this are counted in a vector, longer ones in a map.
const MAX_DENSE_LENGTH: usize = 1 << 16;
/// Reads are split between neighboring GC bins in these units, so
/// counts stay integers and merge exactly.
pub(crate) const GC_WEIGHT: u64 = 1_000;
// Reads scanned to detect the quality offset.
const OFFSET_SAMPLE_READS: usize = 10_000;
const PHRED64_MIN: u8 = 64;

/// Base classes counted for every read.
pub(crate) const BASE_CLASSES: usize = 6;
pub(crate) const BASE_A: usize = 0;
pub(crate) const BASE_C: usize = 1;
pub(crate) const BASE_G: usize = 2;
pub(crate) const BASE_T: usize = 3;
pub(crate) const BASE_N: usize = 4;
const BASE_OTHER: usize = 5;

static BASE_INDEX: [u8; 256] = build_base_index();

const fn build_base_index() -> [u8; 256] {
    let mut table = [BASE_OTHER as u8; 256];
    table[b'A' as usize] = BASE_A as u8;
    table[b'a' as usize] = BASE_A as u8;
    table[b'C' as usize] = BASE_C as u8;
    table[b'c' as usize] = BASE_C as u8;
    table[b'G' as usize] = BASE_G as u8;
    table[b'g' as usize] = BASE_G as u8;
    table[b'T' as usize] = BASE_T as u8;
    table[b't' as usize] = BASE_T as u8;
    table[b'N' as usize] = BASE_N as u8;
    table[b'n' as usize] = BASE_N as u8;
    table
}

/// Adapters and artifacts searched in Complete mode, as (name, sequence).
pub(crate) const ADAPTERS: [(&str, &str); 5] = [
    ("Illumina Universal", "AGATCGGAAGAG"),
    ("Nextera", "CTGTCTCTTATA"),
    ("Illumina Small RNA", "TGGAATTCTCGG"),
    ("PolyA", "AAAAAAAAAAAA"),
    ("PolyG", "GGGGGGGGGGGG"),
];
/// Adapters counted as library adapters, not homopolymer artifacts.
pub(crate) const LIBRARY_ADAPTERS: usize = 3;
pub(crate) const POLY_G: usize = 4;

static ADAPTER_FINDERS: LazyLock<Vec<Finder<'static>>> = LazyLock::new(|| {
    ADAPTERS
        .iter()
        .map(|(_, seq)| Finder::new(seq.as_bytes()).into_owned())
        .collect()
});

/// Bin of a 0-based read position.
pub(crate) fn position_bin(pos: usize) -> usize {
    if pos < EXACT_POSITIONS {
        pos
    } else {
        EXACT_POSITIONS + (pos - EXACT_POSITIONS) / POSITION_BIN_WIDTH
    }
}

/// First and last 1-based read position in a bin.
pub(crate) fn bin_range(bin: usize) -> (usize, usize) {
    if bin < EXACT_POSITIONS {
        (bin + 1, bin + 1)
    } else {
        let start = EXACT_POSITIONS + (bin - EXACT_POSITIONS) * POSITION_BIN_WIDTH;
        (start + 1, start + POSITION_BIN_WIDTH)
    }
}

/// Detect the quality offset from the first reads of a block.
///
/// Like FastQC, any score below '@' means Phred+33, because Phred+64
/// data cannot contain it.
pub(crate) fn detect_offset(data: &[u8]) -> u8 {
    let lowest = RecordParser::new(data, 0)
        .take(OFFSET_SAMPLE_READS)
        .filter_map(|r| r.ok())
        .filter_map(|r| r.quality.iter().min().copied())
        .min();
    match lowest {
        Some(q) if q >= PHRED64_MIN => 64,
        _ => 33,
    }
}

/// Read length counts.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct LengthHistogram {
    dense: Vec<u64>,
    sparse: BTreeMap<usize, u64>,
}

impl LengthHistogram {
    fn add(&mut self, len: usize) {
        if len < MAX_DENSE_LENGTH {
            if len >= self.dense.len() {
                self.dense.resize(len + 1, 0);
            }
            self.dense[len] += 1;
        } else {
            *self.sparse.entry(len).or_insert(0) += 1;
        }
    }

    fn merge(&mut self, other: Self) {
        if other.dense.len() > self.dense.len() {
            self.dense.resize(other.dense.len(), 0);
        }
        self.dense
            .iter_mut()
            .zip(other.dense)
            .for_each(|(a, b)| *a += b);
        other
            .sparse
            .into_iter()
            .for_each(|(len, count)| *self.sparse.entry(len).or_insert(0) += count);
    }

    /// Lengths with at least one read, shortest first, as (length, reads).
    pub(crate) fn counts(&self) -> Vec<(usize, u64)> {
        self.dense
            .iter()
            .enumerate()
            .filter(|(_, c)| **c > 0)
            .map(|(len, c)| (len, *c))
            .chain(self.sparse.iter().map(|(len, c)| (*len, *c)))
            .collect()
    }
}

/// Base and quality counts at one read position or bin.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PositionCounts {
    pub(crate) bases: [u64; BASE_CLASSES],
    pub(crate) qualities: [u64; QSCORE_BINS],
}

impl PositionCounts {
    fn new() -> Self {
        Self {
            bases: [0; BASE_CLASSES],
            qualities: [0; QSCORE_BINS],
        }
    }

    fn merge(&mut self, other: &Self) {
        add_arrays(&mut self.bases, &other.bases);
        add_arrays(&mut self.qualities, &other.qualities);
    }
}

/// Counts collected from FASTQ records.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ReadAccumulator {
    complete: bool,
    offset: u8,
    pub(crate) reads: u64,
    pub(crate) bases: [u64; BASE_CLASSES],
    pub(crate) lengths: LengthHistogram,
    /// Base quality counts. In Complete mode they are summed from `positions`.
    pub(crate) qualities: [u64; QSCORE_BINS],
    /// Reads by their mean quality, rounded down.
    pub(crate) read_qualities: [u64; QSCORE_BINS],
    /// Reads by GC percentage of their A, C, G, and T bases, in units of
    /// [`GC_WEIGHT`]. See [`add_gc`] for how a read is split between bins.
    pub(crate) gc: [u64; 101],
    pub(crate) positions: Vec<PositionCounts>,
    /// Reads by the bin where each adapter first starts.
    pub(crate) adapters: Vec<Vec<u64>>,
    /// Reads with any library adapter.
    pub(crate) adapter_reads: u64,
    /// Most common read starts, for overrepresented sequences.
    pub(crate) heavy_hitters: HeavyHitters,
    invalid_quality: bool,
}

impl ReadAccumulator {
    /// `complete` adds per-position, adapter, and overrepresented counts.
    pub(crate) fn new(complete: bool, offset: u8) -> Self {
        Self {
            complete,
            offset,
            reads: 0,
            bases: [0; BASE_CLASSES],
            lengths: LengthHistogram::default(),
            qualities: [0; QSCORE_BINS],
            read_qualities: [0; QSCORE_BINS],
            gc: [0; 101],
            positions: Vec::new(),
            adapters: vec![Vec::new(); ADAPTERS.len()],
            adapter_reads: 0,
            heavy_hitters: HeavyHitters::new(HEAVY_HITTERS),
            invalid_quality: false,
        }
    }

    pub(crate) fn offset(&self) -> u8 {
        self.offset
    }

    /// Add every record in a block. `first_read` is the index of its first
    /// read. In Complete mode, reads are also added to `duplicates`, which
    /// every worker of the file shares.
    pub(crate) fn add_chunk(
        &mut self,
        data: &[u8],
        first_read: u64,
        duplicates: Option<&DuplicationSketch>,
    ) -> Result<()> {
        for record in RecordParser::new(data, first_read) {
            let record = record?;
            self.add_record(record.sequence, record.quality);
            if let Some(sketch) = duplicates
                && !record.sequence.is_empty()
            {
                sketch.add(record.sequence);
            }
        }
        if self.invalid_quality {
            bail!(
                "Found quality scores outside the Phred+{} range. \
                Try specifying the quality offset.",
                self.offset
            );
        }
        Ok(())
    }

    fn add_record(&mut self, sequence: &[u8], quality: &[u8]) {
        let len = sequence.len();
        self.reads += 1;
        self.lengths.add(len);
        if len == 0 {
            return;
        }
        let (counts, quality_sum) = if self.complete {
            self.count_positions(sequence, quality)
        } else {
            self.count_bases(sequence, quality)
        };
        add_arrays(&mut self.bases, &counts.map(u64::from));
        let mean_quality = (quality_sum / len as u64) as usize;
        self.read_qualities[mean_quality.min(MAX_QSCORE)] += 1;
        let acgt = counts[BASE_A] + counts[BASE_C] + counts[BASE_G] + counts[BASE_T];
        if acgt > 0 {
            add_gc(&mut self.gc, counts[BASE_G] + counts[BASE_C], acgt);
        }
        if self.complete {
            self.find_adapters(sequence);
            self.heavy_hitters.add(sequence);
        }
    }

    /// Count bases and qualities for the whole file only.
    fn count_bases(&mut self, sequence: &[u8], quality: &[u8]) -> ([u32; BASE_CLASSES], u64) {
        let mut counts = [0u32; BASE_CLASSES];
        let mut quality_sum = 0u64;
        let mut invalid = false;
        for (&base, &q) in sequence.iter().zip(quality) {
            let q = q.wrapping_sub(self.offset) as usize;
            invalid |= q > MAX_QSCORE;
            let q = q.min(MAX_QSCORE);
            counts[BASE_INDEX[base as usize] as usize] += 1;
            self.qualities[q] += 1;
            quality_sum += q as u64;
        }
        self.invalid_quality |= invalid;
        (counts, quality_sum)
    }

    /// Count bases and qualities at each read position.
    fn count_positions(&mut self, sequence: &[u8], quality: &[u8]) -> ([u32; BASE_CLASSES], u64) {
        let len = sequence.len();
        let bins = position_bin(len - 1) + 1;
        if self.positions.len() < bins {
            self.positions.resize(bins, PositionCounts::new());
            self.adapters.iter_mut().for_each(|a| a.resize(bins, 0));
        }
        let mut counts = [0u32; BASE_CLASSES];
        let mut quality_sum = 0u64;
        let mut invalid = false;
        let offset = self.offset;
        let mut count = |base: u8, q: u8, position: &mut PositionCounts| {
            let q = q.wrapping_sub(offset) as usize;
            invalid |= q > MAX_QSCORE;
            let q = q.min(MAX_QSCORE);
            let base = BASE_INDEX[base as usize] as usize;
            counts[base] += 1;
            quality_sum += q as u64;
            position.bases[base] += 1;
            position.qualities[q] += 1;
        };
        let exact = len.min(EXACT_POSITIONS);
        // Zipped slices avoid bounds checks in the common short-read case.
        sequence[..exact]
            .iter()
            .zip(&quality[..exact])
            .zip(&mut self.positions[..exact])
            .for_each(|((&base, &q), position)| count(base, q, position));
        for i in exact..len {
            count(
                sequence[i],
                quality[i],
                &mut self.positions[position_bin(i)],
            );
        }
        self.invalid_quality |= invalid;
        (counts, quality_sum)
    }

    fn find_adapters(&mut self, sequence: &[u8]) {
        let mut has_adapter = false;
        ADAPTER_FINDERS
            .iter()
            .zip(self.adapters.iter_mut())
            .enumerate()
            .for_each(|(i, (finder, starts))| {
                if let Some(pos) = finder.find(sequence) {
                    starts[position_bin(pos)] += 1;
                    has_adapter |= i < LIBRARY_ADAPTERS;
                }
            });
        self.adapter_reads += u64::from(has_adapter);
    }

    /// Combine counts from another block of the same file.
    pub(crate) fn merge(mut self, mut other: Self) -> Self {
        self.reads += other.reads;
        add_arrays(&mut self.bases, &other.bases);
        self.lengths.merge(other.lengths);
        add_arrays(&mut self.qualities, &other.qualities);
        add_arrays(&mut self.read_qualities, &other.read_qualities);
        add_arrays(&mut self.gc, &other.gc);
        if other.positions.len() > self.positions.len() {
            std::mem::swap(&mut self.positions, &mut other.positions);
            std::mem::swap(&mut self.adapters, &mut other.adapters);
        }
        self.positions
            .iter_mut()
            .zip(&other.positions)
            .for_each(|(a, b)| a.merge(b));
        self.adapters
            .iter_mut()
            .zip(&other.adapters)
            .for_each(|(a, b)| add_arrays(a, b));
        self.adapter_reads += other.adapter_reads;
        self.heavy_hitters.merge(other.heavy_hitters);
        self.invalid_quality |= other.invalid_quality;
        self
    }

    /// Fill counts derived from per-position counts in Complete mode.
    pub(crate) fn finish(&mut self) {
        if self.complete {
            let mut qualities = [0; QSCORE_BINS];
            self.positions
                .iter()
                .for_each(|p| add_arrays(&mut qualities, &p.qualities));
            self.qualities = qualities;
        }
    }
}

/// Add one read with `gc` G and C bases out of `acgt` to a GC histogram.
///
/// Like FastQC, the read covers the percentages from `gc - 0.5` to
/// `gc + 0.5` bases and is split between bins by overlap. Rounding to the
/// nearest percentage instead leaves some bins empty for 150 bp reads,
/// which draws a saw-tooth curve.
fn add_gc(hist: &mut [u64; 101], gc: u32, acgt: u32) {
    let to_percent = |bases: f64| (100.0 * bases / acgt as f64).clamp(-0.5, 100.5);
    let low = to_percent(gc as f64 - 0.5);
    let high = to_percent(gc as f64 + 0.5);
    let width = high - low;
    let first = (low + 0.5).floor().max(0.0) as usize;
    let last = ((high + 0.5).ceil() as usize).min(hist.len());
    // Cumulative rounding keeps each read's weights summing to GC_WEIGHT.
    let mut assigned = 0;
    for (bin, count) in hist.iter_mut().enumerate().take(last).skip(first) {
        let covered = (high.min(bin as f64 + 0.5) - low).max(0.0) / width;
        let cumulative = (covered.min(1.0) * GC_WEIGHT as f64).round() as u64;
        *count += cumulative - assigned;
        assigned = cumulative;
    }
}

fn add_arrays(target: &mut [u64], values: &[u64]) {
    target.iter_mut().zip(values).for_each(|(a, b)| *a += b);
}

#[cfg(test)]
mod test {
    use super::*;

    const RECORDS: &[u8] = b"@r1\nACGTN\n+\nII5#!\n@r2\nGGCCA\n+\nIIIII\n@r3\nGGCCA\n+\nIIIII\n";

    fn accumulate(data: &[u8], complete: bool) -> ReadAccumulator {
        let mut acc = ReadAccumulator::new(complete, 33);
        acc.add_chunk(data, 0, None).unwrap();
        acc.finish();
        acc
    }

    #[test]
    fn test_counts() {
        let acc = accumulate(RECORDS, false);
        assert_eq!(acc.reads, 3);
        assert_eq!(acc.bases, [3, 5, 5, 1, 1, 0]);
        assert_eq!(acc.lengths.counts(), vec![(5, 3)]);
        // I = Q40, 5 = Q20, # = Q2, ! = Q0
        assert_eq!(acc.qualities[40], 12);
        assert_eq!(acc.qualities[20], 1);
        assert_eq!(acc.qualities[0], 1);
        // Mean of 40, 40, 20, 2, 0 is 20.4.
        assert_eq!(acc.read_qualities[20], 1);
        assert_eq!(acc.read_qualities[40], 2);
        // Three reads, centered on 50% GC (r1) and 80% GC (r2 and r3).
        assert_eq!(acc.gc.iter().sum::<u64>(), 3 * GC_WEIGHT);
        // r1 spreads over 37.5-62.5%, r2 and r3 over 70-90%.
        assert_eq!(acc.gc[50], 40);
        assert_eq!(acc.gc[80], 100);
        assert_eq!(acc.gc[65], 0);
    }

    #[test]
    fn test_complete_matches_default_totals() {
        let default = accumulate(RECORDS, false);
        let complete = accumulate(RECORDS, true);
        assert_eq!(default.qualities, complete.qualities);
        assert_eq!(default.bases, complete.bases);
        assert_eq!(complete.positions.len(), 5);
        assert_eq!(complete.positions[4].bases[BASE_N], 1);
        assert_eq!(complete.positions[4].qualities[0], 1);
        assert_eq!(complete.heavy_hitters.reads(), 3);
    }

    #[test]
    fn test_merge_matches_single_pass() {
        let whole = accumulate(RECORDS, true);
        // Start of the second record.
        let second = 1 + RECORDS[1..].iter().position(|b| *b == b'@').unwrap();
        let mut first = ReadAccumulator::new(true, 33);
        first.add_chunk(&RECORDS[..second], 0, None).unwrap();
        let mut rest = ReadAccumulator::new(true, 33);
        rest.add_chunk(&RECORDS[second..], 1, None).unwrap();
        let mut merged = first.merge(rest);
        merged.finish();
        assert_eq!(merged, whole);
    }

    #[test]
    fn test_gc_spread_evenly() {
        // Every GC count of 150 bp reads gives each bin the same weight.
        let mut hist = [0u64; 101];
        (0..=150).for_each(|gc| add_gc(&mut hist, gc, 150));
        assert_eq!(hist.iter().sum::<u64>(), 151 * GC_WEIGHT);
        let inner = &hist[1..100];
        let (min, max) = (inner.iter().min().unwrap(), inner.iter().max().unwrap());
        assert!(max - min <= 2, "uneven bins: {} to {}", min, max);
        // A 4 bp read with 2 GC covers 37.5% to 62.5%.
        let mut hist = [0u64; 101];
        add_gc(&mut hist, 2, 4);
        assert_eq!(hist[37], 0);
        assert_eq!(hist[38], 40);
        assert_eq!(hist[50], 40);
        assert_eq!(hist.iter().sum::<u64>(), GC_WEIGHT);
    }

    #[test]
    fn test_adapter_position() {
        let data = b"@r1\nTTTAGATCGGAAGAGCC\n+\nIIIIIIIIIIIIIIIII\n";
        let acc = accumulate(data, true);
        assert_eq!(acc.adapters[0][3], 1);
        assert!(acc.adapters[1].iter().all(|c| *c == 0));
        assert_eq!(acc.adapter_reads, 1);
    }

    #[test]
    fn test_long_read_bins() {
        assert_eq!(position_bin(999), 999);
        assert_eq!(position_bin(1_000), 1_000);
        assert_eq!(position_bin(2_999), 1_001);
        assert_eq!(bin_range(1_001), (2_001, 3_000));
        let seq = vec![b'A'; 2_500];
        let qual = vec![b'I'; 2_500];
        let mut acc = ReadAccumulator::new(true, 33);
        acc.add_record(&seq, &qual);
        assert_eq!(acc.positions.len(), 1_002);
        assert_eq!(acc.positions[1_001].bases[BASE_A], 500);
    }

    #[test]
    fn test_detect_offset() {
        assert_eq!(detect_offset(RECORDS), 33);
        assert_eq!(detect_offset(b"@r1\nACGT\n+\nhhB@\n"), 64);
        assert_eq!(detect_offset(b""), 33);
    }

    #[test]
    fn test_invalid_quality() {
        let mut acc = ReadAccumulator::new(false, 64);
        assert!(acc.add_chunk(b"@r1\nACGT\n+\nII!I\n", 0, None).is_err());
    }
}
