//! Single-pass accumulator for FASTQ statistics.
//!
//! Every field is an integer count, fixed-point sum, or histogram, so
//! accumulators built from different blocks of a file merge by addition. The
//! result does not depend on how the file was split between threads.
use std::collections::BTreeMap;
use std::sync::LazyLock;

use anyhow::{Result, bail};
use memchr::memmem::Finder;

use super::header::HeaderAccumulator;
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
/// Number of equal-width relative-position bins along each read (5% each).
pub(crate) const RELATIVE_POSITION_BINS: usize = 20;
/// Histogram bins for normalized dinucleotide sequence complexity (5% each).
pub(crate) const COMPLEXITY_BINS: usize = 20;
/// Reads with normalized dinucleotide entropy below this are low complexity.
pub(crate) const LOW_COMPLEXITY_THRESHOLD: f64 = 0.50;
/// Minimum run length counted as a homopolymer.
pub(crate) const HOMOPOLYMER_MIN_LEN: usize = 6;
/// Minimum 3'-terminal homopolymer run counted as a poly-base tail.
pub(crate) const TERMINAL_HOMOPOLYMER_MIN_LEN: usize = 10;
/// Fixed-point multiplier so per-read float sums stay exact integers on merge.
pub(crate) const FIXED_SCALE: u128 = 1_000_000;
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
pub(crate) const BASE_OTHER: usize = 5;

/// Upper bounds and labels of per-read expected-error buckets.
pub(crate) const EE_BOUNDS: [(f64, &str); 8] = [
    (0.1, "<0.1"),
    (0.25, "0.1-0.25"),
    (0.5, "0.25-0.5"),
    (1.0, "0.5-1"),
    (2.0, "1-2"),
    (5.0, "2-5"),
    (10.0, "5-10"),
    (f64::INFINITY, ">=10"),
];
pub(crate) const EE_BINS: usize = EE_BOUNDS.len();

/// Length ranges used for the joint length x quality summary.
pub(crate) const LENGTH_QUALITY_BOUNDS: [(&str, u64, Option<u64>); 14] = [
    ("1-49", 1, Some(49)),
    ("50-99", 50, Some(99)),
    ("100-149", 100, Some(149)),
    ("150-199", 150, Some(199)),
    ("200-299", 200, Some(299)),
    ("300-499", 300, Some(499)),
    ("500-999", 500, Some(999)),
    ("1k-2k", 1_000, Some(1_999)),
    ("2k-5k", 2_000, Some(4_999)),
    ("5k-10k", 5_000, Some(9_999)),
    ("10k-25k", 10_000, Some(24_999)),
    ("25k-50k", 25_000, Some(49_999)),
    ("50k-100k", 50_000, Some(99_999)),
    ("100k+", 100_000, None),
];
pub(crate) const LENGTH_QUALITY_BINS: usize = LENGTH_QUALITY_BOUNDS.len();

static BASE_INDEX: [u8; 256] = build_base_index();

/// Per-base error probability `10^(-Q/10)` indexed by Phred score.
pub(crate) static ERROR_PROB: LazyLock<[f64; QSCORE_BINS]> = LazyLock::new(|| {
    let mut table = [0.0; QSCORE_BINS];
    for (q, slot) in table.iter_mut().enumerate() {
        *slot = 10f64.powf(-(q as f64) / 10.0);
    }
    table
});

/// `-ln(1 - 10^(-Q/10))` for `Q >= 1` (`0.0` placeholder at `Q == 0`).
static NEG_LN_ONE_MINUS_P: LazyLock<[f64; QSCORE_BINS]> = LazyLock::new(|| {
    let mut table = [0.0; QSCORE_BINS];
    for (q, slot) in table.iter_mut().enumerate().skip(1) {
        let p = 10f64.powf(-(q as f64) / 10.0);
        *slot = -(1.0 - p).ln();
    }
    table
});

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

/// Convert a non-negative floating-point value to fixed-point units for exact merges.
fn to_fixed(value: f64) -> u128 {
    (value.max(0.0) * FIXED_SCALE as f64).round() as u128
}

/// Convert fixed-point units back to floating-point.
pub(crate) fn from_fixed(value: u128) -> f64 {
    value as f64 / FIXED_SCALE as f64
}

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

/// Bin of a read length for the joint length x quality table.
fn length_quality_bin(len: usize) -> usize {
    let len = len as u64;
    LENGTH_QUALITY_BOUNDS
        .iter()
        .position(|(_, _, max)| max.is_none_or(|hi| len <= hi))
        .unwrap_or(LENGTH_QUALITY_BINS - 1)
}

/// Bucket of a per-read expected-error count.
fn expected_error_bin(ee: f64) -> usize {
    EE_BOUNDS
        .iter()
        .position(|(upper, _)| ee < *upper)
        .unwrap_or(EE_BINS - 1)
}

/// Detect the quality offset from the first reads of a block.
///
/// Any score below '@' means Phred+33, because Phred+64 data cannot contain it.
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

/// Read quality counts accumulated within one length range.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct LengthQualityCounts {
    pub(crate) reads: u64,
    pub(crate) bases: u64,
    pub(crate) base_quality_sum: u64,
    pub(crate) expected_errors_sum: u128,
    pub(crate) read_qualities: [u64; QSCORE_BINS],
}

impl LengthQualityCounts {
    fn new() -> Self {
        Self {
            reads: 0,
            bases: 0,
            base_quality_sum: 0,
            expected_errors_sum: 0,
            read_qualities: [0; QSCORE_BINS],
        }
    }

    fn merge(&mut self, other: &Self) {
        self.reads += other.reads;
        self.bases += other.bases;
        self.base_quality_sum += other.base_quality_sum;
        self.expected_errors_sum += other.expected_errors_sum;
        add_arrays(&mut self.read_qualities, &other.read_qualities);
    }
}

/// Single-pass per-read scan outputs returned from base counting.
struct PerReadScan {
    bases: [u32; BASE_CLASSES],
    quality_sum: u64,
    expected_errors: f64,
    error_free_prob: f64,
    complexity: f64,
    has_homopolymer: bool,
    homopolymer_bases: [u32; 4],
    max_homopolymer_len: u32,
    homopolymer_q_sum: u64,
    homopolymer_q_count: u64,
    non_homopolymer_q_sum: u64,
    non_homopolymer_q_count: u64,
    terminal_homopolymer: Option<usize>,
}

/// Running state for per-read base, error, complexity, and homopolymer stats.
struct ReadScanState {
    bases: [u32; BASE_CLASSES],
    quality_sum: u64,
    expected_errors: f64,
    neg_ln_one_minus_p: f64,
    has_q0: bool,
    dinuc: [u32; 16],
    valid_dinucs: u32,
    prev_base: usize,
    run_base: usize,
    run_len: u32,
    run_q_sum: u64,
    has_homopolymer: bool,
    homopolymer_bases: [u32; 4],
    max_homopolymer_len: u32,
    homopolymer_q_sum: u64,
    homopolymer_q_count: u64,
    non_homopolymer_q_sum: u64,
    non_homopolymer_q_count: u64,
    invalid: bool,
}

impl ReadScanState {
    fn new() -> Self {
        Self {
            bases: [0; BASE_CLASSES],
            quality_sum: 0,
            expected_errors: 0.0,
            neg_ln_one_minus_p: 0.0,
            has_q0: false,
            dinuc: [0; 16],
            valid_dinucs: 0,
            prev_base: BASE_OTHER,
            run_base: BASE_OTHER,
            run_len: 0,
            run_q_sum: 0,
            has_homopolymer: false,
            homopolymer_bases: [0; 4],
            max_homopolymer_len: 0,
            homopolymer_q_sum: 0,
            homopolymer_q_count: 0,
            non_homopolymer_q_sum: 0,
            non_homopolymer_q_count: 0,
            invalid: false,
        }
    }

    #[inline(always)]
    fn push(
        &mut self,
        base: u8,
        raw_q: u8,
        offset: u8,
        err_prob: &[f64; QSCORE_BINS],
        neg_ln_lut: &[f64; QSCORE_BINS],
    ) -> (usize, usize) {
        let q = raw_q.wrapping_sub(offset) as usize;
        self.invalid |= q > MAX_QSCORE;
        let q = q.min(MAX_QSCORE);
        let base = BASE_INDEX[base as usize] as usize;
        self.bases[base] += 1;
        self.quality_sum += q as u64;
        self.expected_errors += err_prob[q];
        if q == 0 {
            self.has_q0 = true;
        } else {
            self.neg_ln_one_minus_p += neg_ln_lut[q];
        }
        if self.prev_base < 4 && base < 4 {
            self.dinuc[self.prev_base * 4 + base] += 1;
            self.valid_dinucs += 1;
        }
        self.prev_base = base;
        if base == self.run_base {
            self.run_len += 1;
            self.run_q_sum += q as u64;
        } else {
            self.flush_run();
            self.run_base = base;
            self.run_len = 1;
            self.run_q_sum = q as u64;
        }
        (base, q)
    }

    #[inline(always)]
    fn flush_run(&mut self) {
        if self.run_base < 4 && self.run_len > 0 {
            if self.run_len > self.max_homopolymer_len {
                self.max_homopolymer_len = self.run_len;
            }
            if self.run_len as usize >= HOMOPOLYMER_MIN_LEN {
                self.has_homopolymer = true;
                self.homopolymer_bases[self.run_base] += self.run_len;
                self.homopolymer_q_sum += self.run_q_sum;
                self.homopolymer_q_count += u64::from(self.run_len);
            } else {
                self.non_homopolymer_q_sum += self.run_q_sum;
                self.non_homopolymer_q_count += u64::from(self.run_len);
            }
        }
    }

    fn finish(mut self) -> (PerReadScan, bool) {
        let terminal_homopolymer = (self.run_base < 4
            && self.run_len as usize >= TERMINAL_HOMOPOLYMER_MIN_LEN)
            .then_some(self.run_base);
        self.flush_run();
        let error_free_prob = if self.has_q0 {
            0.0
        } else {
            (-self.neg_ln_one_minus_p).exp()
        };
        let complexity = dinucleotide_complexity(&self.dinuc, self.valid_dinucs);
        let invalid = self.invalid;
        let scan = PerReadScan {
            bases: self.bases,
            quality_sum: self.quality_sum,
            expected_errors: self.expected_errors,
            error_free_prob,
            complexity,
            has_homopolymer: self.has_homopolymer,
            homopolymer_bases: self.homopolymer_bases,
            max_homopolymer_len: self.max_homopolymer_len,
            homopolymer_q_sum: self.homopolymer_q_sum,
            homopolymer_q_count: self.homopolymer_q_count,
            non_homopolymer_q_sum: self.non_homopolymer_q_sum,
            non_homopolymer_q_count: self.non_homopolymer_q_count,
            terminal_homopolymer,
        };
        (scan, invalid)
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
    /// Reads by their mean quality, rounded to the nearest integer.
    pub(crate) read_qualities: [u64; QSCORE_BINS],
    /// Sum of exact per-read arithmetic mean quality scores (scaled by `FIXED_SCALE`).
    pub(crate) read_qscore_sum: u128,
    /// Sum of per-read expected error counts (scaled by `FIXED_SCALE`).
    pub(crate) expected_errors_sum: u128,
    /// Sum of per-read error-free probabilities (scaled by `FIXED_SCALE`).
    pub(crate) error_free_prob_sum: u128,
    /// Reads by expected-error bucket.
    pub(crate) ee_hist: [u64; EE_BINS],
    /// Bases belonging to reads with mean quality at or above 10, 20, and 30.
    pub(crate) q10_read_bases: u64,
    pub(crate) q20_read_bases: u64,
    pub(crate) q30_read_bases: u64,
    /// Reads with exact mean quality below 20.
    pub(crate) low_quality_reads: u64,
    /// Sum of normalized dinucleotide complexity per read (scaled by `FIXED_SCALE`).
    pub(crate) complexity_sum: u128,
    pub(crate) low_complexity_reads: u64,
    pub(crate) complexity_hist: [u64; COMPLEXITY_BINS],
    /// Bases in homopolymer runs of length >= 6, indexed by A, C, G, T.
    pub(crate) homopolymer_bases: [u64; 4],
    pub(crate) homopolymer_reads: u64,
    /// Reads ending in a 3'-terminal homopolymer of length >= 10, indexed by A, C, G, T.
    pub(crate) terminal_homopolymer_reads: [u64; 4],
    pub(crate) homopolymer_q_sum: u64,
    pub(crate) homopolymer_q_count: u64,
    pub(crate) non_homopolymer_q_sum: u64,
    pub(crate) non_homopolymer_q_count: u64,
    pub(crate) max_homopolymer_length: u64,
    /// Metadata parsed from FASTQ header lines.
    pub(crate) header: HeaderAccumulator,
    /// Reads by GC percentage of their A, C, G, and T bases, in units of
    /// [`GC_WEIGHT`]. See [`add_gc`] for how a read is split between bins.
    pub(crate) gc: [u64; 101],
    pub(crate) positions: Vec<PositionCounts>,
    /// Base and quality counts across 20 equal relative-position bins (Complete mode).
    pub(crate) relative_positions: Vec<PositionCounts>,
    /// Quality counts stratified by read-length range (Complete mode).
    pub(crate) length_quality: Vec<LengthQualityCounts>,
    /// Reads by the bin where each adapter first starts.
    pub(crate) adapters: Vec<Vec<u64>>,
    /// Reads with any library adapter.
    pub(crate) adapter_reads: u64,
    /// Most common read starts, for overrepresented sequences.
    pub(crate) heavy_hitters: HeavyHitters,
    invalid_quality: bool,
}

impl ReadAccumulator {
    /// `complete` adds per-position, relative-position, length-quality,
    /// adapter, and overrepresented counts.
    pub(crate) fn new(complete: bool, offset: u8) -> Self {
        let relative_positions = if complete {
            vec![PositionCounts::new(); RELATIVE_POSITION_BINS]
        } else {
            Vec::new()
        };
        let length_quality = if complete {
            vec![LengthQualityCounts::new(); LENGTH_QUALITY_BINS]
        } else {
            Vec::new()
        };
        Self {
            complete,
            offset,
            reads: 0,
            bases: [0; BASE_CLASSES],
            lengths: LengthHistogram::default(),
            qualities: [0; QSCORE_BINS],
            read_qualities: [0; QSCORE_BINS],
            read_qscore_sum: 0,
            expected_errors_sum: 0,
            error_free_prob_sum: 0,
            ee_hist: [0; EE_BINS],
            q10_read_bases: 0,
            q20_read_bases: 0,
            q30_read_bases: 0,
            low_quality_reads: 0,
            complexity_sum: 0,
            low_complexity_reads: 0,
            complexity_hist: [0; COMPLEXITY_BINS],
            homopolymer_bases: [0; 4],
            homopolymer_reads: 0,
            terminal_homopolymer_reads: [0; 4],
            homopolymer_q_sum: 0,
            homopolymer_q_count: 0,
            non_homopolymer_q_sum: 0,
            non_homopolymer_q_count: 0,
            max_homopolymer_length: 0,
            header: HeaderAccumulator::new(),
            gc: [0; 101],
            positions: Vec::new(),
            relative_positions,
            length_quality,
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
            self.header.add(record.header);
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
        let scan = if self.complete {
            self.count_positions(sequence, quality)
        } else {
            self.count_bases(sequence, quality)
        };
        let len_u64 = len as u64;
        add_arrays(&mut self.bases, &scan.bases.map(u64::from));
        let exact_mean_q = scan.quality_sum as f64 / len_u64 as f64;
        let rounded_q = ((scan.quality_sum + len_u64 / 2) / len_u64) as usize;
        let rounded_q = rounded_q.min(MAX_QSCORE);
        self.read_qualities[rounded_q] += 1;
        self.read_qscore_sum += to_fixed(exact_mean_q);
        if scan.quality_sum >= 10 * len_u64 {
            self.q10_read_bases += len_u64;
        }
        if scan.quality_sum >= 20 * len_u64 {
            self.q20_read_bases += len_u64;
        } else {
            self.low_quality_reads += 1;
        }
        if scan.quality_sum >= 30 * len_u64 {
            self.q30_read_bases += len_u64;
        }
        let ee_fixed = to_fixed(scan.expected_errors);
        self.expected_errors_sum += ee_fixed;
        self.error_free_prob_sum += to_fixed(scan.error_free_prob);
        self.ee_hist[expected_error_bin(scan.expected_errors)] += 1;

        self.complexity_sum += to_fixed(scan.complexity);
        if scan.complexity < LOW_COMPLEXITY_THRESHOLD {
            self.low_complexity_reads += 1;
        }
        let cbin =
            ((scan.complexity * COMPLEXITY_BINS as f64).floor() as usize).min(COMPLEXITY_BINS - 1);
        self.complexity_hist[cbin] += 1;

        for i in 0..4 {
            self.homopolymer_bases[i] += u64::from(scan.homopolymer_bases[i]);
        }
        if scan.has_homopolymer {
            self.homopolymer_reads += 1;
        }
        if let Some(base) = scan.terminal_homopolymer {
            self.terminal_homopolymer_reads[base] += 1;
        }
        self.homopolymer_q_sum += scan.homopolymer_q_sum;
        self.homopolymer_q_count += scan.homopolymer_q_count;
        self.non_homopolymer_q_sum += scan.non_homopolymer_q_sum;
        self.non_homopolymer_q_count += scan.non_homopolymer_q_count;
        self.max_homopolymer_length = self
            .max_homopolymer_length
            .max(u64::from(scan.max_homopolymer_len));

        let acgt =
            scan.bases[BASE_A] + scan.bases[BASE_C] + scan.bases[BASE_G] + scan.bases[BASE_T];
        if acgt > 0 {
            add_gc(&mut self.gc, scan.bases[BASE_G] + scan.bases[BASE_C], acgt);
        }
        if self.complete {
            let lq = &mut self.length_quality[length_quality_bin(len)];
            lq.reads += 1;
            lq.bases += len_u64;
            lq.base_quality_sum += scan.quality_sum;
            lq.expected_errors_sum += ee_fixed;
            lq.read_qualities[rounded_q] += 1;
            self.find_adapters(sequence);
            self.heavy_hitters.add(sequence);
        }
    }

    /// Count bases, qualities, errors, complexity, and homopolymers for the whole file.
    fn count_bases(&mut self, sequence: &[u8], quality: &[u8]) -> PerReadScan {
        let err_prob = &*ERROR_PROB;
        let neg_ln_lut = &*NEG_LN_ONE_MINUS_P;
        let offset = self.offset;
        let mut state = ReadScanState::new();
        for (&base, &q) in sequence.iter().zip(quality) {
            let (_, q) = state.push(base, q, offset, err_prob, neg_ln_lut);
            self.qualities[q] += 1;
        }
        let (scan, invalid) = state.finish();
        self.invalid_quality |= invalid;
        scan
    }

    /// Count bases and qualities at each absolute and relative read position.
    fn count_positions(&mut self, sequence: &[u8], quality: &[u8]) -> PerReadScan {
        let len = sequence.len();
        let bins = position_bin(len - 1) + 1;
        if self.positions.len() < bins {
            self.positions.resize(bins, PositionCounts::new());
            self.adapters.iter_mut().for_each(|a| a.resize(bins, 0));
        }
        let err_prob = &*ERROR_PROB;
        let neg_ln_lut = &*NEG_LN_ONE_MINUS_P;
        let offset = self.offset;
        let mut state = ReadScanState::new();
        let exact = len.min(EXACT_POSITIONS);
        let rel = &mut self.relative_positions;
        for (i, ((&base, &q), position)) in sequence[..exact]
            .iter()
            .zip(&quality[..exact])
            .zip(&mut self.positions[..exact])
            .enumerate()
        {
            let (base_idx, q_idx) = state.push(base, q, offset, err_prob, neg_ln_lut);
            position.bases[base_idx] += 1;
            position.qualities[q_idx] += 1;
            let rbin = (i * RELATIVE_POSITION_BINS) / len;
            rel[rbin].bases[base_idx] += 1;
            rel[rbin].qualities[q_idx] += 1;
        }
        for i in exact..len {
            let (base_idx, q_idx) =
                state.push(sequence[i], quality[i], offset, err_prob, neg_ln_lut);
            let position = &mut self.positions[position_bin(i)];
            position.bases[base_idx] += 1;
            position.qualities[q_idx] += 1;
            let rbin = (i * RELATIVE_POSITION_BINS) / len;
            rel[rbin].bases[base_idx] += 1;
            rel[rbin].qualities[q_idx] += 1;
        }
        let (scan, invalid) = state.finish();
        self.invalid_quality |= invalid;
        scan
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
        self.read_qscore_sum += other.read_qscore_sum;
        self.expected_errors_sum += other.expected_errors_sum;
        self.error_free_prob_sum += other.error_free_prob_sum;
        add_arrays(&mut self.ee_hist, &other.ee_hist);
        self.q10_read_bases += other.q10_read_bases;
        self.q20_read_bases += other.q20_read_bases;
        self.q30_read_bases += other.q30_read_bases;
        self.low_quality_reads += other.low_quality_reads;
        self.complexity_sum += other.complexity_sum;
        self.low_complexity_reads += other.low_complexity_reads;
        add_arrays(&mut self.complexity_hist, &other.complexity_hist);
        add_arrays(&mut self.homopolymer_bases, &other.homopolymer_bases);
        self.homopolymer_reads += other.homopolymer_reads;
        add_arrays(
            &mut self.terminal_homopolymer_reads,
            &other.terminal_homopolymer_reads,
        );
        self.homopolymer_q_sum += other.homopolymer_q_sum;
        self.homopolymer_q_count += other.homopolymer_q_count;
        self.non_homopolymer_q_sum += other.non_homopolymer_q_sum;
        self.non_homopolymer_q_count += other.non_homopolymer_q_count;
        self.max_homopolymer_length = self
            .max_homopolymer_length
            .max(other.max_homopolymer_length);
        self.header.merge(other.header);
        add_arrays(&mut self.gc, &other.gc);
        if other.positions.len() > self.positions.len() {
            std::mem::swap(&mut self.positions, &mut other.positions);
            std::mem::swap(&mut self.adapters, &mut other.adapters);
        }
        self.positions
            .iter_mut()
            .zip(&other.positions)
            .for_each(|(a, b)| a.merge(b));
        self.relative_positions
            .iter_mut()
            .zip(&other.relative_positions)
            .for_each(|(a, b)| a.merge(b));
        self.length_quality
            .iter_mut()
            .zip(&other.length_quality)
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

/// Normalized dinucleotide Shannon entropy in `[0.0, 1.0]`.
fn dinucleotide_complexity(dinuc: &[u32; 16], total: u32) -> f64 {
    if total <= 1 {
        return 0.0;
    }
    let total_f = f64::from(total);
    let entropy: f64 = dinuc
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| {
            let p = f64::from(c) / total_f;
            -p * p.log2()
        })
        .sum();
    let max_entropy = total_f.min(16.0).log2();
    if max_entropy <= 0.0 {
        0.0
    } else {
        (entropy / max_entropy).clamp(0.0, 1.0)
    }
}

/// Add one read with `gc` G and C bases out of `acgt` to a GC histogram.
///
/// The read covers the percentages from `gc - 0.5` to `gc + 0.5` bases and is
/// split between bins by overlap. Rounding to the nearest percentage instead
/// leaves some bins empty for 150 bp reads, which draws a saw-tooth curve.
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
        // Mean of 40, 40, 20, 2, 0 is 20.4, rounding to 20.
        assert_eq!(acc.read_qualities[20], 1);
        assert_eq!(acc.read_qualities[40], 2);
        assert_eq!(acc.low_quality_reads, 0);
        assert_eq!(acc.q20_read_bases, 15);
        assert_eq!(acc.q30_read_bases, 10);
        // Three reads, centered on 50% GC (r1) and 80% GC (r2 and r3).
        assert_eq!(acc.gc.iter().sum::<u64>(), 3 * GC_WEIGHT);
        // r1 spreads over 37.5-62.5%, r2 and r3 over 70-90%.
        assert_eq!(acc.gc[50], 40);
        assert_eq!(acc.gc[80], 100);
        assert_eq!(acc.gc[65], 0);
    }

    #[test]
    fn test_complete_matches_default_totals() {
        let basic = accumulate(RECORDS, false);
        let complete = accumulate(RECORDS, true);
        assert_eq!(basic.qualities, complete.qualities);
        assert_eq!(basic.bases, complete.bases);
        assert_eq!(basic.expected_errors_sum, complete.expected_errors_sum);
        assert_eq!(basic.error_free_prob_sum, complete.error_free_prob_sum);
        assert_eq!(basic.complexity_sum, complete.complexity_sum);
        assert_eq!(basic.homopolymer_bases, complete.homopolymer_bases);
        assert_eq!(complete.positions.len(), 5);
        assert_eq!(complete.relative_positions.len(), RELATIVE_POSITION_BINS);
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
    fn test_complexity_and_homopolymers() {
        // 12-base homopolymer tail (G) at Q20 vs 4 non-homopolymer bases at Q40.
        let data = b"@r1\nACGTGGGGGGGGGGGG\n+\nIIII555555555555\n@r2\nATATATATATATATAT\n+\nIIIIIIIIIIIIIIII\n";
        let acc = accumulate(data, true);
        assert_eq!(acc.homopolymer_reads, 1);
        assert_eq!(acc.homopolymer_bases[BASE_G], 12);
        assert_eq!(acc.terminal_homopolymer_reads[BASE_G], 1);
        assert_eq!(acc.max_homopolymer_length, 12);
        assert_eq!(acc.homopolymer_q_sum, 12 * 20);
        assert_eq!(acc.homopolymer_q_count, 12);
        // 4 non-homopolymer bases in r1 + 16 bases in r2, all Q40.
        assert_eq!(acc.non_homopolymer_q_count, 20);
        assert_eq!(acc.non_homopolymer_q_sum, 20 * 40);
        // Both reads are low complexity (homopolymer tail and dinucleotide repeat).
        assert_eq!(acc.low_complexity_reads, 2);
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
        assert_eq!(acc.length_quality[length_quality_bin(2_500)].reads, 1);
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
