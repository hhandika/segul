//! Read statistics derived from accumulated counts.
//!
//! The report is written to `read-report.json` in Complete mode, so
//! applications can draw their own charts without reading the CSV files.
use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use super::accumulate::{
    ADAPTERS, BASE_A, BASE_C, BASE_G, BASE_N, BASE_T, GC_WEIGHT, LIBRARY_ADAPTERS, LOW_QSCORE,
    POLY_G, QSCORE_BINS, ReadAccumulator, bin_range,
};
use super::sketch::DuplicationSketch;

/// Read starts above this share of reads are overrepresented.
const OVERREPRESENTED_FRACTION: f64 = 0.001;
const MAX_OVERREPRESENTED: usize = 100;
const Q20: usize = 20;
const Q30: usize = 30;
/// Lower bound of each duplication level, with its label.
const DUPLICATION_LEVELS: [(u64, &str); 16] = [
    (1, "1"),
    (2, "2"),
    (3, "3"),
    (4, "4"),
    (5, "5"),
    (6, "6"),
    (7, "7"),
    (8, "8"),
    (9, "9"),
    (10, ">10"),
    (50, ">50"),
    (100, ">100"),
    (500, ">500"),
    (1_000, ">1k"),
    (5_000, ">5k"),
    (10_000, ">10k"),
];

/// Read statistics for every sample of a run.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ReadReport {
    pub samples: Vec<SampleReadReport>,
}

impl ReadReport {
    /// Read a report written by `segul read summary --mode complete`.
    pub fn from_path(path: &Path) -> Result<Self> {
        let text = fs::read_to_string(path)
            .with_context(|| format!("Failed reading {}", path.display()))?;
        Self::from_json(&text).with_context(|| format!("Failed parsing {}", path.display()))
    }

    pub fn from_json(json: &str) -> Result<Self> {
        Ok(serde_json::from_str(json)?)
    }

    pub fn to_json(&self) -> Result<String> {
        Ok(serde_json::to_string(self)?)
    }
}

/// Read statistics for one FASTQ file.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct SampleReadReport {
    pub file_path: String,
    pub file_name: String,
    /// Sample name used in per-sample output files (file name without its last extension).
    pub sample: String,
    pub metrics: ReadMetrics,
    /// Reads per length, shortest first.
    pub length_distribution: Vec<LengthCount>,
    /// Reads by mean quality score, rounded down. Index is the quality score.
    pub read_quality_distribution: Vec<u64>,
    /// Reads by GC percentage, 0 to 100. A read between two whole
    /// percentages counts partly toward each, so values can be fractional.
    pub gc_distribution: Vec<f64>,
    /// Normal distribution with the same mean and deviation as `gc_distribution`.
    pub gc_theoretical: Vec<f64>,
    /// Per-position statistics. Empty unless computed in Complete mode.
    pub positions: Vec<PositionSummary>,
    pub duplication: Option<DuplicationSummary>,
    pub overrepresented: Vec<OverrepresentedSequence>,
    /// Cumulative share of reads with each adapter, aligned with `positions`.
    pub adapters: Vec<AdapterContent>,
}

/// Summary statistics for one FASTQ file.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ReadMetrics {
    pub read_count: u64,
    pub base_count: u64,
    pub min_read_length: u64,
    pub max_read_length: u64,
    pub mean_read_length: f64,
    pub median_read_length: f64,
    pub stdev_read_length: f64,
    pub n50_read_length: u64,
    pub g_count: u64,
    pub c_count: u64,
    pub a_count: u64,
    pub t_count: u64,
    pub n_count: u64,
    /// G and C over all bases.
    pub gc_content: f64,
    pub at_content: f64,
    pub n_content: f64,
    /// Mean base quality score.
    pub mean_qscore: f64,
    pub min_qscore: u64,
    pub max_qscore: u64,
    /// Bases with quality below 20.
    pub low_qscore_bases: u64,
    pub q20_bases: u64,
    pub q30_bases: u64,
    pub q20_pct: f64,
    pub q30_pct: f64,
    /// Mean of per-read mean quality scores, each rounded down.
    pub mean_read_qscore: f64,
    /// Reads with mean quality below 20.
    pub low_quality_reads: u64,
    pub phred_offset: u8,
    /// Share of sampled reads that repeat an earlier read, in percent. Complete mode only.
    pub duplicate_pct_est: Option<f64>,
    /// Reads with any library adapter, in percent. Complete mode only.
    pub adapter_reads_pct: Option<f64>,
    pub poly_g_reads_pct: Option<f64>,
    pub overrepresented_count: Option<u64>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LengthCount {
    pub length: u64,
    pub reads: u64,
}

/// Base and quality statistics for one read position or bin of positions.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct PositionSummary {
    /// First and last 1-based position. Equal except for long-read bins.
    pub start: u64,
    pub end: u64,
    pub g_count: u64,
    pub c_count: u64,
    pub a_count: u64,
    pub t_count: u64,
    pub n_count: u64,
    pub mean_qscore: f64,
    pub min_qscore: u64,
    pub max_qscore: u64,
    pub p10_qscore: u64,
    pub q1_qscore: u64,
    pub median_qscore: u64,
    pub q3_qscore: u64,
    pub p90_qscore: u64,
}

impl PositionSummary {
    /// All bases at this position, including ambiguous ones.
    pub fn base_count(&self) -> u64 {
        self.g_count + self.c_count + self.a_count + self.t_count + self.n_count
    }
}

/// Duplication across every read of a file.
///
/// Counted exactly when the file has at most `sketch_size` distinct
/// sequences, or when `sketch_size` is zero. Otherwise it is estimated from
/// a uniform sample of `sketch_size` distinct sequences, with their copies
/// counted over the whole file.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct DuplicationSummary {
    /// Reads counted, which is every read with bases.
    pub reads: u64,
    /// Distinct sequences, estimated unless `exact`.
    pub distinct_sequences: f64,
    /// Reads that repeat another read, in percent.
    pub duplicate_pct: f64,
    /// Distinct sequences kept, or zero when all were kept.
    pub sketch_size: u64,
    pub exact: bool,
    pub levels: Vec<DuplicationLevel>,
}

/// Sequences seen a given number of times.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DuplicationLevel {
    pub label: String,
    /// Share of distinct sequences at this level, in percent.
    pub distinct_pct: f64,
    /// Share of reads at this level, in percent.
    pub reads_pct: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OverrepresentedSequence {
    /// Read start: the whole read up to 75 bp, else its first 50 bp.
    pub sequence: String,
    /// Reads starting with the sequence. A lower bound, short by at most
    /// `max_error_pct` of reads in files with many distinct starts.
    pub count: u64,
    /// Share of all reads, in percent.
    pub percent: f64,
    pub max_error_pct: f64,
    /// Adapter found in the sequence, or "No hit".
    pub source: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AdapterContent {
    pub name: String,
    pub sequence: String,
    /// Share of reads with the adapter starting at or before each position, in percent.
    pub cumulative_pct: Vec<f64>,
}

impl SampleReadReport {
    /// `duplicates` is given in Complete mode, which adds per-position,
    /// adapter, duplication, and overrepresented statistics.
    pub(crate) fn from_accumulator(
        path: &Path,
        acc: &ReadAccumulator,
        duplicates: Option<&DuplicationSketch>,
    ) -> Self {
        let complete = duplicates.is_some();
        let file_name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or_default()
            .to_string();
        let sample = path
            .file_stem()
            .and_then(|n| n.to_str())
            .unwrap_or_default()
            .to_string();
        let lengths = acc.lengths.counts();
        let positions = position_summaries(acc);
        let adapters = if complete {
            adapter_content(acc)
        } else {
            Vec::new()
        };
        let duplication = duplicates.map(|sketch| duplication(acc, sketch));
        let overrepresented = if complete {
            overrepresented(acc)
        } else {
            Vec::new()
        };
        let mut metrics = metrics(acc, &lengths);
        if complete {
            let reads = acc.reads.max(1) as f64;
            metrics.duplicate_pct_est = duplication.as_ref().map(|d| d.duplicate_pct);
            metrics.adapter_reads_pct = Some(100.0 * acc.adapter_reads as f64 / reads);
            let poly_g: u64 = acc.adapters[POLY_G].iter().sum();
            metrics.poly_g_reads_pct = Some(100.0 * poly_g as f64 / reads);
            metrics.overrepresented_count = Some(overrepresented.len() as u64);
        }
        let read_quality_distribution = trim_zeros(&acc.read_qualities);
        Self {
            file_path: path.display().to_string(),
            file_name,
            sample,
            metrics,
            length_distribution: lengths
                .iter()
                .map(|(length, reads)| LengthCount {
                    length: *length as u64,
                    reads: *reads,
                })
                .collect(),
            read_quality_distribution,
            gc_distribution: acc
                .gc
                .iter()
                .map(|w| *w as f64 / GC_WEIGHT as f64)
                .collect(),
            gc_theoretical: normal_curve(&acc.gc, GC_WEIGHT),
            positions,
            duplication,
            overrepresented,
            adapters,
        }
    }
}

fn metrics(acc: &ReadAccumulator, lengths: &[(usize, u64)]) -> ReadMetrics {
    let base_count: u64 = acc.bases.iter().sum();
    let total = base_count as f64;
    let bases = &acc.bases;
    let gc = bases[BASE_G] + bases[BASE_C];
    let at = bases[BASE_A] + bases[BASE_T];
    let quality_count: u64 = acc.qualities.iter().sum();
    let quality_sum: u64 = weighted_sum(&acc.qualities);
    let q20_bases: u64 = acc.qualities[Q20..].iter().sum();
    let q30_bases: u64 = acc.qualities[Q30..].iter().sum();
    let scored_reads: u64 = acc.read_qualities.iter().sum();
    let mean_length = ratio(base_count as f64, acc.reads as f64);
    ReadMetrics {
        read_count: acc.reads,
        base_count,
        min_read_length: lengths.first().map(|l| l.0 as u64).unwrap_or(0),
        max_read_length: lengths.last().map(|l| l.0 as u64).unwrap_or(0),
        mean_read_length: mean_length,
        median_read_length: median_length(lengths, acc.reads),
        stdev_read_length: stdev_length(lengths, acc.reads, mean_length),
        n50_read_length: n50(lengths, base_count),
        g_count: bases[BASE_G],
        c_count: bases[BASE_C],
        a_count: bases[BASE_A],
        t_count: bases[BASE_T],
        n_count: bases[BASE_N],
        gc_content: ratio(gc as f64, total),
        at_content: ratio(at as f64, total),
        n_content: ratio(bases[BASE_N] as f64, total),
        mean_qscore: ratio(quality_sum as f64, quality_count as f64),
        min_qscore: first_nonzero(&acc.qualities).unwrap_or(0) as u64,
        max_qscore: last_nonzero(&acc.qualities).unwrap_or(0) as u64,
        low_qscore_bases: acc.qualities[..LOW_QSCORE].iter().sum(),
        q20_bases,
        q30_bases,
        q20_pct: 100.0 * ratio(q20_bases as f64, quality_count as f64),
        q30_pct: 100.0 * ratio(q30_bases as f64, quality_count as f64),
        mean_read_qscore: ratio(
            weighted_sum(&acc.read_qualities) as f64,
            scored_reads as f64,
        ),
        low_quality_reads: acc.read_qualities[..LOW_QSCORE].iter().sum(),
        phred_offset: acc.offset(),
        ..Default::default()
    }
}

fn position_summaries(acc: &ReadAccumulator) -> Vec<PositionSummary> {
    acc.positions
        .iter()
        .enumerate()
        .map(|(bin, counts)| {
            let (start, end) = bin_range(bin);
            let q = &counts.qualities;
            let n: u64 = q.iter().sum();
            PositionSummary {
                start: start as u64,
                end: end as u64,
                g_count: counts.bases[BASE_G],
                c_count: counts.bases[BASE_C],
                a_count: counts.bases[BASE_A],
                t_count: counts.bases[BASE_T],
                n_count: counts.bases[BASE_N],
                mean_qscore: ratio(weighted_sum(q) as f64, n as f64),
                min_qscore: first_nonzero(q).unwrap_or(0) as u64,
                max_qscore: last_nonzero(q).unwrap_or(0) as u64,
                p10_qscore: quantile(q, 0.10),
                q1_qscore: quantile(q, 0.25),
                median_qscore: quantile(q, 0.50),
                q3_qscore: quantile(q, 0.75),
                p90_qscore: quantile(q, 0.90),
            }
        })
        .collect()
}

fn adapter_content(acc: &ReadAccumulator) -> Vec<AdapterContent> {
    let reads = acc.reads.max(1) as f64;
    ADAPTERS
        .iter()
        .zip(&acc.adapters)
        .map(|((name, sequence), starts)| {
            let mut total = 0;
            let cumulative_pct = starts
                .iter()
                .map(|count| {
                    total += count;
                    100.0 * total as f64 / reads
                })
                .collect();
            AdapterContent {
                name: name.to_string(),
                sequence: sequence.to_string(),
                cumulative_pct,
            }
        })
        .collect()
}

fn duplication(acc: &ReadAccumulator, sketch: &DuplicationSketch) -> DuplicationSummary {
    let mut distinct = [0u64; DUPLICATION_LEVELS.len()];
    let mut copies = [0u64; DUPLICATION_LEVELS.len()];
    sketch.copies().into_iter().for_each(|count| {
        let level = DUPLICATION_LEVELS
            .iter()
            .rposition(|(lower, _)| count >= *lower)
            .unwrap_or(0);
        distinct[level] += 1;
        copies[level] += count;
    });
    let kept: u64 = distinct.iter().sum();
    let kept_copies: u64 = copies.iter().sum();
    let levels = DUPLICATION_LEVELS
        .iter()
        .enumerate()
        .map(|(i, (_, label))| DuplicationLevel {
            label: label.to_string(),
            distinct_pct: 100.0 * ratio(distinct[i] as f64, kept as f64),
            reads_pct: 100.0 * ratio(copies[i] as f64, kept_copies as f64),
        })
        .collect();
    // Reads without bases have no sequence to compare.
    let reads = acc.reads
        - acc
            .lengths
            .counts()
            .first()
            .map_or(0, |(len, count)| if *len == 0 { *count } else { 0 });
    let distinct_sequences = sketch.distinct().min(reads as f64);
    DuplicationSummary {
        reads,
        distinct_sequences,
        duplicate_pct: 100.0 * ratio(reads as f64 - distinct_sequences, reads as f64),
        sketch_size: sketch.capacity() as u64,
        exact: sketch.is_exact(),
        levels,
    }
}

fn overrepresented(acc: &ReadAccumulator) -> Vec<OverrepresentedSequence> {
    let hitters = &acc.heavy_hitters;
    let reads = hitters.reads() as f64;
    let error = hitters.max_error();
    let threshold = (reads * OVERREPRESENTED_FRACTION).max(1.0);
    // Unique sequences are never overrepresented, even in tiny files.
    let mut found: Vec<(&[u8], u64)> = hitters
        .counters()
        .filter(|c| c.count > 1 && (c.count + error) as f64 > threshold)
        .map(|c| (c.sequence(), c.count))
        .collect();
    found.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(b.0)));
    found
        .into_iter()
        .take(MAX_OVERREPRESENTED)
        .map(|(seq, count)| OverrepresentedSequence {
            sequence: String::from_utf8_lossy(seq).to_string(),
            count,
            percent: 100.0 * count as f64 / reads,
            max_error_pct: 100.0 * error as f64 / reads,
            source: adapter_source(seq),
        })
        .collect()
}

fn adapter_source(sequence: &[u8]) -> String {
    ADAPTERS[..LIBRARY_ADAPTERS]
        .iter()
        .find(|(_, adapter)| memchr::memmem::find(sequence, adapter.as_bytes()).is_some())
        .map(|(name, _)| name.to_string())
        .unwrap_or_else(|| String::from("No hit"))
}

/// Normal curve fitted to the GC histogram, scaled to the same number of reads.
///
/// `weight` is the histogram count of one read.
fn normal_curve(gc: &[u64], weight: u64) -> Vec<f64> {
    let total: u64 = gc.iter().sum();
    if total == 0 {
        return vec![0.0; gc.len()];
    }
    let mean = weighted_sum(gc) as f64 / total as f64;
    let variance = gc
        .iter()
        .enumerate()
        .map(|(x, c)| *c as f64 * (x as f64 - mean).powi(2))
        .sum::<f64>()
        / total as f64;
    let sd = variance.sqrt();
    let reads = total as f64 / weight as f64;
    if sd == 0.0 {
        return gc.iter().map(|c| *c as f64 / weight as f64).collect();
    }
    let scale = reads / (sd * (2.0 * std::f64::consts::PI).sqrt());
    (0..gc.len())
        .map(|x| scale * (-(x as f64 - mean).powi(2) / (2.0 * variance)).exp())
        .collect()
}

fn median_length(lengths: &[(usize, u64)], reads: u64) -> f64 {
    if reads == 0 {
        return 0.0;
    }
    let lower = nth_length(lengths, (reads - 1) / 2);
    let upper = nth_length(lengths, reads / 2);
    (lower + upper) as f64 / 2.0
}

/// Length of the read at a 0-based rank, shortest first.
fn nth_length(lengths: &[(usize, u64)], rank: u64) -> usize {
    let mut seen = 0;
    for (len, count) in lengths {
        seen += count;
        if seen > rank {
            return *len;
        }
    }
    0
}

fn stdev_length(lengths: &[(usize, u64)], reads: u64, mean: f64) -> f64 {
    if reads < 2 {
        return 0.0;
    }
    let squares: f64 = lengths
        .iter()
        .map(|(len, count)| *count as f64 * (*len as f64 - mean).powi(2))
        .sum();
    (squares / (reads - 1) as f64).sqrt()
}

/// Shortest read length that covers half the bases with reads at least that long.
fn n50(lengths: &[(usize, u64)], bases: u64) -> u64 {
    let mut covered = 0;
    for (len, count) in lengths.iter().rev() {
        covered += *len as u64 * count;
        if covered * 2 >= bases {
            return *len as u64;
        }
    }
    0
}

/// Quality score at a quantile of a histogram, by nearest rank.
fn quantile(hist: &[u64; QSCORE_BINS], q: f64) -> u64 {
    let n: u64 = hist.iter().sum();
    if n == 0 {
        return 0;
    }
    let rank = ((q * n as f64).ceil() as u64).max(1);
    let mut seen = 0;
    for (score, count) in hist.iter().enumerate() {
        seen += count;
        if seen >= rank {
            return score as u64;
        }
    }
    0
}

fn weighted_sum(hist: &[u64]) -> u64 {
    hist.iter()
        .enumerate()
        .map(|(value, count)| value as u64 * count)
        .sum()
}

fn first_nonzero(hist: &[u64]) -> Option<usize> {
    hist.iter().position(|c| *c > 0)
}

fn last_nonzero(hist: &[u64]) -> Option<usize> {
    hist.iter().rposition(|c| *c > 0)
}

fn trim_zeros(hist: &[u64]) -> Vec<u64> {
    match last_nonzero(hist) {
        Some(last) => hist[..=last].to_vec(),
        None => Vec::new(),
    }
}

fn ratio(numerator: f64, denominator: f64) -> f64 {
    if denominator == 0.0 {
        0.0
    } else {
        numerator / denominator
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_length_statistics() {
        let lengths = vec![(2, 1), (4, 2), (10, 1)];
        assert_eq!(median_length(&lengths, 4), 4.0);
        assert_eq!(median_length(&[(3, 1), (5, 1)], 2), 4.0);
        // 20 bases in total; the 10 bp read covers half.
        assert_eq!(n50(&lengths, 20), 10);
        assert_eq!(n50(&[(5, 4)], 20), 5);
        let mean = 5.0;
        let expected = ((9.0 + 1.0 + 1.0 + 25.0) / 3.0f64).sqrt();
        assert!((stdev_length(&lengths, 4, mean) - expected).abs() < 1e-12);
    }

    #[test]
    fn test_quantile() {
        let mut hist = [0u64; QSCORE_BINS];
        hist[10] = 1;
        hist[20] = 2;
        hist[30] = 1;
        assert_eq!(quantile(&hist, 0.10), 10);
        assert_eq!(quantile(&hist, 0.25), 10);
        assert_eq!(quantile(&hist, 0.50), 20);
        assert_eq!(quantile(&hist, 0.75), 20);
        assert_eq!(quantile(&hist, 0.90), 30);
    }

    #[test]
    fn test_normal_curve_area() {
        let mut gc = vec![0u64; 101];
        gc[40] = 10;
        gc[50] = 20;
        gc[60] = 10;
        let curve = normal_curve(&gc, 1);
        let area: f64 = curve.iter().sum();
        assert!((area - 40.0).abs() < 0.5);
        assert!(curve[50] > curve[40]);
    }
}
