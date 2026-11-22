//! Read statistics derived from accumulated counts.
//!
//! The report is written to `read-report.json` in Complete mode, so
//! applications can draw their own charts without reading the CSV files.
use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use super::accumulate::{
    ADAPTERS, BASE_A, BASE_C, BASE_G, BASE_N, BASE_OTHER, BASE_T, EE_BOUNDS, ERROR_PROB, GC_WEIGHT,
    LENGTH_QUALITY_BOUNDS, LIBRARY_ADAPTERS, LOW_QSCORE, MAX_QSCORE, POLY_G, PositionCounts,
    QSCORE_BINS, RELATIVE_POSITION_BINS, ReadAccumulator, bin_range, from_fixed,
};
use super::sketch::DuplicationSketch;

/// Read starts above this share of reads are overrepresented.
const OVERREPRESENTED_FRACTION: f64 = 0.001;
const MAX_OVERREPRESENTED: usize = 100;
const Q10: usize = 10;
const Q15: usize = 15;
const Q20: usize = 20;
const Q30: usize = 30;
const Q40: usize = 40;
/// Maximum 5' window inspected for leading base-composition bias (bp).
const HEAD_BIAS_WINDOW: usize = 30;
/// Percentage-point deviation from the read body that flags 5' composition bias.
const HEAD_BIAS_THRESHOLD_PCT: f64 = 5.0;
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
    /// Sample name used in per-sample output files (compound FASTQ extensions stripped).
    pub sample: String,
    pub metrics: ReadMetrics,
    /// Reads per length, shortest first.
    pub length_distribution: Vec<LengthCount>,
    /// Reads by mean quality score, rounded to the nearest integer. Index is the quality score.
    pub read_quality_distribution: Vec<u64>,
    /// Reads by GC percentage, 0 to 100. A read between two whole
    /// percentages counts partly toward each, so values can be fractional.
    pub gc_distribution: Vec<f64>,
    /// Normal distribution with the same mean and deviation as `gc_distribution`.
    pub gc_theoretical: Vec<f64>,
    /// Per-position statistics. Empty unless computed in Complete mode.
    pub positions: Vec<PositionSummary>,
    /// Quality and composition across 20 equal relative-position bins (0-100%).
    #[serde(default)]
    pub relative_positions: Vec<RelativePositionSummary>,
    /// Quality statistics stratified by read-length range.
    #[serde(default)]
    pub length_quality: Vec<LengthQualityBin>,
    /// Reads by expected-error bucket.
    #[serde(default)]
    pub expected_error_distribution: Vec<ExpectedErrorBucket>,
    /// Reads by normalized dinucleotide complexity in 20 equal bins (5% each).
    #[serde(default)]
    pub complexity_distribution: Vec<u64>,
    /// Homopolymer run and 3'-tail summary.
    #[serde(default)]
    pub homopolymers: Option<HomopolymerSummary>,
    /// Sequencing platform, filter, and run-group metadata parsed from headers.
    #[serde(default)]
    pub header: Option<HeaderSummary>,
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
    #[serde(default)]
    pub n90_read_length: u64,
    #[serde(default)]
    pub l50_read_length: u64,
    /// Area under the Nx curve (`sum(len^2) / sum(len)`).
    #[serde(default)]
    pub aun_read_length: f64,
    /// Reads with zero bases.
    #[serde(default)]
    pub empty_reads: u64,
    pub g_count: u64,
    pub c_count: u64,
    pub a_count: u64,
    pub t_count: u64,
    pub n_count: u64,
    /// Non-ACGTN bases (such as IUPAC ambiguity codes or RNA U).
    #[serde(default)]
    pub other_count: u64,
    /// G and C over all bases.
    pub gc_content: f64,
    pub at_content: f64,
    pub n_content: f64,
    /// Strand parity `(A - T) / (A + T)`.
    #[serde(default)]
    pub at_skew: f64,
    /// Strand parity `(G - C) / (G + C)`.
    #[serde(default)]
    pub gc_skew: f64,
    /// Standard deviation of per-read GC percentage.
    #[serde(default)]
    pub gc_stdev: f64,
    /// Total-variation divergence between observed and normal GC curves, in percent.
    #[serde(default)]
    pub gc_normal_divergence_pct: f64,
    /// Prominent peaks detected in the smoothed per-read GC distribution.
    #[serde(default)]
    pub gc_peaks: u64,
    /// Arithmetic mean base quality score.
    pub mean_qscore: f64,
    /// Error-rate mean quality score (`-10 * log10(mean_error_rate)`).
    #[serde(default)]
    pub effective_mean_qscore: f64,
    /// Mean per-base sequencing error probability (`mean(10^(-Q/10))`).
    #[serde(default)]
    pub mean_error_rate: f64,
    pub min_qscore: u64,
    pub max_qscore: u64,
    /// Bases with quality below 20.
    pub low_qscore_bases: u64,
    pub q20_bases: u64,
    pub q30_bases: u64,
    #[serde(default)]
    pub q10_pct: f64,
    #[serde(default)]
    pub q15_pct: f64,
    pub q20_pct: f64,
    pub q30_pct: f64,
    #[serde(default)]
    pub q40_pct: f64,
    /// Number of distinct Phred quality values observed.
    #[serde(default)]
    pub distinct_qscores: u64,
    /// Mean of exact per-read arithmetic mean quality scores.
    pub mean_read_qscore: f64,
    /// Reads with exact mean quality below 20.
    pub low_quality_reads: u64,
    /// Mean expected sequencing errors per read.
    #[serde(default)]
    pub mean_expected_errors: f64,
    /// Expected sequencing errors per 100 bases.
    #[serde(default)]
    pub expected_errors_per_100bp: f64,
    /// Estimated share of reads with zero sequencing errors, in percent.
    #[serde(default)]
    pub expected_error_free_reads_pct: f64,
    /// Share of reads with fewer than 1 expected error, in percent.
    #[serde(default)]
    pub reads_ee_lt_1_pct: f64,
    /// Share of bases in reads with mean quality >= 10, in percent.
    #[serde(default)]
    pub q10_read_bases_pct: f64,
    /// Share of bases in reads with mean quality >= 20, in percent.
    #[serde(default)]
    pub q20_read_bases_pct: f64,
    /// Share of bases in reads with mean quality >= 30, in percent.
    #[serde(default)]
    pub q30_read_bases_pct: f64,
    /// Mean normalized dinucleotide Shannon entropy across reads (`0.0` to `1.0`).
    #[serde(default)]
    pub mean_complexity: f64,
    /// Share of reads with dinucleotide complexity below 0.50, in percent.
    #[serde(default)]
    pub low_complexity_reads_pct: f64,
    /// Share of bases in homopolymer runs of 6 bp or longer, in percent.
    #[serde(default)]
    pub homopolymer_bases_pct: f64,
    /// Longest homopolymer run observed in the file.
    #[serde(default)]
    pub max_homopolymer_length: u64,
    /// Mean quality score inside homopolymer runs (>= 6 bp).
    #[serde(default)]
    pub homopolymer_mean_qscore: f64,
    /// Mean quality score of A/C/G/T bases outside homopolymer runs.
    #[serde(default)]
    pub non_homopolymer_mean_qscore: f64,
    /// `non_homopolymer_mean_qscore - homopolymer_mean_qscore`.
    #[serde(default)]
    pub homopolymer_qscore_gap: f64,
    /// Share of reads ending in a 3'-terminal poly-G run (>= 10 bp), in percent.
    #[serde(default)]
    pub poly_g_tail_reads_pct: f64,
    /// Share of reads ending in a 3'-terminal poly-A run (>= 10 bp), in percent.
    #[serde(default)]
    pub poly_a_tail_reads_pct: f64,
    /// Share of reads ending in a 3'-terminal poly-T run (>= 10 bp), in percent.
    #[serde(default)]
    pub poly_t_tail_reads_pct: f64,
    /// Share of reads ending in a 3'-terminal poly-C run (>= 10 bp), in percent.
    #[serde(default)]
    pub poly_c_tail_reads_pct: f64,
    /// Dominant sequencing platform inferred from FASTQ headers, if recognized.
    #[serde(default)]
    pub platform: Option<String>,
    /// Reads flagged as failing the Illumina chastity filter (`:Y:`).
    #[serde(default)]
    pub chastity_failed_reads: u64,
    /// Share of reads failing the Illumina chastity filter, in percent.
    #[serde(default)]
    pub chastity_failed_pct: f64,
    /// Number of distinct flowcell/lane or run groups seen in headers.
    #[serde(default)]
    pub run_group_count: u64,
    pub phred_offset: u8,
    /// Share of sampled reads that repeat an earlier read, in percent. Complete mode only.
    pub duplicate_pct_est: Option<f64>,
    /// Reads with any library adapter, in percent. Complete mode only.
    pub adapter_reads_pct: Option<f64>,
    pub poly_g_reads_pct: Option<f64>,
    pub overrepresented_count: Option<u64>,
    /// First-quarter minus last-quarter effective quality across relative read length. Complete mode only.
    #[serde(default)]
    pub qscore_tail_drop: Option<f64>,
    /// Longest leading span (bp) with per-position median quality >= 20. Complete mode only.
    #[serde(default)]
    pub effective_q20_span: Option<u64>,
    /// Longest leading span (bp) with per-position median quality >= 30. Complete mode only.
    #[serde(default)]
    pub effective_q30_span: Option<u64>,
    /// Leading positions (bp) with 5' base-composition bias > 5% vs read body. Complete mode only.
    #[serde(default)]
    pub head_bias_bp: Option<u64>,
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
    /// Error-rate mean quality score at this position.
    #[serde(default)]
    pub effective_mean_qscore: f64,
    /// Share of reads reaching this position (or average across bin), in percent.
    #[serde(default)]
    pub read_survival_pct: f64,
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

/// Base and quality statistics for one fractional read-length bin (e.g., 0-5%).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct RelativePositionSummary {
    pub start_pct: u64,
    pub end_pct: u64,
    pub g_count: u64,
    pub c_count: u64,
    pub a_count: u64,
    pub t_count: u64,
    pub n_count: u64,
    pub mean_qscore: f64,
    pub effective_mean_qscore: f64,
    pub p10_qscore: u64,
    pub q1_qscore: u64,
    pub median_qscore: u64,
    pub q3_qscore: u64,
    pub p90_qscore: u64,
}

/// Quality summary within one read-length range.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct LengthQualityBin {
    pub label: String,
    pub min_length: u64,
    pub max_length: Option<u64>,
    pub reads: u64,
    pub bases: u64,
    pub mean_qscore: f64,
    pub effective_mean_qscore: f64,
    pub q1_qscore: u64,
    pub median_qscore: u64,
    pub q3_qscore: u64,
    pub mean_expected_errors: f64,
}

/// Reads falling into one expected-error bucket.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ExpectedErrorBucket {
    pub label: String,
    pub reads: u64,
    pub reads_pct: f64,
}

/// Breakdown of homopolymer runs (>= 6 bp) and 3'-terminal tails (>= 10 bp).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct HomopolymerSummary {
    pub a_bases: u64,
    pub c_bases: u64,
    pub g_bases: u64,
    pub t_bases: u64,
    pub homopolymer_reads: u64,
    pub poly_a_tail_reads: u64,
    pub poly_c_tail_reads: u64,
    pub poly_g_tail_reads: u64,
    pub poly_t_tail_reads: u64,
    pub homopolymer_mean_qscore: f64,
    pub non_homopolymer_mean_qscore: f64,
    pub homopolymer_qscore_gap: f64,
    pub max_run_length: u64,
}

/// Sequencing metadata extracted from FASTQ headers.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct HeaderSummary {
    pub platform: Option<String>,
    pub chastity_failed_reads: u64,
    pub chastity_failed_pct: f64,
    pub run_group_count: u64,
    pub run_groups: Vec<String>,
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
    /// relative-position, length-quality, adapter, duplication, and
    /// overrepresented statistics.
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
        let sample = sample_name(path);
        let lengths = acc.lengths.counts();
        let positions = position_summaries(acc);
        let relative_positions = relative_position_summaries(acc);
        let length_quality = length_quality_summaries(acc);
        let gc_distribution: Vec<f64> = acc
            .gc
            .iter()
            .map(|w| *w as f64 / GC_WEIGHT as f64)
            .collect();
        let (gc_theoretical, gc_stdev) = normal_curve(&acc.gc, GC_WEIGHT);
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
        let mut metrics = metrics(acc, &lengths, &gc_distribution, &gc_theoretical, gc_stdev);
        let expected_error_distribution = expected_error_buckets(acc);
        let complexity_distribution = acc.complexity_hist.to_vec();
        let homopolymers = Some(homopolymer_summary(acc, &metrics));
        let header = Some(HeaderSummary {
            platform: metrics.platform.clone(),
            chastity_failed_reads: metrics.chastity_failed_reads,
            chastity_failed_pct: metrics.chastity_failed_pct,
            run_group_count: metrics.run_group_count,
            run_groups: acc.header.run_groups(),
        });
        if complete {
            let reads = acc.reads.max(1) as f64;
            metrics.duplicate_pct_est = duplication.as_ref().map(|d| d.duplicate_pct);
            metrics.adapter_reads_pct = Some(100.0 * acc.adapter_reads as f64 / reads);
            let poly_g: u64 = acc.adapters[POLY_G].iter().sum();
            metrics.poly_g_reads_pct = Some(100.0 * poly_g as f64 / reads);
            metrics.overrepresented_count = Some(overrepresented.len() as u64);
            metrics.qscore_tail_drop = qscore_tail_drop(&acc.relative_positions);
            metrics.effective_q20_span = Some(effective_quality_span(&positions, Q20 as u64));
            metrics.effective_q30_span = Some(effective_quality_span(&positions, Q30 as u64));
            metrics.head_bias_bp = Some(head_bias_bp(&positions));
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
            gc_distribution,
            gc_theoretical,
            positions,
            relative_positions,
            length_quality,
            expected_error_distribution,
            complexity_distribution,
            homopolymers,
            header,
            duplication,
            overrepresented,
            adapters,
        }
    }
}

/// Strip compound FASTQ extensions (`.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`)
/// from a file name, falling back to [`Path::file_stem`].
pub(crate) fn sample_name(path: &Path) -> String {
    let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
        return String::new();
    };
    let lower = name.to_ascii_lowercase();
    for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq"] {
        if lower.ends_with(suffix) && name.len() > suffix.len() {
            return name[..name.len() - suffix.len()].to_string();
        }
    }
    path.file_stem()
        .and_then(|n| n.to_str())
        .unwrap_or_default()
        .to_string()
}

fn metrics(
    acc: &ReadAccumulator,
    lengths: &[(usize, u64)],
    gc_obs: &[f64],
    gc_theo: &[f64],
    gc_stdev: f64,
) -> ReadMetrics {
    let base_count: u64 = acc.bases.iter().sum();
    let total = base_count as f64;
    let bases = &acc.bases;
    let gc = bases[BASE_G] + bases[BASE_C];
    let at = bases[BASE_A] + bases[BASE_T];
    let quality_count: u64 = acc.qualities.iter().sum();
    let quality_sum: u64 = weighted_sum(&acc.qualities);
    let max_qscore = last_nonzero(&acc.qualities).unwrap_or(0) as u64;
    let mean_error_rate = mean_error_rate_from_hist(&acc.qualities);
    let effective_mean_qscore = if quality_count == 0 {
        0.0
    } else {
        effective_qscore_from_error_rate(mean_error_rate, max_qscore)
    };
    let q10_bases: u64 = acc.qualities[Q10..].iter().sum();
    let q15_bases: u64 = acc.qualities[Q15..].iter().sum();
    let q20_bases: u64 = acc.qualities[Q20..].iter().sum();
    let q30_bases: u64 = acc.qualities[Q30..].iter().sum();
    let q40_bases: u64 = acc.qualities[Q40..].iter().sum();
    let distinct_qscores = acc.qualities.iter().filter(|&&c| c > 0).count() as u64;
    let scored_reads: u64 = acc.read_qualities.iter().sum();
    let mean_length = ratio(base_count as f64, acc.reads as f64);
    let total_ee = from_fixed(acc.expected_errors_sum);
    let ee_lt_1: u64 = acc.ee_hist[..4].iter().sum();
    let homopolymer_base_total: u64 = acc.homopolymer_bases.iter().sum();
    let homopolymer_mean_qscore =
        ratio(acc.homopolymer_q_sum as f64, acc.homopolymer_q_count as f64);
    let non_homopolymer_mean_qscore = ratio(
        acc.non_homopolymer_q_sum as f64,
        acc.non_homopolymer_q_count as f64,
    );
    let homopolymer_qscore_gap = if acc.homopolymer_q_count > 0 && acc.non_homopolymer_q_count > 0 {
        non_homopolymer_mean_qscore - homopolymer_mean_qscore
    } else {
        0.0
    };
    let empty_reads = lengths
        .first()
        .map_or(0, |(len, count)| if *len == 0 { *count } else { 0 });
    let non_empty_reads = acc.reads.saturating_sub(empty_reads) as f64;
    ReadMetrics {
        read_count: acc.reads,
        base_count,
        min_read_length: lengths.first().map(|l| l.0 as u64).unwrap_or(0),
        max_read_length: lengths.last().map(|l| l.0 as u64).unwrap_or(0),
        mean_read_length: mean_length,
        median_read_length: median_length(lengths, acc.reads),
        stdev_read_length: stdev_length(lengths, acc.reads, mean_length),
        n50_read_length: nx(lengths, base_count, 50),
        n90_read_length: nx(lengths, base_count, 90),
        l50_read_length: l50(lengths, base_count),
        aun_read_length: aun(lengths, base_count),
        empty_reads,
        g_count: bases[BASE_G],
        c_count: bases[BASE_C],
        a_count: bases[BASE_A],
        t_count: bases[BASE_T],
        n_count: bases[BASE_N],
        other_count: bases[BASE_OTHER],
        gc_content: ratio(gc as f64, total),
        at_content: ratio(at as f64, total),
        n_content: ratio(bases[BASE_N] as f64, total),
        at_skew: skew(bases[BASE_A], bases[BASE_T]),
        gc_skew: skew(bases[BASE_G], bases[BASE_C]),
        gc_stdev,
        gc_normal_divergence_pct: gc_normal_divergence_pct(gc_obs, gc_theo),
        gc_peaks: gc_peak_count(gc_obs),
        mean_qscore: ratio(quality_sum as f64, quality_count as f64),
        effective_mean_qscore,
        mean_error_rate,
        min_qscore: first_nonzero(&acc.qualities).unwrap_or(0) as u64,
        max_qscore,
        low_qscore_bases: acc.qualities[..LOW_QSCORE].iter().sum(),
        q20_bases,
        q30_bases,
        q10_pct: 100.0 * ratio(q10_bases as f64, quality_count as f64),
        q15_pct: 100.0 * ratio(q15_bases as f64, quality_count as f64),
        q20_pct: 100.0 * ratio(q20_bases as f64, quality_count as f64),
        q30_pct: 100.0 * ratio(q30_bases as f64, quality_count as f64),
        q40_pct: 100.0 * ratio(q40_bases as f64, quality_count as f64),
        distinct_qscores,
        mean_read_qscore: ratio(from_fixed(acc.read_qscore_sum), scored_reads as f64),
        low_quality_reads: acc.low_quality_reads,
        mean_expected_errors: ratio(total_ee, scored_reads as f64),
        expected_errors_per_100bp: 100.0 * ratio(total_ee, total),
        expected_error_free_reads_pct: 100.0
            * ratio(from_fixed(acc.error_free_prob_sum), scored_reads as f64),
        reads_ee_lt_1_pct: 100.0 * ratio(ee_lt_1 as f64, scored_reads as f64),
        q10_read_bases_pct: 100.0 * ratio(acc.q10_read_bases as f64, total),
        q20_read_bases_pct: 100.0 * ratio(acc.q20_read_bases as f64, total),
        q30_read_bases_pct: 100.0 * ratio(acc.q30_read_bases as f64, total),
        mean_complexity: ratio(from_fixed(acc.complexity_sum), non_empty_reads),
        low_complexity_reads_pct: 100.0 * ratio(acc.low_complexity_reads as f64, non_empty_reads),
        homopolymer_bases_pct: 100.0 * ratio(homopolymer_base_total as f64, total),
        max_homopolymer_length: acc.max_homopolymer_length,
        homopolymer_mean_qscore,
        non_homopolymer_mean_qscore,
        homopolymer_qscore_gap,
        poly_g_tail_reads_pct: 100.0
            * ratio(
                acc.terminal_homopolymer_reads[BASE_G] as f64,
                non_empty_reads,
            ),
        poly_a_tail_reads_pct: 100.0
            * ratio(
                acc.terminal_homopolymer_reads[BASE_A] as f64,
                non_empty_reads,
            ),
        poly_t_tail_reads_pct: 100.0
            * ratio(
                acc.terminal_homopolymer_reads[BASE_T] as f64,
                non_empty_reads,
            ),
        poly_c_tail_reads_pct: 100.0
            * ratio(
                acc.terminal_homopolymer_reads[BASE_C] as f64,
                non_empty_reads,
            ),
        platform: acc.header.platform().map(String::from),
        chastity_failed_reads: acc.header.chastity_failed,
        chastity_failed_pct: 100.0 * ratio(acc.header.chastity_failed as f64, acc.reads as f64),
        run_group_count: acc.header.run_group_count(),
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
            let width = (end - start + 1) as f64;
            let q = &counts.qualities;
            let n: u64 = q.iter().sum();
            let max_q = last_nonzero(q).unwrap_or(0) as u64;
            let effective_mean_qscore = if n == 0 {
                0.0
            } else {
                effective_qscore_from_error_rate(mean_error_rate_from_hist(q), max_q)
            };
            let avg_reads = n as f64 / width;
            PositionSummary {
                start: start as u64,
                end: end as u64,
                g_count: counts.bases[BASE_G],
                c_count: counts.bases[BASE_C],
                a_count: counts.bases[BASE_A],
                t_count: counts.bases[BASE_T],
                n_count: counts.bases[BASE_N],
                mean_qscore: ratio(weighted_sum(q) as f64, n as f64),
                effective_mean_qscore,
                read_survival_pct: (100.0 * ratio(avg_reads, acc.reads as f64)).min(100.0),
                min_qscore: first_nonzero(q).unwrap_or(0) as u64,
                max_qscore: max_q,
                p10_qscore: quantile(q, 0.10),
                q1_qscore: quantile(q, 0.25),
                median_qscore: quantile(q, 0.50),
                q3_qscore: quantile(q, 0.75),
                p90_qscore: quantile(q, 0.90),
            }
        })
        .collect()
}

fn relative_position_summaries(acc: &ReadAccumulator) -> Vec<RelativePositionSummary> {
    if acc.relative_positions.is_empty() {
        return Vec::new();
    }
    let step = 100 / RELATIVE_POSITION_BINS as u64;
    acc.relative_positions
        .iter()
        .enumerate()
        .map(|(i, counts)| {
            let q = &counts.qualities;
            let n: u64 = q.iter().sum();
            let max_q = last_nonzero(q).unwrap_or(0) as u64;
            let effective_mean_qscore = if n == 0 {
                0.0
            } else {
                effective_qscore_from_error_rate(mean_error_rate_from_hist(q), max_q)
            };
            RelativePositionSummary {
                start_pct: i as u64 * step,
                end_pct: (i as u64 + 1) * step,
                g_count: counts.bases[BASE_G],
                c_count: counts.bases[BASE_C],
                a_count: counts.bases[BASE_A],
                t_count: counts.bases[BASE_T],
                n_count: counts.bases[BASE_N],
                mean_qscore: ratio(weighted_sum(q) as f64, n as f64),
                effective_mean_qscore,
                p10_qscore: quantile(q, 0.10),
                q1_qscore: quantile(q, 0.25),
                median_qscore: quantile(q, 0.50),
                q3_qscore: quantile(q, 0.75),
                p90_qscore: quantile(q, 0.90),
            }
        })
        .collect()
}

fn length_quality_summaries(acc: &ReadAccumulator) -> Vec<LengthQualityBin> {
    acc.length_quality
        .iter()
        .zip(LENGTH_QUALITY_BOUNDS)
        .filter(|(counts, _)| counts.reads > 0)
        .map(|(counts, (label, min_length, max_length))| {
            let total_ee = from_fixed(counts.expected_errors_sum);
            let err_rate = ratio(total_ee, counts.bases as f64);
            let max_q = last_nonzero(&counts.read_qualities).unwrap_or(0) as u64;
            LengthQualityBin {
                label: label.to_string(),
                min_length,
                max_length,
                reads: counts.reads,
                bases: counts.bases,
                mean_qscore: ratio(counts.base_quality_sum as f64, counts.bases as f64),
                effective_mean_qscore: effective_qscore_from_error_rate(err_rate, max_q),
                q1_qscore: quantile(&counts.read_qualities, 0.25),
                median_qscore: quantile(&counts.read_qualities, 0.50),
                q3_qscore: quantile(&counts.read_qualities, 0.75),
                mean_expected_errors: ratio(total_ee, counts.reads as f64),
            }
        })
        .collect()
}

fn expected_error_buckets(acc: &ReadAccumulator) -> Vec<ExpectedErrorBucket> {
    let total: u64 = acc.ee_hist.iter().sum();
    EE_BOUNDS
        .iter()
        .zip(acc.ee_hist)
        .map(|((_, label), reads)| ExpectedErrorBucket {
            label: (*label).to_string(),
            reads,
            reads_pct: 100.0 * ratio(reads as f64, total as f64),
        })
        .collect()
}

fn homopolymer_summary(acc: &ReadAccumulator, metrics: &ReadMetrics) -> HomopolymerSummary {
    HomopolymerSummary {
        a_bases: acc.homopolymer_bases[BASE_A],
        c_bases: acc.homopolymer_bases[BASE_C],
        g_bases: acc.homopolymer_bases[BASE_G],
        t_bases: acc.homopolymer_bases[BASE_T],
        homopolymer_reads: acc.homopolymer_reads,
        poly_a_tail_reads: acc.terminal_homopolymer_reads[BASE_A],
        poly_c_tail_reads: acc.terminal_homopolymer_reads[BASE_C],
        poly_g_tail_reads: acc.terminal_homopolymer_reads[BASE_G],
        poly_t_tail_reads: acc.terminal_homopolymer_reads[BASE_T],
        homopolymer_mean_qscore: metrics.homopolymer_mean_qscore,
        non_homopolymer_mean_qscore: metrics.non_homopolymer_mean_qscore,
        homopolymer_qscore_gap: metrics.homopolymer_qscore_gap,
        max_run_length: acc.max_homopolymer_length,
    }
}

/// First-quarter minus last-quarter effective mean quality across relative read length.
fn qscore_tail_drop(relative: &[PositionCounts]) -> Option<f64> {
    if relative.len() < 4 {
        return None;
    }
    let quarter = relative.len() / 4;
    let first = merge_quality_hists(&relative[..quarter]);
    let last = merge_quality_hists(&relative[relative.len() - quarter..]);
    let first_n: u64 = first.iter().sum();
    let last_n: u64 = last.iter().sum();
    if first_n == 0 || last_n == 0 {
        return None;
    }
    let first_q = effective_qscore_from_error_rate(
        mean_error_rate_from_hist(&first),
        last_nonzero(&first).unwrap_or(0) as u64,
    );
    let last_q = effective_qscore_from_error_rate(
        mean_error_rate_from_hist(&last),
        last_nonzero(&last).unwrap_or(0) as u64,
    );
    Some(first_q - last_q)
}

fn merge_quality_hists(bins: &[PositionCounts]) -> [u64; QSCORE_BINS] {
    let mut merged = [0u64; QSCORE_BINS];
    for bin in bins {
        for (dst, src) in merged.iter_mut().zip(&bin.qualities) {
            *dst += *src;
        }
    }
    merged
}

/// Longest contiguous prefix (bp) where per-position median quality stays >= `threshold`.
fn effective_quality_span(positions: &[PositionSummary], threshold: u64) -> u64 {
    let mut span = 0;
    for pos in positions {
        if pos.median_qscore >= threshold {
            span = pos.end;
        } else {
            break;
        }
    }
    span
}

/// Leading positions (bp) with base composition differing > 5% from the read body.
fn head_bias_bp(positions: &[PositionSummary]) -> u64 {
    if positions.len() < 10 {
        return 0;
    }
    let body_start = positions.len() / 4;
    let body_end = (positions.len() * 3) / 4;
    let mut body = [0u64; 4];
    for pos in &positions[body_start..body_end] {
        body[0] += pos.a_count;
        body[1] += pos.c_count;
        body[2] += pos.g_count;
        body[3] += pos.t_count;
    }
    let body_total: u64 = body.iter().sum();
    if body_total == 0 {
        return 0;
    }
    let body_pct = body.map(|c| 100.0 * c as f64 / body_total as f64);
    let window = HEAD_BIAS_WINDOW.min(positions.len() / 3);
    let mut last_biased = 0u64;
    for pos in &positions[..window] {
        let counts = [pos.a_count, pos.c_count, pos.g_count, pos.t_count];
        let total: u64 = counts.iter().sum();
        if total == 0 {
            continue;
        }
        let max_diff = counts
            .iter()
            .zip(body_pct)
            .map(|(&c, expected)| (100.0 * c as f64 / total as f64 - expected).abs())
            .fold(0.0f64, f64::max);
        if max_diff > HEAD_BIAS_THRESHOLD_PCT {
            last_biased = pos.end;
        }
    }
    last_biased
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

/// Normal curve fitted to the GC histogram, scaled to the same number of reads,
/// along with the histogram standard deviation.
///
/// `weight` is the histogram count of one read.
fn normal_curve(gc: &[u64], weight: u64) -> (Vec<f64>, f64) {
    let total: u64 = gc.iter().sum();
    if total == 0 {
        return (vec![0.0; gc.len()], 0.0);
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
        return (gc.iter().map(|c| *c as f64 / weight as f64).collect(), 0.0);
    }
    let scale = reads / (sd * (2.0 * std::f64::consts::PI).sqrt());
    let curve = (0..gc.len())
        .map(|x| scale * (-(x as f64 - mean).powi(2) / (2.0 * variance)).exp())
        .collect();
    (curve, sd)
}

/// Total-variation divergence between observed and normal GC curves, in percent.
fn gc_normal_divergence_pct(obs: &[f64], theo: &[f64]) -> f64 {
    let obs_sum: f64 = obs.iter().sum();
    let theo_sum: f64 = theo.iter().sum();
    if obs_sum <= 0.0 || theo_sum <= 0.0 {
        return 0.0;
    }
    let l1: f64 = obs
        .iter()
        .zip(theo)
        .map(|(o, t)| (o / obs_sum - t / theo_sum).abs())
        .sum();
    50.0 * l1
}

/// Prominent local maxima in a smoothed per-read GC distribution.
fn gc_peak_count(obs: &[f64]) -> u64 {
    if obs.is_empty() {
        return 0;
    }
    let n = obs.len();
    let mut smooth = vec![0.0; n];
    for i in 0..n {
        let prev = if i > 0 { obs[i - 1] } else { obs[i] };
        let next = if i + 1 < n { obs[i + 1] } else { obs[i] };
        smooth[i] = 0.25 * prev + 0.5 * obs[i] + 0.25 * next;
    }
    let global_max = smooth.iter().copied().fold(0.0f64, f64::max);
    if global_max <= 0.0 {
        return 0;
    }
    let min_height = global_max * 0.05;
    let mut candidates = Vec::new();
    let mut i = 0;
    while i < n {
        // Collapse flat plateaus into their midpoint.
        let mut end = i;
        while end + 1 < n && (smooth[end + 1] - smooth[i]).abs() < 1e-12 {
            end += 1;
        }
        let left_lower = i == 0 || smooth[i - 1] < smooth[i];
        let right_lower = end + 1 == n || smooth[end + 1] < smooth[i];
        if left_lower && right_lower && smooth[i] >= min_height {
            candidates.push(((i + end) / 2, smooth[i]));
        }
        i = end + 1;
    }
    // Keep peaks separated by a real valley (at least 20% drop below the shorter peak).
    let mut kept: Vec<(usize, f64)> = Vec::new();
    for (idx, height) in candidates {
        if let Some((prev_idx, prev_height)) = kept.last_mut() {
            let valley = smooth[*prev_idx..=idx]
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min);
            let shorter = prev_height.min(height);
            if valley > shorter * 0.80 {
                if height > *prev_height {
                    *prev_idx = idx;
                    *prev_height = height;
                }
                continue;
            }
        }
        kept.push((idx, height));
    }
    kept.len() as u64
}

fn mean_error_rate_from_hist(hist: &[u64; QSCORE_BINS]) -> f64 {
    let total: u64 = hist.iter().sum();
    if total == 0 {
        return 0.0;
    }
    let errors: f64 = hist
        .iter()
        .zip(ERROR_PROB.iter())
        .map(|(&count, &p)| count as f64 * p)
        .sum();
    errors / total as f64
}

fn effective_qscore_from_error_rate(error_rate: f64, max_q: u64) -> f64 {
    if error_rate <= 0.0 {
        max_q as f64
    } else {
        (-10.0 * error_rate.log10()).clamp(0.0, MAX_QSCORE as f64)
    }
}

fn skew(a: u64, b: u64) -> f64 {
    let total = a + b;
    if total == 0 {
        0.0
    } else {
        (a as f64 - b as f64) / total as f64
    }
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

/// Shortest read length that covers at least `pct`% of bases with reads at least that long.
fn nx(lengths: &[(usize, u64)], bases: u64, pct: u64) -> u64 {
    if bases == 0 {
        return 0;
    }
    let target = (u128::from(bases) * u128::from(pct)).div_ceil(100);
    let mut covered = 0u128;
    for (len, count) in lengths.iter().rev() {
        covered += *len as u128 * u128::from(*count);
        if covered >= target {
            return *len as u64;
        }
    }
    0
}

/// Minimum number of longest reads covering half the bases.
fn l50(lengths: &[(usize, u64)], bases: u64) -> u64 {
    if bases == 0 {
        return 0;
    }
    let half = u128::from(bases).div_ceil(2);
    let mut covered = 0u128;
    let mut reads = 0u64;
    for &(len, count) in lengths.iter().rev() {
        let len = len as u128;
        if len == 0 {
            continue;
        }
        let group = len * u128::from(count);
        if covered + group >= half {
            let needed = (half - covered).div_ceil(len) as u64;
            return reads + needed;
        }
        covered += group;
        reads += count;
    }
    reads
}

/// Area under the Nx curve (`sum(len^2) / sum(len)`).
fn aun(lengths: &[(usize, u64)], bases: u64) -> f64 {
    if bases == 0 {
        return 0.0;
    }
    let sum_sq: f64 = lengths
        .iter()
        .map(|&(len, count)| (len as f64).powi(2) * count as f64)
        .sum();
    sum_sq / bases as f64
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
    fn test_sample_name_strips_fastq_extensions() {
        assert_eq!(
            sample_name(Path::new("dir/sample_R1.fastq.gz")),
            "sample_R1"
        );
        assert_eq!(sample_name(Path::new("sample_R2.FQ.GZ")), "sample_R2");
        assert_eq!(sample_name(Path::new("read_1.fastq")), "read_1");
        assert_eq!(sample_name(Path::new("read_2.fq")), "read_2");
        assert_eq!(sample_name(Path::new("custom.txt")), "custom");
    }

    #[test]
    fn test_length_statistics() {
        let lengths = vec![(2, 1), (4, 2), (10, 1)];
        assert_eq!(median_length(&lengths, 4), 4.0);
        assert_eq!(median_length(&[(3, 1), (5, 1)], 2), 4.0);
        // 20 bases in total; the 10 bp read covers half.
        assert_eq!(nx(&lengths, 20, 50), 10);
        assert_eq!(nx(&lengths, 20, 90), 4);
        assert_eq!(l50(&lengths, 20), 1);
        assert!((aun(&lengths, 20) - 6.8).abs() < 1e-12);
        assert_eq!(nx(&[(5, 4)], 20, 50), 5);
        assert_eq!(l50(&[(5, 4)], 20), 2);
        let mean = 5.0;
        let expected = ((9.0 + 1.0 + 1.0 + 25.0) / 3.0f64).sqrt();
        assert!((stdev_length(&lengths, 4, mean) - expected).abs() < 1e-12);
    }

    #[test]
    fn test_effective_qscore_exposes_bad_bases() {
        let mut hist = [0u64; QSCORE_BINS];
        hist[40] = 99;
        hist[2] = 1;
        let err = mean_error_rate_from_hist(&hist);
        let eff = effective_qscore_from_error_rate(err, 40);
        // Arithmetic mean is 39.62, while error-rate effective Q is ~21.93.
        assert!((eff - 21.93).abs() < 0.05, "eff={eff}");
    }

    #[test]
    fn test_gc_peaks_unimodal_vs_bimodal() {
        let mut unimodal = vec![0.0; 101];
        unimodal[50] = 100.0;
        unimodal[49] = 80.0;
        unimodal[51] = 80.0;
        assert_eq!(gc_peak_count(&unimodal), 1);

        let mut bimodal = vec![0.0; 101];
        bimodal[35] = 100.0;
        bimodal[34] = 60.0;
        bimodal[36] = 60.0;
        bimodal[65] = 90.0;
        bimodal[64] = 50.0;
        bimodal[66] = 50.0;
        assert_eq!(gc_peak_count(&bimodal), 2);
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
        let (curve, sd) = normal_curve(&gc, 1);
        let area: f64 = curve.iter().sum();
        assert!((area - 40.0).abs() < 0.5);
        assert!(curve[50] > curve[40]);
        assert!(sd > 0.0);
    }
}
