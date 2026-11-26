//! Paired-end integrity and cohort-level read summary diagnostics.
//!
//! Computed from per-sample [`SampleReadReport`] entries after single-pass
//! summarization, requiring no extra file I/O.
use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use super::report::SampleReadReport;

/// Minimum share of cohort median reads before a non-empty sample is flagged as low yield.
const LOW_YIELD_FRACTION: f64 = 0.10;
/// Modified Z-score threshold (Iglewicz & Hoaglin) for cohort-relative outliers.
const OUTLIER_MODIFIED_Z: f64 = 3.5;
/// Minimum non-empty samples needed to compute cohort-relative outliers.
const MIN_OUTLIER_SAMPLES: usize = 4;

/// Paired-end consistency and quality comparison between R1 and R2 mates.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct PairedReadSummary {
    /// Shared sample identifier with the mate marker removed.
    pub pair_name: String,
    pub r1_sample: String,
    pub r2_sample: String,
    pub r1_reads: u64,
    pub r2_reads: u64,
    /// `r1_reads - r2_reads` as a signed count.
    pub read_count_diff: i64,
    /// True when R1 and R2 have the exact same read count.
    pub in_sync: bool,
    /// `r1_q30_pct - r2_q30_pct` (positive means R2 has lower Q30%).
    pub q30_drop_pct: f64,
    /// `r1_effective_q - r2_effective_q` (positive means R2 has lower effective Q).
    pub effective_q_drop: f64,
    /// `r1_mean_length - r2_mean_length`.
    pub mean_length_diff: f64,
    /// `100 * (r1_gc - r2_gc)` in percentage points.
    pub gc_diff_pct: f64,
}

/// Cohort-level pooling balance and outlier summary across all samples in a run.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct CohortSummary {
    pub sample_count: u64,
    pub empty_sample_count: u64,
    pub total_reads: u64,
    pub total_bases: u64,
    pub median_reads: f64,
    pub p10_reads: f64,
    pub p90_reads: f64,
    /// Ratio of P90 to P10 read count across non-empty samples (0 when P10 is 0).
    pub depth_p90_p10_ratio: f64,
    /// Coefficient of variation (`stdev / mean`) of read counts across non-empty samples.
    pub depth_cv: f64,
    /// Gini coefficient (`0.0` = even pooling, `1.0` = concentrated in one sample).
    pub depth_gini: f64,
    /// Non-empty samples with fewer than 10% of the cohort median read count.
    pub low_yield_samples: Vec<String>,
    /// Samples strongly deviating from the cohort median on key quality or composition metrics.
    pub outliers: Vec<CohortOutlier>,
}

/// One sample-metric pair flagged as a cohort-relative outlier.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct CohortOutlier {
    pub sample: String,
    pub metric: String,
    pub value: f64,
    pub cohort_median: f64,
    pub modified_z: f64,
}

/// Pair R1 and R2 samples by standard naming conventions (`_R1`/`_R2`,
/// `_R1_001`/`_R2_001`, `_1`/`_2`, `.R1`/`.R2`, `.1`/`.2`, `-R1`/`-R2`, `-1`/`-2`).
pub fn pair_samples(samples: &[SampleReadReport]) -> Vec<PairedReadSummary> {
    let mut grouped: BTreeMap<String, (Vec<&SampleReadReport>, Vec<&SampleReadReport>)> =
        BTreeMap::new();
    for sample in samples {
        if let Some((pair_name, mate)) = split_mate_name(&sample.sample) {
            let entry = grouped.entry(pair_name).or_default();
            if mate == 1 {
                entry.0.push(sample);
            } else {
                entry.1.push(sample);
            }
        }
    }
    let mut pairs = Vec::new();
    for (pair_name, (r1s, r2s)) in grouped {
        // Only pair unambiguous 1-to-1 matches.
        if r1s.len() != 1 || r2s.len() != 1 {
            continue;
        }
        let r1 = r1s[0];
        let r2 = r2s[0];
        let m1 = &r1.metrics;
        let m2 = &r2.metrics;
        let read_count_diff = m1.read_count as i64 - m2.read_count as i64;
        pairs.push(PairedReadSummary {
            pair_name,
            r1_sample: r1.sample.clone(),
            r2_sample: r2.sample.clone(),
            r1_reads: m1.read_count,
            r2_reads: m2.read_count,
            read_count_diff,
            in_sync: m1.read_count == m2.read_count,
            q30_drop_pct: m1.q30_pct - m2.q30_pct,
            effective_q_drop: m1.effective_mean_qscore - m2.effective_mean_qscore,
            mean_length_diff: m1.mean_read_length - m2.mean_read_length,
            gc_diff_pct: 100.0 * (m1.gc_content - m2.gc_content),
        });
    }
    pairs
}

/// Summarize sequencing depth balance and robust cohort outliers across `samples`.
pub fn summarize_cohort(samples: &[SampleReadReport]) -> Option<CohortSummary> {
    if samples.is_empty() {
        return None;
    }
    let non_empty: Vec<&SampleReadReport> = samples
        .iter()
        .filter(|s| s.metrics.read_count > 0)
        .collect();
    let sample_count = samples.len() as u64;
    let empty_sample_count = sample_count - non_empty.len() as u64;
    let total_reads: u64 = non_empty.iter().map(|s| s.metrics.read_count).sum();
    let total_bases: u64 = non_empty.iter().map(|s| s.metrics.base_count).sum();
    if non_empty.is_empty() {
        return Some(CohortSummary {
            sample_count,
            empty_sample_count,
            ..Default::default()
        });
    }
    let mut read_counts: Vec<f64> = non_empty
        .iter()
        .map(|s| s.metrics.read_count as f64)
        .collect();
    read_counts.sort_by(f64::total_cmp);
    let median_reads = percentile(&read_counts, 0.50);
    let p10_reads = percentile(&read_counts, 0.10);
    let p90_reads = percentile(&read_counts, 0.90);
    let depth_p90_p10_ratio = if p10_reads > 0.0 {
        p90_reads / p10_reads
    } else {
        0.0
    };
    let mean_reads = total_reads as f64 / non_empty.len() as f64;
    let depth_cv = if non_empty.len() > 1 && mean_reads > 0.0 {
        let var = read_counts
            .iter()
            .map(|&r| (r - mean_reads).powi(2))
            .sum::<f64>()
            / (non_empty.len() - 1) as f64;
        var.sqrt() / mean_reads
    } else {
        0.0
    };
    let depth_gini = gini(&read_counts);
    let low_threshold = median_reads * LOW_YIELD_FRACTION;
    let low_yield_samples: Vec<String> = if non_empty.len() >= 2 && low_threshold > 0.0 {
        non_empty
            .iter()
            .filter(|s| (s.metrics.read_count as f64) < low_threshold)
            .map(|s| s.sample.clone())
            .collect()
    } else {
        Vec::new()
    };
    let outliers = detect_outliers(&non_empty);
    Some(CohortSummary {
        sample_count,
        empty_sample_count,
        total_reads,
        total_bases,
        median_reads,
        p10_reads,
        p90_reads,
        depth_p90_p10_ratio,
        depth_cv,
        depth_gini,
        low_yield_samples,
        outliers,
    })
}

/// Extract `(pair_name, mate)` where `mate` is `1` or `2`, if `name` follows a
/// standard paired-end naming pattern.
fn split_mate_name(name: &str) -> Option<(String, u8)> {
    // 1. Illumina-style `_R1_001` / `_R2_001` (or `.R1_001`, `-R1_001`).
    if let Some((prefix, tail)) = name.rsplit_once('_')
        && !tail.is_empty()
        && tail.bytes().all(|b| b.is_ascii_digit())
        && let Some((base, mate)) = strip_r_marker(prefix)
    {
        return Some((format!("{base}_{tail}"), mate));
    }
    // 2. Trailing `_R1`/`_R2`, `.R1`/`.R2`, `-R1`/`-R2`.
    if let Some(pair) = strip_r_marker(name) {
        return Some(pair);
    }
    // 3. Trailing `_1`/`_2`, `.1`/`.2`, `-1`/`-2`.
    for sep in ['_', '.', '-'] {
        for (suffix, mate) in [('1', 1u8), ('2', 2u8)] {
            let marker = format!("{sep}{suffix}");
            if let Some(base) = name.strip_suffix(&marker)
                && !base.is_empty()
            {
                return Some((base.to_string(), mate));
            }
        }
    }
    None
}

fn strip_r_marker(name: &str) -> Option<(String, u8)> {
    let bytes = name.as_bytes();
    if bytes.len() < 4 {
        return None;
    }
    let sep = bytes[bytes.len() - 3];
    let r = bytes[bytes.len() - 2];
    let num = bytes[bytes.len() - 1];
    if matches!(sep, b'_' | b'.' | b'-') && matches!(r, b'R' | b'r') && matches!(num, b'1' | b'2') {
        let base = &name[..bytes.len() - 3];
        if !base.is_empty() {
            return Some((base.to_string(), if num == b'1' { 1 } else { 2 }));
        }
    }
    None
}

type SampleMetricFn = fn(&SampleReadReport) -> Option<f64>;

fn detect_outliers(samples: &[&SampleReadReport]) -> Vec<CohortOutlier> {
    if samples.len() < MIN_OUTLIER_SAMPLES {
        return Vec::new();
    }
    let specs: [(&str, f64, SampleMetricFn); 5] = [
        ("GC (%)", 3.0, |s| Some(100.0 * s.metrics.gc_content)),
        ("Effective Q", 3.0, |s| {
            let q = s.metrics.effective_mean_qscore;
            (q > 0.0).then_some(q)
        }),
        ("Mean length (bp)", 10.0, |s| {
            Some(s.metrics.mean_read_length)
        }),
        ("Low complexity (%)", 5.0, |s| {
            Some(s.metrics.low_complexity_reads_pct)
        }),
        ("Duplication (%)", 10.0, |s| s.metrics.duplicate_pct_est),
    ];
    let mut found = Vec::new();
    for (metric, min_delta, extract) in specs {
        let pairs: Vec<(&str, f64)> = samples
            .iter()
            .filter_map(|s| extract(s).map(|v| (s.sample.as_str(), v)))
            .collect();
        if pairs.len() < MIN_OUTLIER_SAMPLES {
            continue;
        }
        let mut sorted: Vec<f64> = pairs.iter().map(|&(_, v)| v).collect();
        sorted.sort_by(f64::total_cmp);
        let median = percentile(&sorted, 0.50);
        let mut abs_devs: Vec<f64> = sorted.iter().map(|&v| (v - median).abs()).collect();
        abs_devs.sort_by(f64::total_cmp);
        let mad = percentile(&abs_devs, 0.50);
        if mad > 1e-9 {
            for (sample, value) in pairs {
                let delta = (value - median).abs();
                let z = 0.6745 * (value - median) / mad;
                if z.abs() >= OUTLIER_MODIFIED_Z && delta >= min_delta {
                    found.push(CohortOutlier {
                        sample: sample.to_string(),
                        metric: metric.to_string(),
                        value,
                        cohort_median: median,
                        modified_z: z,
                    });
                }
            }
        } else {
            // Flat majority (MAD == 0): flag samples well beyond a meaningful floor.
            for (sample, value) in pairs {
                let delta = value - median;
                if delta.abs() >= min_delta * 2.0 {
                    found.push(CohortOutlier {
                        sample: sample.to_string(),
                        metric: metric.to_string(),
                        value,
                        cohort_median: median,
                        modified_z: (delta / min_delta) * OUTLIER_MODIFIED_Z,
                    });
                }
            }
        }
    }
    found.sort_by(|a, b| {
        b.modified_z
            .abs()
            .total_cmp(&a.modified_z.abs())
            .then_with(|| a.sample.cmp(&b.sample))
    });
    found
}

/// Gini coefficient of sorted non-negative values in `[0.0, 1.0]`.
fn gini(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    if n <= 1 {
        return 0.0;
    }
    let sum: f64 = sorted.iter().sum();
    if sum <= 0.0 {
        return 0.0;
    }
    let weighted: f64 = sorted
        .iter()
        .enumerate()
        .map(|(i, &x)| (2.0 * (i + 1) as f64 - n as f64 - 1.0) * x)
        .sum();
    (weighted / (n as f64 * sum)).clamp(0.0, 1.0)
}

/// Linearly interpolated percentile of sorted values.
fn percentile(sorted: &[f64], level: f64) -> f64 {
    let rank = (sorted.len() - 1) as f64 * level;
    let low = rank.floor() as usize;
    let high = rank.ceil() as usize;
    sorted[low] + (sorted[high] - sorted[low]) * (rank - low as f64)
}

#[cfg(test)]
mod test {
    use super::*;

    fn sample(name: &str, reads: u64, gc: f64, q30: f64, eff_q: f64) -> SampleReadReport {
        let mut s = SampleReadReport {
            sample: name.to_string(),
            ..Default::default()
        };
        s.metrics.read_count = reads;
        s.metrics.base_count = reads * 150;
        s.metrics.mean_read_length = 150.0;
        s.metrics.gc_content = gc;
        s.metrics.q30_pct = q30;
        s.metrics.effective_mean_qscore = eff_q;
        s
    }

    #[test]
    fn test_pair_samples_conventions_and_sync() {
        let samples = vec![
            sample("taxA_S1_L001_R1_001", 1000, 0.45, 92.0, 31.0),
            sample("taxA_S1_L001_R2_001", 1000, 0.46, 84.0, 27.0),
            sample("taxB_1", 500, 0.40, 90.0, 30.0),
            sample("taxB_2", 495, 0.40, 88.0, 29.0),
            sample("single_end", 800, 0.45, 91.0, 30.0),
        ];
        let pairs = pair_samples(&samples);
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].pair_name, "taxA_S1_L001_001");
        assert!(pairs[0].in_sync);
        assert!((pairs[0].q30_drop_pct - 8.0).abs() < 1e-9);
        assert!((pairs[0].effective_q_drop - 4.0).abs() < 1e-9);
        assert_eq!(pairs[1].pair_name, "taxB");
        assert!(!pairs[1].in_sync);
        assert_eq!(pairs[1].read_count_diff, 5);
    }

    #[test]
    fn test_cohort_balance_and_outliers() {
        let samples = vec![
            sample("s1", 1_000, 0.42, 92.0, 31.0),
            sample("s2", 1_050, 0.43, 91.0, 30.5),
            sample("s3", 980, 0.41, 93.0, 31.5),
            sample("s4", 1_020, 0.42, 92.0, 31.0),
            sample("dropout", 50, 0.68, 70.0, 18.0),
        ];
        let cohort = summarize_cohort(&samples).unwrap();
        assert_eq!(cohort.sample_count, 5);
        assert_eq!(cohort.empty_sample_count, 0);
        assert_eq!(cohort.total_reads, 4_100);
        assert_eq!(cohort.low_yield_samples, vec!["dropout"]);
        assert!(cohort.depth_gini > 0.1);
        assert!(
            cohort
                .outliers
                .iter()
                .any(|o| o.sample == "dropout" && o.metric == "GC (%)")
        );
        assert!(
            cohort
                .outliers
                .iter()
                .any(|o| o.sample == "dropout" && o.metric == "Effective Q")
        );
    }
}
