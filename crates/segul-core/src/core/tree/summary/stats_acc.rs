//! Online statistics, exact Chen-Shao HPD intervals, and BEAST/NHX annotation accumulation.

use ahash::AHashMap;
use indexmap::IndexSet;

use crate::parser::tree::{Annotation, AnnotationFormat, AnnotationValue};

/// Welford online summary statistics (count, mean, sample SD, min, max).
#[derive(Debug, Clone, Default)]
pub(crate) struct OnlineStats {
    pub(crate) count: u64,
    pub(crate) mean: f64,
    pub(crate) m2: f64,
    pub(crate) min: f64,
    pub(crate) max: f64,
}

impl OnlineStats {
    #[inline]
    pub(crate) fn push(&mut self, value: f64) {
        if !value.is_finite() {
            return;
        }
        if self.count == 0 {
            self.count = 1;
            self.mean = value;
            self.m2 = 0.0;
            self.min = value;
            self.max = value;
            return;
        }
        self.count += 1;
        let delta = value - self.mean;
        self.mean += delta / (self.count as f64);
        let delta2 = value - self.mean;
        self.m2 += delta * delta2;
        if value < self.min {
            self.min = value;
        }
        if value > self.max {
            self.max = value;
        }
    }

    pub(crate) fn merge(&mut self, other: &Self) {
        if other.count == 0 {
            return;
        }
        if self.count == 0 {
            *self = other.clone();
            return;
        }
        let n_a = self.count as f64;
        let n_b = other.count as f64;
        let n_total = n_a + n_b;
        let delta = other.mean - self.mean;
        self.mean += delta * (n_b / n_total);
        self.m2 += other.m2 + delta * delta * (n_a * n_b / n_total);
        self.count += other.count;
        if other.min < self.min {
            self.min = other.min;
        }
        if other.max > self.max {
            self.max = other.max;
        }
    }
}

/// Finalized summary statistics for a continuous node/branch metric.
#[derive(Debug, Clone, Default)]
pub(crate) struct MetricSummary {
    pub(crate) count: u64,
    pub(crate) mean: Option<f64>,
    pub(crate) median: Option<f64>,
    pub(crate) min: Option<f64>,
    pub(crate) max: Option<f64>,
    pub(crate) hpd_low: Option<f64>,
    pub(crate) hpd_high: Option<f64>,
}

/// Accumulator combining O(1)-memory Welford statistics with optional sample storage.
#[derive(Debug, Clone)]
pub(crate) struct MetricAccumulator {
    online: OnlineStats,
    samples: Option<Vec<f64>>,
}

impl MetricAccumulator {
    pub(crate) fn new(store_samples: bool) -> Self {
        Self {
            online: OnlineStats::default(),
            samples: if store_samples {
                Some(Vec::new())
            } else {
                None
            },
        }
    }

    #[inline]
    pub(crate) fn push(&mut self, value: f64) {
        if !value.is_finite() {
            return;
        }
        self.online.push(value);
        if let Some(samples) = &mut self.samples {
            samples.push(value);
        }
    }

    pub(crate) fn merge(&mut self, other: Self) {
        self.online.merge(&other.online);
        if let (Some(dst), Some(src)) = (&mut self.samples, other.samples) {
            dst.extend(src);
        }
    }

    #[cfg(test)]
    #[inline]
    pub(crate) fn mean(&self) -> Option<f64> {
        if self.online.count > 0 {
            Some(self.online.mean)
        } else {
            None
        }
    }

    /// Finalize statistics, sorting any retained samples for exact median and Chen-Shao HPD.
    pub(crate) fn finalize(&mut self, hpd_level: f64) -> MetricSummary {
        if self.online.count == 0 {
            return MetricSummary::default();
        }
        let (median, hpd) = if let Some(samples) = &mut self.samples {
            samples.sort_unstable_by(f64::total_cmp);
            (compute_median(samples), compute_hpd(samples, hpd_level))
        } else {
            (None, None)
        };
        MetricSummary {
            count: self.online.count,
            mean: Some(self.online.mean),
            median,
            min: Some(self.online.min),
            max: Some(self.online.max),
            hpd_low: hpd.map(|(lo, _)| lo),
            hpd_high: hpd.map(|(_, hi)| hi),
        }
    }
}

/// Compute exact sample median from a sorted finite slice.
pub(crate) fn compute_median(sorted: &[f64]) -> Option<f64> {
    let n = sorted.len();
    if n == 0 {
        return None;
    }
    if n % 2 == 1 {
        Some(sorted[n / 2])
    } else {
        Some(0.5 * (sorted[n / 2 - 1] + sorted[n / 2]))
    }
}

/// Compute the shortest empirical Highest Posterior Density (Chen-Shao) interval.
pub(crate) fn compute_hpd(sorted: &[f64], level: f64) -> Option<(f64, f64)> {
    if !(level > 0.0 && level < 1.0) || sorted.len() < 2 {
        return None;
    }
    let n = sorted.len();
    let span = ((n as f64) * level).floor() as usize;
    let span = span.clamp(1, n - 1);

    let mut best_idx = 0usize;
    let mut best_width = f64::INFINITY;
    for i in 0..=(n - 1 - span) {
        let width = sorted[i + span] - sorted[i];
        if width < best_width {
            best_width = width;
            best_idx = i;
        }
    }
    Some((sorted[best_idx], sorted[best_idx + span]))
}

/// Interned annotation key table so per-clade maps store compact `usize` keys.
#[derive(Debug, Clone, Default)]
pub(crate) struct AnnotationKeyTable {
    keys: IndexSet<String>,
}

impl AnnotationKeyTable {
    pub(crate) fn intern(&mut self, key: &str) -> Option<usize> {
        let trimmed = key.trim();
        if trimmed.is_empty() || is_synthetic_summary_key(trimmed) {
            return None;
        }
        if let Some(idx) = self.keys.get_index_of(trimmed) {
            Some(idx)
        } else {
            let (idx, _) = self.keys.insert_full(trimmed.to_string());
            Some(idx)
        }
    }

    pub(crate) fn get(&self, idx: usize) -> Option<&str> {
        self.keys.get_index(idx).map(String::as_str)
    }
}

fn is_synthetic_summary_key(key: &str) -> bool {
    let lower = key.to_ascii_lowercase();
    matches!(
        lower.as_str(),
        "posterior"
            | "prob"
            | "support"
            | "height"
            | "height_median"
            | "height_range"
            | "length"
            | "length_median"
            | "length_range"
            | "ca_height"
            | "decisive"
            | "supporting"
    ) || lower.ends_with("_hpd")
        || lower.ends_with("_median")
        || lower.ends_with("_range")
        || lower.ends_with(".prob")
}

/// Per-key accumulator supporting both numeric and categorical annotation tags.
#[derive(Debug, Clone)]
pub(crate) enum AnnotationValueAcc {
    Numeric(MetricAccumulator),
    Categorical(AHashMap<String, u64>),
}

/// Sparse per-clade collection of summarized BEAST / FigTree / NHX annotations.
#[derive(Debug, Clone, Default)]
pub(crate) struct CladeAnnotationsAcc {
    entries: AHashMap<usize, AnnotationValueAcc>,
}

impl CladeAnnotationsAcc {
    pub(crate) fn record_annotation(
        &mut self,
        annot: &Annotation,
        keys: &mut AnnotationKeyTable,
        store_samples: bool,
    ) {
        if !matches!(
            annot.format,
            AnnotationFormat::Beast | AnnotationFormat::Nhx
        ) {
            return;
        }
        for (raw_key, val) in &annot.entries {
            let Some(key_id) = keys.intern(raw_key) else {
                continue;
            };
            match extract_numeric_or_cat(val) {
                ParsedTag::Numeric(num) => match self.entries.get_mut(&key_id) {
                    Some(AnnotationValueAcc::Numeric(acc)) => acc.push(num),
                    Some(AnnotationValueAcc::Categorical(_)) => {}
                    None => {
                        let mut acc = MetricAccumulator::new(store_samples);
                        acc.push(num);
                        self.entries
                            .insert(key_id, AnnotationValueAcc::Numeric(acc));
                    }
                },
                ParsedTag::Categorical(label) => match self.entries.get_mut(&key_id) {
                    Some(AnnotationValueAcc::Categorical(counts)) => {
                        *counts.entry(label).or_insert(0) += 1;
                    }
                    Some(AnnotationValueAcc::Numeric(_)) => {}
                    None => {
                        let mut counts = AHashMap::new();
                        counts.insert(label, 1);
                        self.entries
                            .insert(key_id, AnnotationValueAcc::Categorical(counts));
                    }
                },
                ParsedTag::Skip => {}
            }
        }
    }

    pub(crate) fn merge(
        &mut self,
        other: Self,
        other_keys: &AnnotationKeyTable,
        global_keys: &mut AnnotationKeyTable,
    ) {
        for (other_kid, other_val) in other.entries {
            let Some(name) = other_keys.get(other_kid) else {
                continue;
            };
            let Some(global_kid) = global_keys.intern(name) else {
                continue;
            };
            match (self.entries.get_mut(&global_kid), other_val) {
                (Some(AnnotationValueAcc::Numeric(dst)), AnnotationValueAcc::Numeric(src)) => {
                    dst.merge(src);
                }
                (
                    Some(AnnotationValueAcc::Categorical(dst)),
                    AnnotationValueAcc::Categorical(src),
                ) => {
                    for (cat, cnt) in src {
                        *dst.entry(cat).or_insert(0) += cnt;
                    }
                }
                (None, val) => {
                    self.entries.insert(global_kid, val);
                }
                _ => {}
            }
        }
    }

    pub(crate) fn iter_mut(&mut self) -> impl Iterator<Item = (usize, &mut AnnotationValueAcc)> {
        self.entries.iter_mut().map(|(&k, v)| (k, v))
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

enum ParsedTag {
    Numeric(f64),
    Categorical(String),
    Skip,
}

fn extract_numeric_or_cat(val: &AnnotationValue) -> ParsedTag {
    match val {
        AnnotationValue::Float(v) if v.is_finite() => ParsedTag::Numeric(*v),
        AnnotationValue::Float(_) => ParsedTag::Skip,
        AnnotationValue::Int(v) => ParsedTag::Numeric(*v as f64),
        AnnotationValue::Bool(b) => ParsedTag::Categorical(b.to_string()),
        AnnotationValue::Str(s) => {
            let trimmed = s.trim();
            if trimmed.is_empty() {
                ParsedTag::Skip
            } else if let Ok(v) = trimmed.parse::<f64>()
                && v.is_finite()
            {
                ParsedTag::Numeric(v)
            } else {
                ParsedTag::Categorical(trimmed.to_string())
            }
        }
        AnnotationValue::Array(_) => ParsedTag::Skip,
    }
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;

    use super::*;

    #[test]
    fn test_online_stats_merge_and_hpd() {
        let mut a = MetricAccumulator::new(true);
        let mut b = MetricAccumulator::new(true);
        for x in 1..=50 {
            a.push(x as f64);
        }
        for x in 51..=100 {
            b.push(x as f64);
        }
        a.merge(b);
        let summary = a.finalize(0.95);
        assert_eq!(summary.count, 100);
        assert_approx_eq!(summary.mean.unwrap(), 50.5, 1e-9);
        assert_approx_eq!(summary.median.unwrap(), 50.5, 1e-9);
        assert_approx_eq!(summary.min.unwrap(), 1.0, 1e-9);
        assert_approx_eq!(summary.max.unwrap(), 100.0, 1e-9);
        assert_eq!(summary.hpd_low, Some(1.0));
        assert_eq!(summary.hpd_high, Some(96.0));
    }

    #[test]
    fn test_hpd_skewed_distribution_picks_narrowest_window() {
        // Tight cluster near 1.0..1.4 plus one outlier at 100.0 -> 80% HPD excludes 100.0
        let mut samples = vec![1.0, 1.1, 1.2, 1.3, 1.4, 100.0];
        samples.sort_unstable_by(f64::total_cmp);
        let (lo, hi) = compute_hpd(&samples, 0.80).unwrap();
        assert_approx_eq!(lo, 1.0, 1e-9);
        assert_approx_eq!(hi, 1.4, 1e-9);
    }
}
