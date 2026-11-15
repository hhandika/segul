//! Samples of a read report with their checks, the metrics shown for
//! them, and the order of the sample list.

use std::cmp::Ordering;

use segul_core::stats::fastq::qc::{
    ADAPTERS_RULE, DUPLICATION_RULE, MEAN_QUALITY_RULE, Q30_RULE, QcCheck, QcRule, QcStatus,
    overall_status, qc_checks,
};
use segul_core::stats::fastq::report::{ReadMetrics, ReadReport, SampleReadReport};

use crate::cli::view::thousands;

/// A sample with its quality checks.
#[derive(Debug, Clone)]
pub(super) struct ReadSample {
    pub(super) report: SampleReadReport,
    pub(super) checks: Vec<QcCheck>,
    pub(super) status: QcStatus,
}

impl ReadSample {
    pub(super) fn new(report: SampleReadReport) -> Self {
        let checks = qc_checks(&report);
        let status = overall_status(&checks);
        Self {
            report,
            checks,
            status,
        }
    }

    pub(super) fn name(&self) -> &str {
        &self.report.sample
    }

    pub(super) fn metrics(&self) -> &ReadMetrics {
        &self.report.metrics
    }
}

/// Samples with at least one read, which are the only ones with charts.
pub(super) fn samples_with_reads(report: ReadReport) -> Vec<ReadSample> {
    report
        .samples
        .into_iter()
        .filter(|s| s.metrics.read_count > 0)
        .map(ReadSample::new)
        .collect()
}

/// A per-sample statistic shown on the summary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ReadMetric {
    Reads,
    Bases,
    MeanLength,
    N50,
    Q30,
    MeanQuality,
    Gc,
    Duplication,
    Adapters,
}

impl ReadMetric {
    pub(super) const ALL: [Self; 9] = [
        Self::Reads,
        Self::Bases,
        Self::MeanLength,
        Self::N50,
        Self::Q30,
        Self::MeanQuality,
        Self::Gc,
        Self::Duplication,
        Self::Adapters,
    ];

    /// Columns of the summary table, in the order of the GUI.
    pub(super) const TABLE: [Self; 9] = [
        Self::Reads,
        Self::Bases,
        Self::MeanLength,
        Self::N50,
        Self::Gc,
        Self::Q30,
        Self::MeanQuality,
        Self::Duplication,
        Self::Adapters,
    ];

    /// Name above the bars.
    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Reads => "Reads",
            Self::Bases => "Bases",
            Self::MeanLength => "Mean length (bp)",
            Self::N50 => "N50 (bp)",
            Self::Q30 => "Q30 bases (%)",
            Self::MeanQuality => "Mean quality",
            Self::Gc => "GC (%)",
            Self::Duplication => "Duplication (%)",
            Self::Adapters => "Adapters (%)",
        }
    }

    /// Heading of the table column.
    pub(super) fn column_label(self) -> &'static str {
        match self {
            Self::MeanLength => "Mean length",
            Self::N50 => "N50",
            Self::Q30 => "Q30 (%)",
            Self::MeanQuality => "Mean Q",
            _ => self.label(),
        }
    }

    /// The metric of [metrics], or none when the run did not compute it.
    pub(super) fn value_of(self, metrics: &ReadMetrics) -> Option<f64> {
        match self {
            Self::Reads => Some(metrics.read_count as f64),
            Self::Bases => Some(metrics.base_count as f64),
            Self::MeanLength => Some(metrics.mean_read_length),
            Self::N50 => Some(metrics.n50_read_length as f64),
            Self::Q30 => Some(metrics.q30_pct),
            Self::MeanQuality => Some(metrics.mean_qscore),
            Self::Gc => Some(100.0 * metrics.gc_content),
            Self::Duplication => metrics.duplicate_pct_est,
            Self::Adapters => metrics.adapter_reads_pct,
        }
    }

    /// Thresholds of the check on this metric, if any.
    pub(super) fn rule(self) -> Option<QcRule> {
        match self {
            Self::Q30 => Some(Q30_RULE),
            Self::MeanQuality => Some(MEAN_QUALITY_RULE),
            Self::Duplication => Some(DUPLICATION_RULE),
            Self::Adapters => Some(ADAPTERS_RULE),
            _ => None,
        }
    }

    /// [value] beside a bar, with counts shortened.
    pub(super) fn short(self, value: f64) -> String {
        match self {
            Self::Reads | Self::Bases => compact(value),
            Self::N50 => format!("{}", value.round()),
            _ => format!("{value:.1}"),
        }
    }

    /// [value] in a table cell.
    pub(super) fn exact(self, value: f64) -> String {
        match self {
            Self::Reads | Self::Bases | Self::N50 => thousands(value.round() as i128),
            _ => format!("{value:.1}"),
        }
    }

    /// Metrics computed for every sample.
    pub(super) fn available_for(samples: &[ReadSample]) -> Vec<Self> {
        Self::ALL
            .into_iter()
            .filter(|m| samples.iter().all(|s| m.value_of(s.metrics()).is_some()))
            .collect()
    }
}

/// A count with a K, M, or G suffix, such as 4.1M.
pub(super) fn compact(value: f64) -> String {
    for (size, suffix) in [(1e9, "G"), (1e6, "M"), (1e3, "K")] {
        if value.abs() >= size {
            return format!("{:.1}{suffix}", value / size);
        }
    }
    format!("{}", value.round())
}

/// Orders of the sample list.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SampleSort {
    Name,
    Status,
    Reads,
    Q30,
}

impl SampleSort {
    pub(super) fn next(self) -> Self {
        match self {
            Self::Name => Self::Status,
            Self::Status => Self::Reads,
            Self::Reads => Self::Q30,
            Self::Q30 => Self::Name,
        }
    }

    /// Name of the order and its direction.
    pub(super) fn describe(self, ascending: bool) -> String {
        let (label, up, down) = match self {
            Self::Name => ("Name", "A to Z", "Z to A"),
            Self::Status => ("Status", "failing last", "failing first"),
            Self::Reads => ("Reads", "fewest first", "most first"),
            Self::Q30 => ("Q30 bases", "lowest first", "highest first"),
        };
        format!("{label}, {}", if ascending { up } else { down })
    }
}

/// Positions in [samples] of those whose name contains [query], ignoring
/// case, ordered by [sort]. Names sort naturally, so `s2` comes before
/// `s10`, and break ties in the other orders.
pub(super) fn sample_view(
    samples: &[ReadSample],
    query: &str,
    sort: SampleSort,
    ascending: bool,
) -> Vec<usize> {
    let needle = query.trim().to_lowercase();
    let mut view: Vec<usize> = (0..samples.len())
        .filter(|&i| samples[i].name().to_lowercase().contains(&needle))
        .collect();
    let by_value = |a: &ReadSample, b: &ReadSample| match sort {
        SampleSort::Name => Ordering::Equal,
        SampleSort::Status => a.status.cmp(&b.status),
        SampleSort::Reads => a.metrics().read_count.cmp(&b.metrics().read_count),
        SampleSort::Q30 => a.metrics().q30_pct.total_cmp(&b.metrics().q30_pct),
    };
    view.sort_by(|&a, &b| {
        let (a, b) = (&samples[a], &samples[b]);
        let ordered =
            by_value(a, b).then_with(|| alphanumeric_sort::compare_str(a.name(), b.name()));
        if ascending {
            ordered
        } else {
            ordered.reverse()
        }
    });
    view
}

/// Sample [delta] steps from [index] through [view], wrapping at both
/// ends. Starts at either end when [index] is filtered out.
pub(super) fn step_through(view: &[usize], index: usize, delta: isize) -> usize {
    if view.is_empty() {
        return index;
    }
    match view.iter().position(|&i| i == index) {
        Some(position) => {
            let len = view.len() as isize;
            view[(position as isize + delta).rem_euclid(len) as usize]
        }
        None if delta > 0 => view[0],
        None => view[view.len() - 1],
    }
}

#[cfg(test)]
pub(super) mod test {
    use super::*;

    /// Samples named [names] with [reads] reads and [q30] percent Q30
    /// bases each.
    pub(in super::super) fn samples(names: &[&str], reads: &[u64], q30: &[f64]) -> Vec<ReadSample> {
        names
            .iter()
            .zip(reads.iter().zip(q30))
            .map(|(name, (&reads, &q30))| {
                let mut report = SampleReadReport {
                    sample: name.to_string(),
                    ..Default::default()
                };
                report.metrics.read_count = reads;
                report.metrics.q30_pct = q30;
                report.metrics.mean_qscore = 35.0;
                ReadSample::new(report)
            })
            .collect()
    }

    #[test]
    fn test_sample_view() {
        let samples = samples(&["s10", "s2", "a1"], &[5, 9, 1], &[90.0, 60.0, 75.0]);
        assert_eq!(sample_view(&samples, "", SampleSort::Name, true), [2, 1, 0]);
        assert_eq!(sample_view(&samples, "S", SampleSort::Name, false), [0, 1]);
        assert_eq!(
            sample_view(&samples, "", SampleSort::Reads, false),
            [1, 0, 2]
        );
        // Pass, fail, warn by Q30.
        assert_eq!(
            sample_view(&samples, "", SampleSort::Status, true),
            [0, 2, 1]
        );
        assert_eq!(sample_view(&samples, "", SampleSort::Q30, true), [1, 2, 0]);
    }

    #[test]
    fn test_step_through() {
        let view = [2, 0, 1];
        assert_eq!(step_through(&view, 0, 1), 1);
        assert_eq!(step_through(&view, 1, 1), 2);
        assert_eq!(step_through(&view, 2, -1), 1);
        assert_eq!(step_through(&[0, 1], 5, -1), 1);
        assert_eq!(step_through(&[], 5, 1), 5);
    }

    #[test]
    fn test_metrics() {
        let samples = samples(&["a"], &[4_100_000], &[90.0]);
        assert_eq!(ReadMetric::Reads.short(4_100_000.0), "4.1M");
        assert_eq!(ReadMetric::Reads.exact(4_100_000.0), "4,100,000");
        assert_eq!(compact(950.0), "950");
        let available = ReadMetric::available_for(&samples);
        assert!(available.contains(&ReadMetric::Q30));
        assert!(!available.contains(&ReadMetric::Duplication));
        assert_eq!(ReadMetric::Q30.rule(), Some(Q30_RULE));
    }
}
