//! Contig summary statistics as one JSON report.
//!
//! The report is written to `contig-report.json` next to the contig
//! summary CSV, so applications can draw their own charts without parsing
//! the CSV.
use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use crate::stats::contigs::{ContigSummary, LENGTH_BIN_EDGES, proportion};

/// Every statistic of a contig summary.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ContigReport {
    pub summary: ContigOverview,
    /// Contig length bins, shortest first. Each sample's
    /// `length_distribution` aligns with these.
    pub length_bins: Vec<LengthBin>,
    /// One entry per contig file, in input order.
    pub samples: Vec<SampleContigReport>,
}

/// Totals across all contig files.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ContigOverview {
    pub samples: usize,
    pub contigs: usize,
    /// A, C, G, T, and N.
    pub bases: usize,
    /// A, C, G, and T.
    pub nucleotides: usize,
    pub n_count: usize,
    /// Characters other than A, C, G, T, and N.
    pub unknown: usize,
    pub gc: f64,
    pub at: f64,
    /// N over all bases.
    pub prop_n: f64,
    /// Fewest, most, and mean contigs in one sample.
    pub min_contigs: usize,
    pub max_contigs: usize,
    pub mean_contigs: f64,
    /// Lowest, highest, and median N50 across samples.
    pub min_n50: usize,
    pub max_n50: usize,
    pub median_n50: f64,
    /// Length of the longest contig.
    pub longest: usize,
}

/// Contigs from `min` bp up to, not including, `max` bp.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LengthBin {
    pub min: usize,
    /// `None` for the last bin, which has no upper bound.
    pub max: Option<usize>,
}

/// Contigs longer than `length` bp.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LongContigCount {
    pub length: usize,
    pub contigs: usize,
}

/// Statistics of one contig file.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SampleContigReport {
    pub path: String,
    pub name: String,
    pub contigs: usize,
    pub bases: usize,
    pub nucleotides: usize,
    pub g_count: usize,
    pub c_count: usize,
    pub a_count: usize,
    pub t_count: usize,
    pub n_count: usize,
    pub unknown: usize,
    pub gc: f64,
    pub at: f64,
    pub prop_n: f64,
    /// Contig length statistics.
    pub min: usize,
    pub max: usize,
    pub mean: f64,
    pub median: f64,
    pub sd: f64,
    pub n50: usize,
    pub n75: usize,
    pub n90: usize,
    /// Contigs longer than 750, 1000, and 1500 bp.
    pub long_contigs: Vec<LongContigCount>,
    /// Contigs per length bin, aligned with [`ContigReport::length_bins`].
    pub length_distribution: Vec<usize>,
}

impl LengthBin {
    /// Short label, such as `500-1k` or `100k+`.
    pub fn label(&self) -> String {
        match self.max {
            Some(max) => format!("{}-{}", short_length(self.min), short_length(max)),
            None => format!("{}+", short_length(self.min)),
        }
    }
}

impl ContigReport {
    /// Build the report from the summaries of each contig file.
    pub fn new(summaries: &[ContigSummary]) -> Self {
        let samples: Vec<SampleContigReport> =
            summaries.iter().map(SampleContigReport::new).collect();
        let length_bins = LENGTH_BIN_EDGES
            .iter()
            .enumerate()
            .map(|(i, &min)| LengthBin {
                min,
                max: LENGTH_BIN_EDGES.get(i + 1).copied(),
            })
            .collect();
        Self {
            summary: ContigOverview::new(&samples),
            length_bins,
            samples,
        }
    }

    /// Read a report written by `segul contig summary`.
    pub fn from_path(path: &Path) -> Result<Self> {
        let text = fs::read_to_string(path)
            .with_context(|| format!("Failed reading {}", path.display()))?;
        Self::from_json(&text).with_context(|| format!("Failed parsing {}", path.display()))
    }

    pub fn from_json(json: &str) -> Result<Self> {
        crate::report::compat::parse_report_payload(json)
    }

    pub fn to_json(&self) -> Result<String> {
        Ok(serde_json::to_string(self)?)
    }

    /// Write the report as JSON to `path`.
    pub fn write(&self, path: &Path) -> Result<()> {
        fs::write(path, self.to_json()?)
            .with_context(|| format!("Failed writing {}", path.display()))
    }
}

impl ContigOverview {
    fn new(samples: &[SampleContigReport]) -> Self {
        let total = |value: fn(&SampleContigReport) -> usize| samples.iter().map(value).sum();
        let bases: usize = total(|s| s.bases);
        let n_count: usize = total(|s| s.n_count);
        let gc: usize = total(|s| s.g_count + s.c_count);
        let at: usize = total(|s| s.a_count + s.t_count);
        let contigs: usize = total(|s| s.contigs);
        let mut n50: Vec<usize> = samples.iter().map(|s| s.n50).collect();
        n50.sort_unstable();
        Self {
            samples: samples.len(),
            contigs,
            bases,
            nucleotides: total(|s| s.nucleotides),
            n_count,
            unknown: total(|s| s.unknown),
            gc: proportion(gc, bases),
            at: proportion(at, bases),
            prop_n: proportion(n_count, bases),
            min_contigs: samples.iter().map(|s| s.contigs).min().unwrap_or(0),
            max_contigs: samples.iter().map(|s| s.contigs).max().unwrap_or(0),
            mean_contigs: proportion(contigs, samples.len()),
            min_n50: n50.first().copied().unwrap_or(0),
            max_n50: n50.last().copied().unwrap_or(0),
            median_n50: median(&n50),
            longest: samples.iter().map(|s| s.max).max().unwrap_or(0),
        }
    }
}

impl SampleContigReport {
    fn new(summary: &ContigSummary) -> Self {
        Self {
            path: summary.file_path.clone(),
            name: summary.contig_name.clone(),
            contigs: summary.contig_count,
            bases: summary.base_count,
            nucleotides: summary.nucleotide,
            g_count: summary.g_count,
            c_count: summary.c_count,
            a_count: summary.a_count,
            t_count: summary.t_count,
            n_count: summary.n_count,
            unknown: summary.unknown,
            gc: summary.gc_content,
            at: summary.at_content,
            prop_n: proportion(summary.n_count, summary.base_count),
            min: summary.stats.min,
            max: summary.stats.max,
            mean: summary.stats.mean,
            median: summary.stats.median,
            sd: summary.stats.stdev,
            n50: summary.nstats.n50,
            n75: summary.nstats.n75,
            n90: summary.nstats.n90,
            long_contigs: vec![
                LongContigCount {
                    length: 750,
                    contigs: summary.contig750,
                },
                LongContigCount {
                    length: 1000,
                    contigs: summary.contig1000,
                },
                LongContigCount {
                    length: 1500,
                    contigs: summary.contig1500,
                },
            ],
            length_distribution: summary.length_distribution.clone(),
        }
    }
}

/// Median of sorted values, 0 when empty.
fn median(sorted: &[usize]) -> f64 {
    let n = sorted.len();
    match n {
        0 => 0.0,
        _ if n.is_multiple_of(2) => (sorted[n / 2 - 1] + sorted[n / 2]) as f64 / 2.0,
        _ => sorted[n / 2] as f64,
    }
}

/// Length in bp, with whole thousands as `k`.
fn short_length(bp: usize) -> String {
    if bp >= 1_000 && bp.is_multiple_of(1_000) {
        format!("{}k", bp / 1_000)
    } else {
        bp.to_string()
    }
}

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use super::*;
    use crate::helper::types::ContigFmt;

    fn summaries() -> Vec<ContigSummary> {
        ["contigs1.fa", "contigs2.fa"]
            .iter()
            .map(|f| {
                let mut summary = ContigSummary::new();
                summary.summarize(
                    &PathBuf::from("tests/files/contigs").join(f),
                    &ContigFmt::Fasta,
                );
                summary
            })
            .collect()
    }

    #[test]
    fn test_contig_report() {
        let report = ContigReport::new(&summaries());
        let s = &report.summary;
        assert_eq!(s.samples, 2);
        assert_eq!(s.contigs, 4);
        assert_eq!(s.bases, 485 + 552);
        assert_eq!((s.min_n50, s.max_n50, s.median_n50), (258, 316, 287.0));
        assert_eq!(s.longest, 316);
        assert_eq!(report.samples[0].name, "contigs1");
        assert_eq!(report.samples[0].long_contigs[1].length, 1000);
        assert_eq!(report.length_bins.len(), LENGTH_BIN_EDGES.len());
        assert_eq!(
            report.samples[0].length_distribution.len(),
            report.length_bins.len()
        );
        let binned: usize = report
            .samples
            .iter()
            .flat_map(|s| &s.length_distribution)
            .sum();
        assert_eq!(binned, s.contigs);
    }

    #[test]
    fn test_round_trip() {
        let report = ContigReport::new(&summaries());
        let json = report.to_json().unwrap();
        let parsed = ContigReport::from_json(&json).unwrap();
        assert_eq!(parsed.samples[0].name, report.samples[0].name);
        assert_eq!(parsed.summary.contigs, report.summary.contigs);
    }

    #[test]
    fn test_empty_report() {
        let report = ContigReport::new(&[ContigSummary::new()]);
        let json = report.to_json().unwrap();
        assert_eq!(ContigReport::from_json(&json).unwrap(), report);
        assert_eq!(report.summary.gc, 0.0);
    }

    #[test]
    fn test_bin_labels() {
        let report = ContigReport::new(&[]);
        let labels: Vec<String> = report.length_bins.iter().map(|b| b.label()).collect();
        assert_eq!(labels[0], "0-200");
        assert_eq!(labels[2], "500-1k");
        assert_eq!(labels.last().unwrap(), "100k+");
    }
}
