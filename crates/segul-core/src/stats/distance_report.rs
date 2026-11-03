//! Pairwise genetic distances as one JSON report.
//!
//! The report is written to `distance-report.json` next to the distance
//! CSV files, so applications can draw their own charts without parsing
//! the per-locus matrices.
use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

/// Every statistic of a pairwise distance run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DistanceReport {
    /// `dna` or `aa`.
    pub datatype: String,
    /// Distance model, such as `K2P`.
    pub model: String,
    pub summary: DistanceOverview,
    /// All defined pairwise distances across loci.
    pub histogram: DistanceHistogram,
    /// One entry per alignment, sorted by path.
    pub loci: Vec<LocusDistanceReport>,
    /// One entry per sample, sorted by name.
    pub samples: Vec<SampleDistanceReport>,
    /// Mean distance between samples across loci. `None` when there are
    /// too many samples.
    pub matrix: Option<SampleMatrix>,
    /// Distances within and between groups. `None` without grouping.
    pub groups: Option<GroupReport>,
}

/// Totals across all loci.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DistanceOverview {
    pub loci: usize,
    pub samples: usize,
    /// Sample pairs summed over loci.
    pub pairs: usize,
    /// Pairs without a distance: no compared sites or saturated.
    pub undefined_pairs: usize,
    /// Statistics of every defined pairwise distance across loci.
    pub min: Option<f64>,
    pub max: Option<f64>,
    pub mean: Option<f64>,
    pub sd: Option<f64>,
    /// Median of the per-locus mean distances.
    pub median_locus_mean: Option<f64>,
    /// Mean compared sites per pair, across loci.
    pub mean_compared_sites: Option<f64>,
    /// Number of groups. `None` without grouping.
    pub groups: Option<usize>,
    /// Samples that match no group.
    pub unassigned: usize,
}

/// Counts of pairwise distances in bins of equal width, starting at zero.
/// Distances beyond the last bin are counted in it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DistanceHistogram {
    pub bin_width: f64,
    pub counts: Vec<usize>,
}

/// Pairwise distances of one alignment.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LocusDistanceReport {
    pub path: String,
    pub name: String,
    pub ntax: usize,
    pub nchar: usize,
    pub pairs: usize,
    pub undefined_pairs: usize,
    pub min: Option<f64>,
    pub max: Option<f64>,
    pub mean: Option<f64>,
    pub median: Option<f64>,
    pub sd: Option<f64>,
    pub mean_compared_sites: Option<f64>,
}

/// Distances of one sample to every other sample, across loci.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SampleDistanceReport {
    pub name: String,
    /// `None` without grouping or if the sample matches no group.
    pub group: Option<String>,
    /// Loci the sample appears in.
    pub loci: usize,
    pub stats: PairStats,
    /// Sample with the lowest mean distance in the sample matrix.
    pub nearest: Option<String>,
    pub nearest_distance: Option<f64>,
    /// Highest mean distance to a sample of the same group, and lowest to
    /// a sample of another group, from the sample matrix.
    pub max_within: Option<f64>,
    pub min_between: Option<f64>,
}

/// Square matrix of mean distances across the loci two samples share.
/// `mean[i][j]` is `None` if samples `i` and `j` share no locus with a
/// defined distance.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SampleMatrix {
    pub names: Vec<String>,
    pub mean: Vec<Vec<Option<f64>>>,
}

/// Distances within and between groups, pooled across loci.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GroupReport {
    /// One entry per group, sorted by name.
    pub groups: Vec<GroupSummary>,
    /// Every pair of groups, including each group with itself.
    pub pairs: Vec<GroupPairReport>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GroupSummary {
    pub name: String,
    pub samples: usize,
    pub within: PairStats,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GroupPairReport {
    pub group_1: String,
    pub group_2: String,
    pub stats: PairStats,
}

/// Statistics of a set of pairwise distances.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct PairStats {
    /// Pairs with a defined distance.
    pub pairs: usize,
    pub mean: Option<f64>,
    pub sd: Option<f64>,
    pub min: Option<f64>,
    pub max: Option<f64>,
}

impl DistanceReport {
    /// Read a report written by `segul align distance`.
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

    /// Write the report as JSON to `path`.
    pub fn write(&self, path: &Path) -> Result<()> {
        fs::write(path, self.to_json()?)
            .with_context(|| format!("Failed writing {}", path.display()))
    }
}
