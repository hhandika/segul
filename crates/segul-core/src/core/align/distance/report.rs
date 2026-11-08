//! Pool pairwise distances across loci into a [`DistanceReport`].
use std::collections::HashMap;

use super::group::{GroupDistances, GroupStats, SampleGrouping};
use super::{LocusDistance, PairDistances};
use crate::helper::types::DataType;
use crate::stats::distance::{DistanceModel, DistanceSummary};
use crate::stats::distance_report::{
    DistanceHistogram, DistanceOverview, DistanceReport, GroupPairReport, GroupReport,
    GroupSummary, LocusDistanceReport, PairStats, SampleDistanceReport, SampleMatrix,
};

/// Width of the pairwise distance histogram bins.
const BIN_WIDTH: f64 = 0.001;
/// Distances of 10 or more go in the last bin.
const MAX_BINS: usize = 10_000;
/// The sample matrix is left out above this many samples,
/// since its size grows with the square of the sample count.
pub(super) const MAX_MATRIX_SAMPLES: usize = 1_000;

/// Running statistics of the distances of every locus.
#[derive(Debug)]
pub(super) struct ReportAccumulator {
    /// Sample names in the order first seen.
    names: Vec<String>,
    index: HashMap<String, usize>,
    /// Loci each sample appears in.
    loci: Vec<usize>,
    /// Distances of each sample to every other sample.
    samples: Vec<GroupStats>,
    /// `matrix[i][j]` for `j < i` holds the distance sum and locus count of
    /// samples `i` and `j`. `None` once there are more than `max_matrix` samples.
    matrix: Option<Vec<Vec<(f64, u32)>>>,
    max_matrix: usize,
    all: GroupStats,
    histogram: Vec<usize>,
}

impl ReportAccumulator {
    pub(super) fn new(max_matrix: usize) -> Self {
        Self {
            names: Vec::new(),
            index: HashMap::new(),
            loci: Vec::new(),
            samples: Vec::new(),
            matrix: Some(Vec::new()),
            max_matrix,
            all: GroupStats::default(),
            histogram: Vec::new(),
        }
    }

    /// Add the pairwise distances of a locus. `ids[i]` is sample `i` of `distances`.
    pub(super) fn add(&mut self, ids: &[String], distances: &PairDistances) {
        let positions: Vec<usize> = ids.iter().map(|id| self.position(id)).collect();
        positions.iter().for_each(|&p| self.loci[p] += 1);
        for i in 0..distances.ntax {
            for j in i + 1..distances.ntax {
                let Some(d) = distances.get(i, j) else {
                    continue;
                };
                let (a, b) = (positions[i], positions[j]);
                self.samples[a].add(d);
                self.samples[b].add(d);
                self.all.add(d);
                let bin = (((d / BIN_WIDTH) + 1e-9).floor() as usize).min(MAX_BINS - 1);
                if bin >= self.histogram.len() {
                    self.histogram.resize(bin + 1, 0);
                }
                self.histogram[bin] += 1;
                if let Some(matrix) = &mut self.matrix {
                    let cell = &mut matrix[a.max(b)][a.min(b)];
                    cell.0 += d;
                    cell.1 += 1;
                }
            }
        }
    }

    /// Index of a sample, adding it if new.
    fn position(&mut self, id: &str) -> usize {
        if let Some(&p) = self.index.get(id) {
            return p;
        }
        let p = self.names.len();
        self.names.push(id.to_string());
        self.index.insert(id.to_string(), p);
        self.loci.push(0);
        self.samples.push(GroupStats::default());
        if self.names.len() > self.max_matrix {
            self.matrix = None;
        }
        if let Some(matrix) = &mut self.matrix {
            matrix.push(vec![(0.0, 0); p]);
        }
        p
    }

    /// Whether the sample matrix was left out for too many samples.
    pub(super) fn has_no_matrix(&self) -> bool {
        self.matrix.is_none()
    }

    /// Mean distance of samples `i` and `j`, or `None` if they share no locus.
    fn mean(matrix: &[Vec<(f64, u32)>], i: usize, j: usize) -> Option<f64> {
        if i == j {
            return Some(0.0);
        }
        let (sum, n) = matrix[i.max(j)][i.min(j)];
        (n > 0).then(|| sum / n as f64)
    }

    /// Build the report. `groups` is the pooled group statistics with grouping.
    pub(super) fn into_report(
        self,
        datatype: &DataType,
        model: &DistanceModel,
        loci: &[&LocusDistance],
        grouping: Option<(&SampleGrouping, &GroupDistances)>,
    ) -> DistanceReport {
        let mut order: Vec<usize> = (0..self.names.len()).collect();
        order.sort_by(|&a, &b| alphanumeric_sort::compare_str(&self.names[a], &self.names[b]));
        let group_of: Vec<Option<String>> = self
            .names
            .iter()
            .map(|n| grouping.and_then(|(g, _)| g.group_of(n)))
            .collect();
        let samples = order
            .iter()
            .map(|&i| {
                let mut sample = SampleDistanceReport {
                    name: self.names[i].clone(),
                    group: group_of[i].clone(),
                    loci: self.loci[i],
                    stats: pair_stats(&self.samples[i]),
                    nearest: None,
                    nearest_distance: None,
                    max_within: None,
                    min_between: None,
                };
                if let Some(matrix) = &self.matrix {
                    let others = (0..self.names.len())
                        .filter(|&j| j != i)
                        .filter_map(|j| Some((j, Self::mean(matrix, i, j)?)));
                    for (j, d) in others {
                        if sample.nearest_distance.is_none_or(|n| d < n) {
                            sample.nearest = Some(self.names[j].clone());
                            sample.nearest_distance = Some(d);
                        }
                        if let (Some(a), Some(b)) = (&group_of[i], &group_of[j]) {
                            if a == b {
                                sample.max_within = Some(sample.max_within.map_or(d, |m| m.max(d)));
                            } else {
                                sample.min_between =
                                    Some(sample.min_between.map_or(d, |m| m.min(d)));
                            }
                        }
                    }
                }
                sample
            })
            .collect();
        let matrix = self.matrix.as_ref().map(|matrix| SampleMatrix {
            names: order.iter().map(|&i| self.names[i].clone()).collect(),
            mean: order
                .iter()
                .map(|&i| order.iter().map(|&j| Self::mean(matrix, i, j)).collect())
                .collect(),
        });
        let groups = grouping.map(|(_, pooled)| {
            let names = pooled.names();
            GroupReport {
                groups: names
                    .iter()
                    .map(|name| GroupSummary {
                        name: name.to_string(),
                        samples: group_of
                            .iter()
                            .filter(|g| g.as_deref() == Some(*name))
                            .count(),
                        within: pair_stats(&pooled.get(name, name)),
                    })
                    .collect(),
                pairs: names
                    .iter()
                    .enumerate()
                    .flat_map(|(i, a)| {
                        names[i..].iter().map(|b| GroupPairReport {
                            group_1: a.to_string(),
                            group_2: b.to_string(),
                            stats: pair_stats(&pooled.get(a, b)),
                        })
                    })
                    .collect(),
            }
        });
        let pairs: usize = loci.iter().map(|l| l.pairs).sum();
        let compared: f64 = loci
            .iter()
            .filter_map(|l| Some(l.mean_compared_sites? * l.pairs as f64))
            .sum();
        let locus_means: Vec<f64> = loci.iter().filter_map(|l| l.mean).collect();
        let summary = DistanceOverview {
            loci: loci.len(),
            samples: self.names.len(),
            pairs,
            undefined_pairs: loci.iter().map(|l| l.undefined_pairs).sum(),
            min: self.all.min,
            max: self.all.max,
            mean: self.all.mean(),
            sd: self.all.sd(),
            median_locus_mean: DistanceSummary::new(&locus_means).median,
            mean_compared_sites: (pairs > 0).then(|| compared / pairs as f64),
            groups: grouping.map(|(_, pooled)| pooled.names().len()),
            unassigned: if grouping.is_some() {
                group_of.iter().filter(|g| g.is_none()).count()
            } else {
                0
            },
        };
        DistanceReport {
            datatype: if *datatype == DataType::Aa {
                "aa"
            } else {
                "dna"
            }
            .to_string(),
            model: model.to_string(),
            summary,
            histogram: DistanceHistogram {
                bin_width: BIN_WIDTH,
                counts: self.histogram,
            },
            loci: loci.iter().map(|l| locus_report(l)).collect(),
            samples,
            matrix,
            groups,
        }
    }
}

fn pair_stats(stats: &GroupStats) -> PairStats {
    PairStats {
        pairs: stats.n_pairs,
        mean: stats.mean(),
        sd: stats.sd(),
        min: stats.min,
        max: stats.max,
    }
}

fn locus_report(locus: &LocusDistance) -> LocusDistanceReport {
    LocusDistanceReport {
        path: locus.path.display().to_string(),
        name: locus.locus.clone(),
        ntax: locus.ntax,
        nchar: locus.nchar,
        pairs: locus.pairs,
        undefined_pairs: locus.undefined_pairs,
        min: locus.min,
        max: locus.max,
        mean: locus.mean,
        median: locus.median,
        sd: locus.sd,
        mean_compared_sites: locus.mean_compared_sites,
    }
}
