//! Charts that compare samples: series of points, a matrix of samples on
//! shared x values, and the spread of samples at each x.

/// Reports with more samples than this summarize them as a band or a
/// heatmap instead of one line per sample.
pub(super) const MAX_OVERLAID: usize = 10;
/// Most outlying samples drawn on a band.
const MAX_OUTLIERS: usize = 8;

/// One series of (x, y) points.
#[derive(Debug, Clone, Default, PartialEq)]
pub(super) struct Points {
    pub(super) x: Vec<f64>,
    pub(super) y: Vec<f64>,
}

impl Points {
    /// Points from a histogram where the index is the x value.
    pub(super) fn indexed(counts: &[f64]) -> Self {
        Self {
            x: (0..counts.len()).map(|i| i as f64).collect(),
            y: counts.to_vec(),
        }
    }

    pub(super) fn pairs(&self) -> Vec<(f64, f64)> {
        self.x.iter().copied().zip(self.y.iter().copied()).collect()
    }

    /// The points as a share of their total, in percent, so samples of
    /// different depths compare.
    pub(super) fn shares(&self) -> Self {
        let total: f64 = self.y.iter().sum();
        Self {
            x: self.x.clone(),
            y: self
                .y
                .iter()
                .map(|&y| if total > 0.0 { 100.0 * y / total } else { 0.0 })
                .collect(),
        }
    }
}

/// Values of many samples on a shared x grid: one row per sample, one
/// column per x value.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct SampleMatrix {
    /// Every x value of any sample, ascending.
    pub(super) x: Vec<f64>,
    pub(super) values: Vec<Vec<Option<f64>>>,
}

impl SampleMatrix {
    /// Aligns [samples] on the x values of all of them. With
    /// [missing_is_zero], x values a sample lacks count as zero, as in
    /// histograms; otherwise they stay empty, as past the end of shorter
    /// reads.
    pub(super) fn of(samples: &[Points], missing_is_zero: bool) -> Self {
        let mut x: Vec<f64> = samples.iter().flat_map(|s| s.x.iter().copied()).collect();
        x.sort_by(f64::total_cmp);
        x.dedup();
        let fill = missing_is_zero.then_some(0.0);
        let values = samples
            .iter()
            .map(|s| {
                let mut row = vec![fill; x.len()];
                for (&sx, &sy) in s.x.iter().zip(&s.y) {
                    if let Ok(col) = x.binary_search_by(|v| v.total_cmp(&sx)) {
                        row[col] = Some(sy);
                    }
                }
                row
            })
            .collect();
        Self { x, values }
    }

    /// Largest value of any cell, or zero.
    pub(super) fn max_value(&self) -> f64 {
        self.values
            .iter()
            .flatten()
            .flatten()
            .fold(0.0, |a: f64, &b| a.max(b))
    }
}

/// How samples spread at each x: percentiles across samples, and the
/// samples furthest from the median.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct SampleSpread {
    pub(super) x: Vec<f64>,
    pub(super) p10: Vec<f64>,
    pub(super) q1: Vec<f64>,
    pub(super) median: Vec<f64>,
    pub(super) q3: Vec<f64>,
    pub(super) p90: Vec<f64>,
    /// Rows of the matrix far from the median, most distant first.
    pub(super) outliers: Vec<usize>,
}

impl SampleSpread {
    /// Percentiles of each column of [matrix]. Outliers are samples whose
    /// mean distance from the median is above the Tukey fence of all
    /// samples, at most [`MAX_OUTLIERS`] of them.
    pub(super) fn of(matrix: &SampleMatrix) -> Self {
        let levels = [0.1, 0.25, 0.5, 0.75, 0.9];
        let mut x = Vec::new();
        let mut rows: [Vec<f64>; 5] = Default::default();
        // Median of each column, or none for empty columns.
        let mut median_at = vec![None; matrix.x.len()];
        for (col, &cx) in matrix.x.iter().enumerate() {
            let mut column: Vec<f64> = matrix.values.iter().filter_map(|r| r[col]).collect();
            if column.is_empty() {
                continue;
            }
            column.sort_by(f64::total_cmp);
            x.push(cx);
            for (row, &level) in rows.iter_mut().zip(&levels) {
                row.push(percentile(&column, level));
            }
            median_at[col] = Some(percentile(&column, 0.5));
        }
        let outliers = outliers(matrix, &median_at);
        let [p10, q1, median, q3, p90] = rows;
        Self {
            x,
            p10,
            q1,
            median,
            q3,
            p90,
            outliers,
        }
    }
}

fn outliers(matrix: &SampleMatrix, median_at: &[Option<f64>]) -> Vec<usize> {
    let distance: Vec<f64> = matrix
        .values
        .iter()
        .map(|row| {
            let gaps: Vec<f64> = row
                .iter()
                .zip(median_at)
                .filter_map(|(v, m)| Some((v.as_ref()? - m.as_ref()?).abs()))
                .collect();
            if gaps.is_empty() {
                0.0
            } else {
                gaps.iter().sum::<f64>() / gaps.len() as f64
            }
        })
        .collect();
    if distance.len() < 4 {
        return Vec::new();
    }
    let mut sorted = distance.clone();
    sorted.sort_by(f64::total_cmp);
    let q1 = percentile(&sorted, 0.25);
    let q3 = percentile(&sorted, 0.75);
    let fence = q3 + 1.5 * (q3 - q1);
    let mut outliers: Vec<usize> = (0..distance.len())
        .filter(|&i| distance[i] > fence)
        .collect();
    outliers.sort_by(|&a, &b| distance[b].total_cmp(&distance[a]));
    outliers.truncate(MAX_OUTLIERS);
    outliers
}

/// Linear interpolation between the closest ranks of [sorted].
fn percentile(sorted: &[f64], level: f64) -> f64 {
    let rank = (sorted.len() - 1) as f64 * level;
    let low = rank.floor() as usize;
    let high = rank.ceil() as usize;
    sorted[low] + (sorted[high] - sorted[low]) * (rank - low as f64)
}

#[cfg(test)]
mod test {
    use super::*;

    fn points(x: &[f64], y: &[f64]) -> Points {
        Points {
            x: x.to_vec(),
            y: y.to_vec(),
        }
    }

    #[test]
    fn test_shares() {
        let shares = points(&[0.0, 1.0], &[1.0, 3.0]).shares();
        assert_eq!(shares.y, vec![25.0, 75.0]);
        assert_eq!(points(&[0.0], &[0.0]).shares().y, vec![0.0]);
        assert_eq!(Points::indexed(&[4.0, 5.0]).x, vec![0.0, 1.0]);
    }

    #[test]
    fn test_matrix() {
        let samples = [
            points(&[1.0, 2.0], &[5.0, 6.0]),
            points(&[2.0, 3.0], &[7.0, 8.0]),
        ];
        let matrix = SampleMatrix::of(&samples, true);
        assert_eq!(matrix.x, vec![1.0, 2.0, 3.0]);
        assert_eq!(matrix.values[0], vec![Some(5.0), Some(6.0), Some(0.0)]);
        let matrix = SampleMatrix::of(&samples, false);
        assert_eq!(matrix.values[1], vec![None, Some(7.0), Some(8.0)]);
        assert_eq!(matrix.max_value(), 8.0);
    }

    #[test]
    fn test_spread() {
        let mut samples: Vec<Points> = (0..6)
            .map(|i| points(&[0.0, 1.0], &[10.0 + i as f64 * 0.1, 10.0]))
            .collect();
        samples.push(points(&[0.0, 1.0], &[90.0, 90.0]));
        let spread = SampleSpread::of(&SampleMatrix::of(&samples, false));
        assert_eq!(spread.x, vec![0.0, 1.0]);
        assert_eq!(spread.median[1], 10.0);
        assert!(spread.p90[0] > spread.q3[0]);
        assert_eq!(spread.outliers, vec![6]);
        assert_eq!(percentile(&[1.0, 2.0, 3.0, 4.0], 0.5), 2.5);
    }
}
