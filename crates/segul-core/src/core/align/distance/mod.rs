//! Compute pairwise genetic distances for each locus alignment.
//!
//! Each alignment gets its own square distance matrix,
//! written to `output_dir/per_locus/<locus>.csv`.
//! A summary of each locus's pairwise distances
//! is written to `output_dir/<prefix>_summary.csv`.
//!
//! With a [`SampleGrouping`], the mean distance within and between groups,
//! such as populations or species, is written for each locus to
//! `output_dir/per_locus_group/<locus>.csv` and pooled across loci to
//! `output_dir/<prefix>_group_mean.csv`, with pair counts in
//! `output_dir/<prefix>_group_pairs.csv`. The diagonal of the mean matrices
//! follows [`GroupDiagonal`] and defaults to the within-group mean.
//!
//! With long format enabled, all pairwise distances across loci are written to
//! `output_dir/<prefix>_pairwise_long.csv`, and with grouping, all group
//! statistics to `output_dir/<prefix>_group_long.csv`.
pub mod group;

use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc;

use anyhow::{Context, Result};
use colored::Colorize;
use rayon::prelude::*;
use serde::Serialize;

use crate::helper::sequence::SeqParser;
use crate::helper::types::{DataType, InputFmt, SeqMatrix};
use crate::helper::utils;
use crate::stats::distance::{self, BaseFreqs, DistanceModel, DistanceSummary, PairCounts};

use group::{GroupAssignment, GroupDistances, SampleGrouping};

const PER_LOCUS_DIR: &str = "per_locus";
const PER_LOCUS_GROUP_DIR: &str = "per_locus_group";
const LONG_TMP_DIR: &str = ".pairwise_long_tmp";
const DEFAULT_PREFIX: &str = "distance";
/// Locus name of the pooled rows in the long group table.
const ALL_LOCI: &str = "all";

/// Which part of the distance matrix to write.
/// The diagonal is always written, see [`MatrixDiagonal`].
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum MatrixTriangle {
    /// Write both triangles.
    #[default]
    Full,
    /// Write the lower triangle only. Cells above the diagonal are empty.
    Lower,
    /// Write the upper triangle only. Cells below the diagonal are empty.
    Upper,
}

impl MatrixTriangle {
    /// Check whether the cell at row `i` and column `j` is written.
    fn includes(&self, i: usize, j: usize) -> bool {
        match self {
            Self::Full => true,
            Self::Lower => j <= i,
            Self::Upper => j >= i,
        }
    }
}

impl std::str::FromStr for MatrixTriangle {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "full" => Ok(Self::Full),
            "lower" => Ok(Self::Lower),
            "upper" => Ok(Self::Upper),
            _ => Err(format!("{} is not a valid matrix triangle", s)),
        }
    }
}

/// What to write on the diagonal of the distance matrix.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum MatrixDiagonal {
    /// Write zero, the distance of a sample to itself.
    #[default]
    Zero,
    /// Leave the diagonal empty.
    Blank,
    /// Write the highest distance of the sample to any other sample (the row maximum).
    /// Written as `NA` if the sample has no defined distances.
    RowMax,
}

impl std::str::FromStr for MatrixDiagonal {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "zero" => Ok(Self::Zero),
            "blank" => Ok(Self::Blank),
            "max" => Ok(Self::RowMax),
            _ => Err(format!("{} is not a valid matrix diagonal", s)),
        }
    }
}

/// What to write on the diagonal of the group mean matrices.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum GroupDiagonal {
    /// Write the mean distance within each group.
    #[default]
    Within,
    /// Write zero.
    Zero,
    /// Leave the diagonal empty.
    Blank,
    /// Write the highest mean distance of the group to any other group (the row maximum).
    /// Written as `NA` if the group has no defined mean to another group.
    RowMax,
}

impl std::str::FromStr for GroupDiagonal {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "within" => Ok(Self::Within),
            "zero" => Ok(Self::Zero),
            "blank" => Ok(Self::Blank),
            "max" => Ok(Self::RowMax),
            _ => Err(format!("{} is not a valid group matrix diagonal", s)),
        }
    }
}

/// Compute pairwise distances for alignments.
pub struct AlignmentDistance<'a> {
    /// The input format of the alignment files.
    input_fmt: &'a InputFmt,
    /// The data type of the alignment files. Use DNA or amino acid.
    datatype: &'a DataType,
    /// The distance model.
    model: &'a DistanceModel,
    /// The output directory.
    output_dir: &'a Path,
    /// Which part of the distance matrix to write. Defaults to the full matrix.
    triangle: MatrixTriangle,
    /// What to write on the diagonal. Defaults to zero.
    diagonal: MatrixDiagonal,
    /// How to group samples for group mean distances. Defaults to no grouping.
    grouping: Option<SampleGrouping>,
    /// What to write on the diagonal of the group mean matrices.
    /// Defaults to the within-group mean.
    group_diagonal: GroupDiagonal,
    /// Whether to write long-format tables. Defaults to false.
    long_format: bool,
}

impl<'a> AlignmentDistance<'a> {
    pub fn new(
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        model: &'a DistanceModel,
        output_dir: &'a Path,
    ) -> Self {
        Self {
            input_fmt,
            datatype,
            model,
            output_dir,
            triangle: MatrixTriangle::Full,
            diagonal: MatrixDiagonal::Zero,
            grouping: None,
            group_diagonal: GroupDiagonal::Within,
            long_format: false,
        }
    }

    /// Write only the lower or upper triangle of each distance matrix,
    /// leaving the other triangle empty.
    pub fn with_triangle(mut self, triangle: MatrixTriangle) -> Self {
        self.triangle = triangle;
        self
    }

    /// Set what to write on the diagonal of each distance matrix.
    pub fn with_diagonal(mut self, diagonal: MatrixDiagonal) -> Self {
        self.diagonal = diagonal;
        self
    }

    /// Compute the mean distance within and between groups of samples.
    pub fn with_grouping(mut self, grouping: SampleGrouping) -> Self {
        self.grouping = Some(grouping);
        self
    }

    /// Set what to write on the diagonal of the group mean matrices.
    pub fn with_group_diagonal(mut self, diagonal: GroupDiagonal) -> Self {
        self.group_diagonal = diagonal;
        self
    }

    /// Write all pairwise distances, and group statistics when grouping is set,
    /// as long-format tables with one row per pair across all loci.
    pub fn with_long_format(mut self, long_format: bool) -> Self {
        self.long_format = long_format;
        self
    }

    /// Compute pairwise distances for each alignment.
    /// Panics if the model is not valid for the data type.
    /// Example:
    /// ```rust
    /// use std::path::PathBuf;
    /// use segul_core::core::align::distance::AlignmentDistance;
    /// use segul_core::helper::types::{DataType, InputFmt};
    /// use segul_core::stats::distance::DistanceModel;
    /// use tempdir::TempDir;
    ///
    /// let files = vec![PathBuf::from("tests/files/alignments/gene_2.nex")];
    /// let output = TempDir::new("temp").unwrap();
    /// let handle = AlignmentDistance::new(
    ///     &InputFmt::Nexus,
    ///     &DataType::Dna,
    ///     &DistanceModel::Uncorrected,
    ///     output.path(),
    /// );
    /// handle.compute(&files, None);
    /// assert!(output.path().join("per_locus").join("gene_2.csv").exists());
    /// assert!(output.path().join("distance_summary.csv").exists());
    /// ```
    pub fn compute(&self, files: &[PathBuf], prefix: Option<&str>) {
        assert!(
            self.model.is_valid_for(self.datatype),
            "The {} model is not valid for {} data",
            self.model,
            self.datatype
        );
        let prefix = prefix.unwrap_or(DEFAULT_PREFIX);
        let spin = utils::set_spinner();
        spin.set_message("Computing pairwise distances...");
        let mut results = self.par_compute(files);
        results.sort_by(|a, b| alphanumeric_sort::compare_path(&a.summary.path, &b.summary.path));
        spin.set_message("Writing summary...");
        let summaries: Vec<&LocusDistance> = results.iter().map(|r| &r.summary).collect();
        let mut outputs = vec![("Summary", self.write_summary(&summaries, prefix))];
        if self.grouping.is_some() {
            spin.set_message("Writing group mean distances...");
            outputs.extend(
                self.write_group_outputs(&results, prefix)
                    .expect("Failed to write group mean distances"),
            );
        }
        if self.long_format {
            spin.set_message("Writing long-format pairwise distances...");
            let path = self
                .write_pairwise_long(&results, prefix)
                .expect("Failed to write long-format pairwise distances");
            outputs.push(("Pairwise long", path));
        }
        spin.finish_with_message("Finished computing pairwise distances!\n");
        self.log_unassigned(&results);
        self.print_output_info(results.len(), &outputs);
    }

    fn par_compute(&self, files: &[PathBuf]) -> Vec<LocusResult> {
        let matrix_dir = self.output_dir.join(PER_LOCUS_DIR);
        fs::create_dir_all(&matrix_dir).expect("Failed to create output directory");
        if self.grouping.is_some() {
            fs::create_dir_all(self.output_dir.join(PER_LOCUS_GROUP_DIR))
                .expect("Failed to create output directory");
        }
        if self.long_format {
            fs::create_dir_all(self.output_dir.join(LONG_TMP_DIR))
                .expect("Failed to create temporary directory");
        }
        let (tx, rx) = mpsc::channel();
        files
            .par_iter()
            .enumerate()
            .for_each_with(tx, |tx, (index, file)| {
                let result = self.compute_locus(index, file, &matrix_dir);
                tx.send(result).expect("Failed to send summary");
            });
        rx.iter().collect()
    }

    fn compute_locus(&self, index: usize, file: &Path, matrix_dir: &Path) -> LocusResult {
        let (matrix, header) = SeqParser::new(file, self.datatype).get_alignment(self.input_fmt);
        let distances = self.compute_distances(&matrix);
        let locus = self.get_locus_id(file);
        let output_path = matrix_dir.join(&locus).with_extension("csv");
        self.write_matrix(&output_path, &matrix, &distances)
            .unwrap_or_else(|e| {
                log::error!(
                    "Failed to write distance matrix for {}: {}",
                    file.display(),
                    e
                )
            });
        let ids: Vec<String> = matrix.keys().cloned().collect();
        let assignment = self.grouping.as_ref().map(|g| g.assign(&ids));
        let groups = assignment
            .as_ref()
            .map(|a| self.compute_group_distances(a, &distances));
        if let Some(groups) = &groups {
            let path = self
                .output_dir
                .join(PER_LOCUS_GROUP_DIR)
                .join(&locus)
                .with_extension("csv");
            self.write_group_mean_matrix(&path, groups)
                .unwrap_or_else(|e| {
                    log::error!(
                        "Failed to write group distance matrix for {}: {}",
                        file.display(),
                        e
                    )
                });
        }
        let long_path = if self.long_format {
            let path = self.long_tmp_path(index);
            self.write_locus_long(&path, &locus, &ids, assignment.as_ref(), &distances)
                .unwrap_or_else(|e| {
                    log::error!(
                        "Failed to write long-format distances for {}: {}",
                        file.display(),
                        e
                    )
                });
            Some(path)
        } else {
            None
        };
        LocusResult {
            summary: LocusDistance::new(file, self.model, header.nchar, &distances),
            unassigned: assignment.map(|a| a.unassigned()).unwrap_or_default(),
            groups,
            long_path,
        }
    }

    fn compute_distances(&self, matrix: &SeqMatrix) -> PairDistances {
        let seqs: Vec<Vec<u8>> = matrix
            .values()
            .map(|seq| distance::encode(seq, self.datatype))
            .collect();
        let freqs = match self.model {
            DistanceModel::TamuraNei => BaseFreqs::from_encoded(&seqs),
            _ => BaseFreqs::default(),
        };
        // Upper triangle of the matrix, row by row
        let pairs = (0..seqs.len())
            .into_par_iter()
            .map(|i| {
                (i + 1..seqs.len())
                    .map(|j| {
                        let counts = PairCounts::new(&seqs[i], &seqs[j], self.datatype);
                        let distance = self.model.distance(&counts, self.datatype, &freqs);
                        (counts.sites, distance)
                    })
                    .collect::<Vec<_>>()
            })
            .collect();
        PairDistances {
            ntax: seqs.len(),
            pairs,
        }
    }

    fn write_matrix(
        &self,
        output_path: &Path,
        matrix: &SeqMatrix,
        distances: &PairDistances,
    ) -> Result<()> {
        let mut writer = csv::Writer::from_path(output_path)
            .with_context(|| format!("Failed to create {}", output_path.display()))?;
        let ids: Vec<&String> = matrix.keys().collect();
        let mut header = vec![""];
        header.extend(ids.iter().map(|id| id.as_str()));
        writer.write_record(&header)?;
        for (i, id) in ids.iter().enumerate() {
            let mut row = vec![id.to_string()];
            row.extend((0..ids.len()).map(|j| {
                if i == j {
                    self.fmt_diagonal(distances, i)
                } else if self.triangle.includes(i, j) {
                    fmt_distance(distances.get(i, j))
                } else {
                    String::new()
                }
            }));
            writer.write_record(&row)?;
        }
        writer.flush()?;
        Ok(())
    }

    fn fmt_diagonal(&self, distances: &PairDistances, i: usize) -> String {
        match self.diagonal {
            MatrixDiagonal::Zero => fmt_distance(Some(0.0)),
            MatrixDiagonal::Blank => String::new(),
            MatrixDiagonal::RowMax => fmt_distance(distances.row_max(i)),
        }
    }

    fn compute_group_distances(
        &self,
        assignment: &GroupAssignment,
        distances: &PairDistances,
    ) -> GroupDistances {
        let mut groups = GroupDistances::default();
        assignment
            .groups
            .iter()
            .flatten()
            .for_each(|g| groups.add_group(g));
        for i in 0..distances.ntax {
            for j in i + 1..distances.ntax {
                let pair = (&assignment.groups[i], &assignment.groups[j]);
                if let ((Some(a), Some(b)), Some(d)) = (pair, distances.get(i, j)) {
                    groups.add(a, b, d);
                }
            }
        }
        groups
    }

    /// Write a group by group matrix of mean distances.
    /// The diagonal follows [`GroupDiagonal`].
    fn write_group_mean_matrix(&self, output_path: &Path, groups: &GroupDistances) -> Result<()> {
        self.write_group_matrix(output_path, groups, |a, b| {
            if a == b {
                self.fmt_group_diagonal(groups, a)
            } else {
                fmt_distance(groups.get(a, b).mean())
            }
        })
    }

    fn fmt_group_diagonal(&self, groups: &GroupDistances, group: &str) -> String {
        match self.group_diagonal {
            GroupDiagonal::Within => fmt_distance(groups.get(group, group).mean()),
            GroupDiagonal::Zero => fmt_distance(Some(0.0)),
            GroupDiagonal::Blank => String::new(),
            GroupDiagonal::RowMax => fmt_distance(groups.row_max(group)),
        }
    }

    /// Write a group by group matrix with `cell` giving the value for each pair of groups.
    fn write_group_matrix(
        &self,
        output_path: &Path,
        groups: &GroupDistances,
        cell: impl Fn(&str, &str) -> String,
    ) -> Result<()> {
        let mut writer = csv::Writer::from_path(output_path)
            .with_context(|| format!("Failed to create {}", output_path.display()))?;
        let names = groups.names();
        let mut header = vec![""];
        header.extend(names.iter());
        writer.write_record(&header)?;
        for (i, a) in names.iter().enumerate() {
            let mut row = vec![a.to_string()];
            row.extend(names.iter().enumerate().map(|(j, b)| {
                if self.triangle.includes(i, j) {
                    cell(a, b)
                } else {
                    String::new()
                }
            }));
            writer.write_record(&row)?;
        }
        writer.flush()?;
        Ok(())
    }

    fn write_group_outputs(
        &self,
        results: &[LocusResult],
        prefix: &str,
    ) -> Result<Vec<(&'static str, PathBuf)>> {
        let mut pooled = GroupDistances::default();
        results
            .iter()
            .filter_map(|r| r.groups.as_ref())
            .for_each(|g| pooled.merge(g));
        let mean_path = self.prefixed_path(prefix, "group_mean");
        self.write_group_mean_matrix(&mean_path, &pooled)?;
        let pairs_path = self.prefixed_path(prefix, "group_pairs");
        self.write_group_matrix(&pairs_path, &pooled, |a, b| groups_pairs(&pooled, a, b))?;
        let mut outputs = vec![("Group mean", mean_path), ("Group pairs", pairs_path)];
        if self.long_format {
            let long_path = self.prefixed_path(prefix, "group_long");
            self.write_group_long(&long_path, results, &pooled)?;
            outputs.push(("Group long", long_path));
        }
        Ok(outputs)
    }

    fn write_group_long(
        &self,
        output_path: &Path,
        results: &[LocusResult],
        pooled: &GroupDistances,
    ) -> Result<()> {
        let mut writer = csv::Writer::from_path(output_path)
            .with_context(|| format!("Failed to create {}", output_path.display()))?;
        writer.write_record([
            "locus", "group_1", "group_2", "n_pairs", "mean", "sd", "min", "max",
        ])?;
        let loci = results
            .iter()
            .filter_map(|r| Some((r.summary.locus.as_str(), r.groups.as_ref()?)));
        for (locus, groups) in loci.chain(std::iter::once((ALL_LOCI, pooled))) {
            let names = groups.names();
            for (i, a) in names.iter().enumerate() {
                for b in &names[i..] {
                    let stats = groups.get(a, b);
                    writer.write_record([
                        locus.to_string(),
                        a.to_string(),
                        b.to_string(),
                        stats.n_pairs.to_string(),
                        fmt_distance(stats.mean()),
                        fmt_distance(stats.sd()),
                        fmt_distance(stats.min),
                        fmt_distance(stats.max),
                    ])?;
                }
            }
        }
        writer.flush()?;
        Ok(())
    }

    /// Write the pairwise distances of a locus as long-format rows without a header.
    fn write_locus_long(
        &self,
        output_path: &Path,
        locus: &str,
        ids: &[String],
        assignment: Option<&GroupAssignment>,
        distances: &PairDistances,
    ) -> Result<()> {
        let mut writer = csv::WriterBuilder::new()
            .has_headers(false)
            .from_path(output_path)
            .with_context(|| format!("Failed to create {}", output_path.display()))?;
        for i in 0..distances.ntax {
            for j in i + 1..distances.ntax {
                let (sites, distance) = distances.pairs[i][j - i - 1];
                let mut row = vec![locus.to_string(), ids[i].to_string(), ids[j].to_string()];
                if let Some(assignment) = assignment {
                    row.push(assignment.groups[i].clone().unwrap_or_default());
                    row.push(assignment.groups[j].clone().unwrap_or_default());
                }
                row.push(sites.to_string());
                row.push(fmt_distance(distance));
                writer.write_record(&row)?;
            }
        }
        writer.flush()?;
        Ok(())
    }

    /// Join the per-locus long-format rows in locus order and remove the temporary files.
    fn write_pairwise_long(&self, results: &[LocusResult], prefix: &str) -> Result<PathBuf> {
        let output_path = self.prefixed_path(prefix, "pairwise_long");
        let file = File::create(&output_path)
            .with_context(|| format!("Failed to create {}", output_path.display()))?;
        let mut writer = BufWriter::new(file);
        let mut header = vec!["locus", "sample_1", "sample_2"];
        if self.grouping.is_some() {
            header.extend(["group_1", "group_2"]);
        }
        header.extend(["compared_sites", "distance"]);
        writeln!(writer, "{}", header.join(","))?;
        for path in results.iter().filter_map(|r| r.long_path.as_ref()) {
            let mut rows =
                File::open(path).with_context(|| format!("Failed to open {}", path.display()))?;
            io::copy(&mut rows, &mut writer)?;
        }
        writer.flush()?;
        fs::remove_dir_all(self.output_dir.join(LONG_TMP_DIR))
            .context("Failed to remove temporary long-format files")?;
        Ok(output_path)
    }

    fn long_tmp_path(&self, index: usize) -> PathBuf {
        self.output_dir
            .join(LONG_TMP_DIR)
            .join(index.to_string())
            .with_extension("csv")
    }

    fn prefixed_path(&self, prefix: &str, name: &str) -> PathBuf {
        self.output_dir
            .join(format!("{}_{}", prefix, name))
            .with_extension("csv")
    }

    fn log_unassigned(&self, results: &[LocusResult]) {
        let unassigned: HashSet<&str> = results
            .iter()
            .flat_map(|r| r.unassigned.iter().map(|id| id.as_str()))
            .collect();
        if !unassigned.is_empty() {
            log::warn!(
                "{} samples match no group and are left out of group means",
                utils::fmt_num(&unassigned.len())
            );
        }
    }

    fn write_summary(&self, summary: &[&LocusDistance], prefix: &str) -> PathBuf {
        let output_path = self.prefixed_path(prefix, "summary");
        let mut writer = csv::Writer::from_path(&output_path).expect("Failed to create CSV writer");
        summary.iter().for_each(|s| {
            writer.serialize(s).expect("Failed to write summary");
        });
        writer.flush().expect("Failed to write summary");
        output_path
    }

    fn get_locus_id(&self, file: &Path) -> String {
        file.file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string()
    }

    fn print_output_info(&self, locus_counts: usize, outputs: &[(&str, PathBuf)]) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Distance model", self.model);
        log::info!("{:18}: {}", "Locus counts", utils::fmt_num(&locus_counts));
        log::info!("{:18}: {}", "Output dir", self.output_dir.display());
        outputs
            .iter()
            .for_each(|(label, path)| log::info!("{:18}: {}", label, path.display()));
        log::info!("");
    }
}

fn groups_pairs(groups: &GroupDistances, a: &str, b: &str) -> String {
    groups.get(a, b).n_pairs.to_string()
}

fn fmt_distance(distance: Option<f64>) -> String {
    match distance {
        Some(d) => format!("{:.6}", d),
        None => String::from("NA"),
    }
}

/// Pairwise distances of a locus stored as the upper triangle of the matrix.
struct PairDistances {
    ntax: usize,
    /// `pairs[i][j - i - 1]` holds (compared sites, distance) for `i < j`.
    pairs: Vec<Vec<(usize, Option<f64>)>>,
}

impl PairDistances {
    fn get(&self, i: usize, j: usize) -> Option<f64> {
        match i.cmp(&j) {
            std::cmp::Ordering::Equal => Some(0.0),
            std::cmp::Ordering::Less => self.pairs[i][j - i - 1].1,
            std::cmp::Ordering::Greater => self.pairs[j][i - j - 1].1,
        }
    }

    /// The highest defined distance between sample `i` and any other sample.
    fn row_max(&self, i: usize) -> Option<f64> {
        (0..self.ntax)
            .filter(|j| *j != i)
            .filter_map(|j| self.get(i, j))
            .max_by(|a, b| a.total_cmp(b))
    }

    fn iter(&self) -> impl Iterator<Item = &(usize, Option<f64>)> {
        self.pairs.iter().flatten()
    }
}

/// Everything computed for a locus.
struct LocusResult {
    summary: LocusDistance,
    /// Samples without a group. Empty without grouping.
    unassigned: Vec<String>,
    groups: Option<GroupDistances>,
    /// Temporary file with the long-format rows of the locus.
    long_path: Option<PathBuf>,
}

/// Summary of pairwise distances for a locus.
#[derive(Debug, Serialize)]
struct LocusDistance {
    locus: String,
    path: PathBuf,
    model: String,
    ntax: usize,
    nchar: usize,
    pairs: usize,
    undefined_pairs: usize,
    min: Option<f64>,
    max: Option<f64>,
    mean: Option<f64>,
    median: Option<f64>,
    sd: Option<f64>,
    mean_compared_sites: Option<f64>,
}

impl LocusDistance {
    fn new(path: &Path, model: &DistanceModel, nchar: usize, distances: &PairDistances) -> Self {
        let values: Vec<f64> = distances.iter().filter_map(|(_, d)| *d).collect();
        let pairs = distances.iter().count();
        let total_sites: usize = distances.iter().map(|(sites, _)| sites).sum();
        let mean_compared_sites = if pairs > 0 {
            Some(total_sites as f64 / pairs as f64)
        } else {
            None
        };
        let summary = DistanceSummary::new(&values);
        Self {
            locus: path
                .file_stem()
                .unwrap_or_default()
                .to_string_lossy()
                .to_string(),
            path: path.to_path_buf(),
            model: model.to_string(),
            ntax: distances.ntax,
            nchar,
            pairs,
            undefined_pairs: pairs - values.len(),
            min: summary.min,
            max: summary.max,
            mean: summary.mean,
            median: summary.median,
            sd: summary.sd,
            mean_compared_sites,
        }
    }
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;
    use tempdir::TempDir;

    use super::*;

    /// Write a FASTA alignment to `dir` and return its path.
    /// ABCD vs ABCE: 1 diff over 6 sites. ABCF has a gap, so
    /// ABCD vs ABCF: 1 diff over 5 sites.
    fn write_test_locus(dir: &Path) -> PathBuf {
        let path = dir.join("locus_1.fas");
        fs::write(&path, ">ABCD\nACGTAC\n>ABCE\nACGTAT\n>ABCF\nAC-TTC\n").unwrap();
        path
    }

    /// Write two loci for two species. Locus 2 has a single sp2 sample.
    /// Locus 1: within sp1 1/6, within sp2 1/6, between (1 + 2 + 2 + 3) / 24.
    /// Locus 2: within sp1 1/4, between (2 + 3) / 8.
    fn write_group_loci(dir: &Path) -> Vec<PathBuf> {
        let locus_1 = dir.join("locus_1.fas");
        fs::write(
            &locus_1,
            ">sp1_a\nACGTAC\n>sp1_b\nACGTAT\n>sp2_a\nTCGTAC\n>sp2_b\nTCGAAC\n",
        )
        .unwrap();
        let locus_2 = dir.join("locus_2.fas");
        fs::write(&locus_2, ">sp1_a\nAAAA\n>sp1_b\nAAAT\n>sp2_a\nTTAA\n").unwrap();
        vec![locus_2, locus_1]
    }

    fn compute_grouped(output: &Path, long_format: bool) {
        let input = TempDir::new("temp").unwrap();
        let files = write_group_loci(input.path());
        AlignmentDistance::new(
            &InputFmt::Fasta,
            &DataType::Dna,
            &DistanceModel::Uncorrected,
            output,
        )
        .with_grouping(SampleGrouping::from_delimiter('_', 1).unwrap())
        .with_long_format(long_format)
        .compute(&files, None);
    }

    fn read_lines(path: &Path) -> Vec<String> {
        fs::read_to_string(path)
            .unwrap()
            .lines()
            .map(|l| l.to_string())
            .collect()
    }

    #[test]
    fn test_group_means() {
        let output = TempDir::new("temp").unwrap();
        compute_grouped(output.path(), false);
        let locus_1 = output.path().join(PER_LOCUS_GROUP_DIR).join("locus_1.csv");
        assert_eq!(
            read_lines(&locus_1),
            [",sp1,sp2", "sp1,0.166667,0.333333", "sp2,0.333333,0.166667"]
        );
        let locus_2 = output.path().join(PER_LOCUS_GROUP_DIR).join("locus_2.csv");
        assert_eq!(
            read_lines(&locus_2),
            [",sp1,sp2", "sp1,0.250000,0.625000", "sp2,0.625000,NA"]
        );
        // Pooled over all pairs: within sp1 (1/6 + 1/4) / 2,
        // between (8/6 + 5/4) / 6, and within sp2 a single pair of 1/6.
        let pooled = output.path().join("distance_group_mean.csv");
        assert_eq!(
            read_lines(&pooled),
            [",sp1,sp2", "sp1,0.208333,0.430556", "sp2,0.430556,0.166667"]
        );
        let pairs = output.path().join("distance_group_pairs.csv");
        assert_eq!(read_lines(&pairs), [",sp1,sp2", "sp1,2,6", "sp2,6,1"]);
        assert!(!output.path().join("distance_pairwise_long.csv").exists());
        assert!(!output.path().join("distance_group_long.csv").exists());
    }

    #[test]
    fn test_group_matrix_triangle() {
        let input = TempDir::new("temp").unwrap();
        let output = TempDir::new("temp").unwrap();
        let files = write_group_loci(input.path());
        AlignmentDistance::new(
            &InputFmt::Fasta,
            &DataType::Dna,
            &DistanceModel::Uncorrected,
            output.path(),
        )
        .with_triangle(MatrixTriangle::Lower)
        .with_diagonal(MatrixDiagonal::Blank)
        .with_grouping(SampleGrouping::from_delimiter('_', 1).unwrap())
        .compute(&files, Some("test"));
        let pooled = output.path().join("test_group_mean.csv");
        assert_eq!(
            read_lines(&pooled),
            [",sp1,sp2", "sp1,0.208333,", "sp2,0.430556,0.166667"]
        );
    }

    #[test]
    fn test_group_diagonal() {
        // Locus 1 means: within sp1 and sp2 1/6, between 1/3.
        // Pooled: within sp1 0.208333, within sp2 0.166667, between 0.430556.
        for (diagonal, expected) in [
            (
                GroupDiagonal::Zero,
                ["sp1,0.000000,0.430556", "sp2,0.430556,0.000000"],
            ),
            (GroupDiagonal::Blank, ["sp1,,0.430556", "sp2,0.430556,"]),
            (
                GroupDiagonal::RowMax,
                ["sp1,0.430556,0.430556", "sp2,0.430556,0.430556"],
            ),
        ] {
            let input = TempDir::new("temp").unwrap();
            let output = TempDir::new("temp").unwrap();
            let files = write_group_loci(input.path());
            AlignmentDistance::new(
                &InputFmt::Fasta,
                &DataType::Dna,
                &DistanceModel::Uncorrected,
                output.path(),
            )
            .with_grouping(SampleGrouping::from_delimiter('_', 1).unwrap())
            .with_group_diagonal(diagonal)
            .compute(&files, None);
            let pooled = read_lines(&output.path().join("distance_group_mean.csv"));
            assert_eq!(pooled[1..], expected, "{:?}", diagonal);
            // Pair counts always keep the within-group counts on the diagonal.
            let pairs = read_lines(&output.path().join("distance_group_pairs.csv"));
            assert_eq!(pairs[1..], ["sp1,2,6", "sp2,6,1"]);
        }
    }

    #[test]
    fn test_group_row_max_without_other_groups() {
        let input = TempDir::new("temp").unwrap();
        let output = TempDir::new("temp").unwrap();
        let file = write_test_locus(input.path());
        AlignmentDistance::new(
            &InputFmt::Fasta,
            &DataType::Dna,
            &DistanceModel::Uncorrected,
            output.path(),
        )
        .with_grouping(SampleGrouping::from_regex("^ABC").unwrap())
        .with_group_diagonal(GroupDiagonal::RowMax)
        .compute(&[file], None);
        let pooled = read_lines(&output.path().join("distance_group_mean.csv"));
        assert_eq!(pooled, [",ABC", "ABC,NA"]);
    }

    #[test]
    fn test_group_diagonal_from_str() {
        assert_eq!("within".parse::<GroupDiagonal>(), Ok(GroupDiagonal::Within));
        assert_eq!("max".parse::<GroupDiagonal>(), Ok(GroupDiagonal::RowMax));
        assert!("mean".parse::<GroupDiagonal>().is_err());
    }

    #[test]
    fn test_long_format() {
        let output = TempDir::new("temp").unwrap();
        compute_grouped(output.path(), true);
        let lines = read_lines(&output.path().join("distance_pairwise_long.csv"));
        assert_eq!(lines.len(), 10);
        assert_eq!(
            lines[0],
            "locus,sample_1,sample_2,group_1,group_2,compared_sites,distance"
        );
        assert_eq!(lines[1], "locus_1,sp1_a,sp1_b,sp1,sp1,6,0.166667");
        assert_eq!(lines[6], "locus_1,sp2_a,sp2_b,sp2,sp2,6,0.166667");
        assert_eq!(lines[9], "locus_2,sp1_b,sp2_a,sp1,sp2,4,0.750000");

        let lines = read_lines(&output.path().join("distance_group_long.csv"));
        assert_eq!(lines.len(), 10);
        assert_eq!(lines[0], "locus,group_1,group_2,n_pairs,mean,sd,min,max");
        assert_eq!(lines[6], "locus_2,sp2,sp2,0,NA,NA,NA,NA");
        assert!(lines[8].starts_with("all,sp1,sp2,6,0.430556,"));
        assert!(!output.path().join(LONG_TMP_DIR).exists());
    }

    #[test]
    fn test_long_format_without_grouping() {
        let input = TempDir::new("temp").unwrap();
        let output = TempDir::new("temp").unwrap();
        let file = write_test_locus(input.path());
        AlignmentDistance::new(
            &InputFmt::Fasta,
            &DataType::Dna,
            &DistanceModel::Uncorrected,
            output.path(),
        )
        .with_long_format(true)
        .compute(&[file], None);
        let lines = read_lines(&output.path().join("distance_pairwise_long.csv"));
        assert_eq!(
            lines,
            [
                "locus,sample_1,sample_2,compared_sites,distance",
                "locus_1,ABCD,ABCE,6,0.166667",
                "locus_1,ABCD,ABCF,5,0.200000",
                "locus_1,ABCE,ABCF,5,0.400000",
            ]
        );
        assert!(!output.path().join("distance_group_long.csv").exists());
        assert!(!output.path().join(PER_LOCUS_GROUP_DIR).exists());
    }

    #[test]
    fn test_unassigned_samples_in_long_format() {
        let input = TempDir::new("temp").unwrap();
        let output = TempDir::new("temp").unwrap();
        let file = write_test_locus(input.path());
        AlignmentDistance::new(
            &InputFmt::Fasta,
            &DataType::Dna,
            &DistanceModel::Uncorrected,
            output.path(),
        )
        .with_grouping(SampleGrouping::from_regex("^ABC[DE]").unwrap())
        .with_long_format(true)
        .compute(&[file], None);
        let lines = read_lines(&output.path().join("distance_pairwise_long.csv"));
        assert_eq!(lines[2], "locus_1,ABCD,ABCF,ABCD,,5,0.200000");
        let pooled = output.path().join("distance_group_mean.csv");
        assert_eq!(
            read_lines(&pooled),
            [",ABCD,ABCE", "ABCD,NA,0.166667", "ABCE,0.166667,NA"]
        );
    }

    #[test]
    fn test_compute_distances() {
        let input = TempDir::new("temp").unwrap();
        let file = write_test_locus(input.path());
        let (matrix, _) = SeqParser::new(&file, &DataType::Dna).get_alignment(&InputFmt::Fasta);
        let handle = AlignmentDistance::new(
            &InputFmt::Fasta,
            &DataType::Dna,
            &DistanceModel::Uncorrected,
            Path::new("."),
        );
        let distances = handle.compute_distances(&matrix);
        assert_eq!(distances.ntax, 3);
        assert_eq!(distances.get(0, 0), Some(0.0));
        assert_approx_eq!(distances.get(0, 1).unwrap(), 1.0 / 6.0);
        assert_approx_eq!(distances.get(0, 2).unwrap(), 1.0 / 5.0);
        assert_eq!(distances.get(2, 0), distances.get(0, 2));
        assert_eq!(distances.pairs[0][1].0, 5);
    }

    #[test]
    fn test_compute_output() {
        let input = TempDir::new("temp").unwrap();
        let output = TempDir::new("temp").unwrap();
        let file = write_test_locus(input.path());
        let handle = AlignmentDistance::new(
            &InputFmt::Fasta,
            &DataType::Dna,
            &DistanceModel::Uncorrected,
            output.path(),
        );
        handle.compute(&[file], Some("test"));

        let matrix_path = output.path().join(PER_LOCUS_DIR).join("locus_1.csv");
        let content = fs::read_to_string(matrix_path).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines[0], ",ABCD,ABCE,ABCF");
        assert_eq!(lines[1], "ABCD,0.000000,0.166667,0.200000");
        assert_eq!(lines[3], "ABCF,0.200000,0.400000,0.000000");

        let summary_path = output.path().join("test_summary.csv");
        let mut reader = csv::Reader::from_path(summary_path).unwrap();
        let headers = reader.headers().unwrap().clone();
        let record = reader.records().next().unwrap().unwrap();
        let field = |name: &str| {
            let idx = headers.iter().position(|h| h == name).unwrap();
            record.get(idx).unwrap().to_string()
        };
        assert_eq!(field("locus"), "locus_1");
        assert_eq!(field("model"), "p-distance");
        assert_eq!(field("ntax"), "3");
        assert_eq!(field("nchar"), "6");
        assert_eq!(field("pairs"), "3");
        assert_eq!(field("undefined_pairs"), "0");
        assert_approx_eq!(field("max").parse::<f64>().unwrap(), 0.4);
        assert_approx_eq!(field("median").parse::<f64>().unwrap(), 0.2);
    }

    #[test]
    fn test_matrix_triangle() {
        let input = TempDir::new("temp").unwrap();
        let file = write_test_locus(input.path());
        for (triangle, expected) in [
            (
                MatrixTriangle::Lower,
                [
                    "ABCD,0.000000,,",
                    "ABCE,0.166667,0.000000,",
                    "ABCF,0.200000,0.400000,0.000000",
                ],
            ),
            (
                MatrixTriangle::Upper,
                [
                    "ABCD,0.000000,0.166667,0.200000",
                    "ABCE,,0.000000,0.400000",
                    "ABCF,,,0.000000",
                ],
            ),
        ] {
            let output = TempDir::new("temp").unwrap();
            AlignmentDistance::new(
                &InputFmt::Fasta,
                &DataType::Dna,
                &DistanceModel::Uncorrected,
                output.path(),
            )
            .with_triangle(triangle)
            .compute(std::slice::from_ref(&file), None);
            let matrix_path = output.path().join(PER_LOCUS_DIR).join("locus_1.csv");
            let content = fs::read_to_string(matrix_path).unwrap();
            let lines: Vec<&str> = content.lines().skip(1).collect();
            assert_eq!(lines, expected, "{:?}", triangle);
        }
    }

    #[test]
    fn test_matrix_diagonal() {
        let input = TempDir::new("temp").unwrap();
        let file = write_test_locus(input.path());
        for (diagonal, expected) in [
            (
                MatrixDiagonal::Blank,
                [",,", "0.166667,,", "0.200000,0.400000,"],
            ),
            (
                MatrixDiagonal::RowMax,
                [
                    "0.200000,,",
                    "0.166667,0.400000,",
                    "0.200000,0.400000,0.400000",
                ],
            ),
        ] {
            let output = TempDir::new("temp").unwrap();
            AlignmentDistance::new(
                &InputFmt::Fasta,
                &DataType::Dna,
                &DistanceModel::Uncorrected,
                output.path(),
            )
            .with_triangle(MatrixTriangle::Lower)
            .with_diagonal(diagonal)
            .compute(std::slice::from_ref(&file), None);
            let matrix_path = output.path().join(PER_LOCUS_DIR).join("locus_1.csv");
            let content = fs::read_to_string(matrix_path).unwrap();
            // Drop the header row and the sample ID column
            let lines: Vec<&str> = content
                .lines()
                .skip(1)
                .map(|line| line.split_once(',').unwrap().1)
                .collect();
            assert_eq!(lines, expected, "{:?}", diagonal);
        }
    }

    #[test]
    fn test_row_max_without_pairs() {
        let distances = PairDistances {
            ntax: 2,
            pairs: vec![vec![(0, None)], vec![]],
        };
        assert_eq!(distances.row_max(0), None);
        assert_eq!(fmt_distance(distances.row_max(0)), "NA");
    }

    #[test]
    fn test_matrix_diagonal_from_str() {
        assert_eq!("max".parse::<MatrixDiagonal>(), Ok(MatrixDiagonal::RowMax));
        assert!("one".parse::<MatrixDiagonal>().is_err());
    }

    #[test]
    fn test_matrix_triangle_from_str() {
        assert_eq!("lower".parse::<MatrixTriangle>(), Ok(MatrixTriangle::Lower));
        assert!("diagonal".parse::<MatrixTriangle>().is_err());
    }

    #[test]
    fn test_single_taxon_locus() {
        let output = TempDir::new("temp").unwrap();
        let file = PathBuf::from("tests/files/alignments/gene_1.nex");
        let handle = AlignmentDistance::new(
            &InputFmt::Nexus,
            &DataType::Dna,
            &DistanceModel::Kimura2P,
            output.path(),
        );
        let results = handle.par_compute(&[file]);
        assert_eq!(results[0].summary.pairs, 0);
        assert_eq!(results[0].summary.mean, None);
    }

    #[test]
    #[should_panic(expected = "not valid")]
    fn test_invalid_model() {
        let output = TempDir::new("temp").unwrap();
        AlignmentDistance::new(
            &InputFmt::Nexus,
            &DataType::Aa,
            &DistanceModel::Kimura2P,
            output.path(),
        )
        .compute(&[], None);
    }
}
