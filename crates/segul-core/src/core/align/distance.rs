//! Compute pairwise genetic distances for each locus alignment.
//!
//! Each alignment gets its own square distance matrix,
//! written to `output_dir/per_locus/<locus>.csv`.
//! A summary of each locus's pairwise distances
//! is written to `output_dir/<prefix>_summary.csv`.
use std::fs;
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

const PER_LOCUS_DIR: &str = "per_locus";
const DEFAULT_PREFIX: &str = "distance";

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
        }
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
        let spin = utils::set_spinner();
        spin.set_message("Computing pairwise distances...");
        let mut summary = self.par_compute(files);
        summary.sort_by(|a, b| alphanumeric_sort::compare_path(&a.path, &b.path));
        spin.set_message("Writing summary...");
        let summary_path = self.write_summary(&summary, prefix);
        spin.finish_with_message("Finished computing pairwise distances!\n");
        self.print_output_info(summary.len(), &summary_path);
    }

    fn par_compute(&self, files: &[PathBuf]) -> Vec<LocusDistance> {
        let matrix_dir = self.output_dir.join(PER_LOCUS_DIR);
        fs::create_dir_all(&matrix_dir).expect("Failed to create output directory");
        let (tx, rx) = mpsc::channel();
        files.par_iter().for_each_with(tx, |tx, file| {
            let (matrix, header) =
                SeqParser::new(file, self.datatype).get_alignment(self.input_fmt);
            let distances = self.compute_distances(&matrix);
            let output_path = matrix_dir
                .join(self.get_locus_id(file))
                .with_extension("csv");
            self.write_matrix(&output_path, &matrix, &distances)
                .unwrap_or_else(|e| {
                    log::error!(
                        "Failed to write distance matrix for {}: {}",
                        file.display(),
                        e
                    )
                });
            let summary = LocusDistance::new(file, self.model, header.nchar, &distances);
            tx.send(summary).expect("Failed to send summary");
        });
        rx.iter().collect()
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
            row.extend((0..ids.len()).map(|j| fmt_distance(distances.get(i, j))));
            writer.write_record(&row)?;
        }
        writer.flush()?;
        Ok(())
    }

    fn write_summary(&self, summary: &[LocusDistance], prefix: Option<&str>) -> PathBuf {
        let fname = format!("{}_summary", prefix.unwrap_or(DEFAULT_PREFIX));
        let output_path = self.output_dir.join(fname).with_extension("csv");
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

    fn print_output_info(&self, locus_counts: usize, summary_path: &Path) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Distance model", self.model);
        log::info!("{:18}: {}", "Locus counts", utils::fmt_num(&locus_counts));
        log::info!("{:18}: {}", "Output dir", self.output_dir.display());
        log::info!("{:18}: {}\n", "Summary", summary_path.display());
    }
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

    fn iter(&self) -> impl Iterator<Item = &(usize, Option<f64>)> {
        self.pairs.iter().flatten()
    }
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
    fn test_single_taxon_locus() {
        let output = TempDir::new("temp").unwrap();
        let file = PathBuf::from("tests/files/alignments/gene_1.nex");
        let handle = AlignmentDistance::new(
            &InputFmt::Nexus,
            &DataType::Dna,
            &DistanceModel::Kimura2P,
            output.path(),
        );
        let summary = handle.par_compute(&[file]);
        assert_eq!(summary[0].pairs, 0);
        assert_eq!(summary[0].mean, None);
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
