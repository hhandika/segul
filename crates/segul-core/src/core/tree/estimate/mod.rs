//! Phylogenetic tree estimation using distance-based (NJ, BIONJ, UPGMA)
//! and Maximum Parsimony methods.

pub mod bootstrap;
pub mod distance;
pub mod matrix;
pub mod parsimony;
pub mod rng;
pub mod rooting;

use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use colored::Colorize;
use serde::Serialize;

use crate::core::tree::convert::output_tree_path;
use crate::helper::sequence::SeqParser;
use crate::helper::types::{DataType, InputFmt, SeqMatrix};
use crate::helper::utils;
use crate::parser::tree::{Rootedness, Tree, TreeOutputFmt};
use crate::report::spec::{TaskKind, TreeInferParams, canonical};
use crate::report::stats::{LocusTreeInferStats, TreeInferStats};
use crate::report::{DirSnapshot, ReportContext, RunReport, RunReportBuilder, report_path};
use crate::stats::distance::DistanceModel;
use crate::writer::tree::{TreeWriteOptions, TreeWriter};

pub use matrix::DistanceMatrix;
pub use parsimony::{
    BranchSwap, GapMode, ParsimonyBranchLengths, ParsimonyOptions, ParsimonyResult,
    ParsimonySearch, estimate_parsimony,
};
pub use rooting::{RootingMethod, reroot_midpoint, reroot_outgroup};

const DEFAULT_PREFIX: &str = "tree_infer";

/// Supported tree estimation methods.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum TreeMethod {
    /// Neighbor-Joining (Saitou & Nei 1987).
    #[default]
    Nj,
    /// Variance-weighted Neighbor-Joining (Gascuel 1997).
    BioNj,
    /// Unweighted Pair Group Method with Arithmetic Mean.
    Upgma,
    /// Maximum Parsimony (Fitch).
    Parsimony,
}

impl TreeMethod {
    /// Canonical lowercase identifier recorded in `report.json`.
    pub fn canonical_str(self) -> &'static str {
        match self {
            Self::Nj => "nj",
            Self::BioNj => "bionj",
            Self::Upgma => "upgma",
            Self::Parsimony => "parsimony",
        }
    }

    /// Return true for distance-based methods (`Nj`, `BioNj`, `Upgma`).
    pub fn is_distance(self) -> bool {
        matches!(self, Self::Nj | Self::BioNj | Self::Upgma)
    }
}

impl std::fmt::Display for TreeMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Nj => write!(f, "NJ"),
            Self::BioNj => write!(f, "BIONJ"),
            Self::Upgma => write!(f, "UPGMA"),
            Self::Parsimony => write!(f, "Maximum Parsimony"),
        }
    }
}

impl std::str::FromStr for TreeMethod {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "nj" => Ok(Self::Nj),
            "bionj" => Ok(Self::BioNj),
            "upgma" => Ok(Self::Upgma),
            "mp" | "parsimony" => Ok(Self::Parsimony),
            _ => Err(format!("{} is not a valid tree estimation method", s)),
        }
    }
}

/// Configuration for [`estimate_matrix`] and [`TreeEstimation`].
#[derive(Debug, Clone)]
pub struct TreeEstimateOptions {
    pub method: TreeMethod,
    pub model: DistanceModel,
    pub parsimony: ParsimonyOptions,
    pub bootstrap: usize,
    pub seed: u64,
    pub rooting: RootingMethod,
    pub output_fmt: TreeOutputFmt,
    pub single_file: bool,
    pub write_opts: TreeWriteOptions,
}

impl Default for TreeEstimateOptions {
    fn default() -> Self {
        Self {
            method: TreeMethod::Nj,
            model: DistanceModel::JukesCantor,
            parsimony: ParsimonyOptions::default(),
            bootstrap: 0,
            seed: 42,
            rooting: RootingMethod::None,
            output_fmt: TreeOutputFmt::Newick,
            single_file: false,
            write_opts: TreeWriteOptions::default(),
        }
    }
}

/// Result of estimating a single phylogenetic tree from an in-memory [`SeqMatrix`].
#[derive(Debug, Clone)]
pub struct TreeEstimateResult {
    pub tree: Tree,
    pub ntax: usize,
    pub nchar: usize,
    pub parsimony_score: Option<usize>,
    pub consistency_index: Option<f64>,
    pub retention_index: Option<f64>,
    pub tree_length: Option<f64>,
    pub undefined_pairs: usize,
}

/// Files written by [`TreeEstimation`].
#[derive(Debug, Clone)]
pub struct TreeEstimateFiles {
    /// Path to the per-locus summary CSV.
    pub summary: PathBuf,
    /// Written tree file paths (one per locus, or a single combined file).
    pub tree_files: Vec<PathBuf>,
    /// Path to unified `report.json`, if written.
    pub report: Option<PathBuf>,
}

#[derive(Debug, Serialize)]
struct CsvSummaryRow {
    locus: String,
    ntax: usize,
    nchar: usize,
    method: String,
    model: String,
    parsimony_score: Option<usize>,
    ci: Option<f64>,
    ri: Option<f64>,
    tree_length: Option<f64>,
    bootstrap_replicates: usize,
    rooted: bool,
    undefined_pairs: usize,
    output_tree: String,
}

/// Estimate a phylogenetic tree from an in-memory alignment [`SeqMatrix`].
pub fn estimate_matrix(
    matrix: &SeqMatrix,
    datatype: &DataType,
    options: &TreeEstimateOptions,
) -> Result<TreeEstimateResult> {
    if matches!(datatype, DataType::Ignore) {
        bail!("Tree estimation requires DNA or amino acid datatype, not 'ignore'");
    }
    let ntax = matrix.len();
    if ntax < 3 {
        bail!(
            "At least 3 taxa are required for tree estimation (found {})",
            ntax
        );
    }
    let nchar = matrix.values().next().map(<[u8]>::len).unwrap_or(0);

    let (mut tree, parsimony_score, ci, ri, undefined_pairs) = match options.method {
        TreeMethod::Nj | TreeMethod::BioNj | TreeMethod::Upgma => {
            let dm = DistanceMatrix::from_seq_matrix(matrix, datatype, &options.model)?;
            if dm.undefined_pairs > 0 {
                log::warn!(
                    "{} pairwise distance(s) were undefined or saturated and imputed",
                    dm.undefined_pairs
                );
            }
            let t = match options.method {
                TreeMethod::Nj => distance::estimate_nj(&dm)?,
                TreeMethod::BioNj => distance::estimate_bionj(&dm)?,
                TreeMethod::Upgma => distance::estimate_upgma(&dm)?,
                TreeMethod::Parsimony => unreachable!(),
            };
            (t, None, None, None, dm.undefined_pairs)
        }
        TreeMethod::Parsimony => {
            let mut p_opts = options.parsimony.clone();
            p_opts.seed = options.seed;
            let res = estimate_parsimony(matrix, datatype, &p_opts)?;
            (
                res.tree,
                Some(res.score),
                res.consistency_index,
                res.retention_index,
                0,
            )
        }
    };

    tree = match &options.rooting {
        RootingMethod::None => tree,
        RootingMethod::Midpoint => reroot_midpoint(&tree)?,
        RootingMethod::Outgroup(outgroup) => reroot_outgroup(&tree, outgroup)?,
    };

    if options.bootstrap > 0 {
        tree = bootstrap::run_bootstrap(&tree, matrix, datatype, options)?;
    }

    let tree_length = compute_tree_length(&tree);

    Ok(TreeEstimateResult {
        tree,
        ntax,
        nchar,
        parsimony_score,
        consistency_index: ci,
        retention_index: ri,
        tree_length,
        undefined_pairs,
    })
}

fn compute_tree_length(tree: &Tree) -> Option<f64> {
    let mut has_any = false;
    let mut sum = 0.0f64;
    for id in tree.preorder() {
        if let Some(len) = tree.branch_length(id)
            && len.is_finite()
        {
            has_any = true;
            sum += len;
        }
    }
    if has_any { Some(sum) } else { None }
}

/// Batch phylogenetic tree estimation task across alignment files.
pub struct TreeEstimation<'a> {
    input_fmt: &'a InputFmt,
    datatype: &'a DataType,
    options: TreeEstimateOptions,
    report_context: ReportContext,
}

impl<'a> TreeEstimation<'a> {
    /// Create a new batch tree estimation task.
    pub fn new(
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        options: TreeEstimateOptions,
    ) -> Self {
        Self {
            input_fmt,
            datatype,
            options,
            report_context: ReportContext::default(),
        }
    }

    /// Attach frontend provenance and SHA-256 checksum settings for `report.json`.
    pub fn with_report_context(mut self, context: ReportContext) -> Self {
        self.report_context = context;
        self
    }

    /// Estimate trees for `files` into `output` without writing `report.json`.
    pub fn estimate(
        &self,
        files: &[PathBuf],
        output: &Path,
        prefix: Option<&str>,
    ) -> Result<TreeEstimateFiles> {
        let (files_out, _) = self.run_inner(files, output, prefix, false)?;
        Ok(files_out)
    }

    /// Estimate trees for `files` into `output` and write a unified `report.json`.
    pub fn estimate_with_report(
        &self,
        files: &[PathBuf],
        output: &Path,
        prefix: Option<&str>,
    ) -> Result<(TreeEstimateFiles, RunReport)> {
        let (files_out, report) = self.run_inner(files, output, prefix, true)?;
        Ok((files_out, report.expect("report requested")))
    }

    fn run_inner(
        &self,
        files: &[PathBuf],
        output: &Path,
        prefix: Option<&str>,
        write_report: bool,
    ) -> Result<(TreeEstimateFiles, Option<RunReport>)> {
        if files.is_empty() {
            bail!("No alignment files provided for tree estimation");
        }
        if matches!(self.datatype, DataType::Ignore) {
            bail!("Tree estimation requires DNA or amino acid datatype, not 'ignore'");
        }
        if self.options.method.is_distance() && !self.options.model.is_valid_for(self.datatype) {
            bail!(
                "The {} model is not valid for {} data. \
                Use p, jc69, k2p, or tn93 for DNA, and p, jc69, or poisson for amino acid.",
                self.options.model,
                self.datatype
            );
        }

        let mut sorted = files.to_vec();
        alphanumeric_sort::sort_path_slice(&mut sorted);
        fs::create_dir_all(output)
            .with_context(|| format!("Failed creating output directory {}", output.display()))?;

        let snapshot = DirSnapshot::capture(output);

        let (rooting_str, outgroup_opt) = match &self.options.rooting {
            RootingMethod::None => ("none".to_string(), None),
            RootingMethod::Midpoint => ("midpoint".to_string(), None),
            RootingMethod::Outgroup(taxa) => ("outgroup".to_string(), Some(taxa.clone())),
        };

        let params = TreeInferParams {
            input_format: canonical::input_fmt(self.input_fmt).to_string(),
            datatype: canonical::datatype(self.datatype).to_string(),
            method: self.options.method.canonical_str().to_string(),
            model: if self.options.method.is_distance() {
                Some(canonical::distance_model(&self.options.model).to_string())
            } else {
                None
            },
            parsimony_search: if !self.options.method.is_distance() {
                Some(self.options.parsimony.search.to_string())
            } else {
                None
            },
            swap: if !self.options.method.is_distance() {
                Some(self.options.parsimony.swap.to_string())
            } else {
                None
            },
            replicates: if !self.options.method.is_distance() {
                Some(self.options.parsimony.replicates)
            } else {
                None
            },
            gap_mode: if !self.options.method.is_distance() {
                Some(self.options.parsimony.gap_mode.to_string())
            } else {
                None
            },
            mp_branch_lengths: if !self.options.method.is_distance() {
                Some(self.options.parsimony.branch_lengths.to_string())
            } else {
                None
            },
            bootstrap: self.options.bootstrap,
            seed: self.options.seed,
            rooting: rooting_str,
            outgroup: outgroup_opt,
            output_format: canonical::tree_output_fmt(&self.options.output_fmt).to_string(),
            single_file: self.options.single_file,
            translate: self.options.write_opts.use_translate,
            gzip: self.options.write_opts.gzip,
            prefix: prefix.map(str::to_string),
        };

        let report_builder = if write_report {
            Some(
                RunReportBuilder::new(TaskKind::TreeInfer, &self.report_context)
                    .primary_inputs(&sorted)
                    .parameters(&params)?,
            )
        } else {
            None
        };

        let spin = utils::set_spinner();
        spin.set_message(format!(
            "Estimating phylogenetic trees ({})...",
            self.options.method
        ));

        let mut estimated = Vec::with_capacity(sorted.len());
        for file in &sorted {
            if !file.is_file() {
                bail!("Input alignment file not found: {}", file.display());
            }
            let (matrix, _header) = SeqParser::new(file, self.datatype).parse(self.input_fmt);
            let mut res = estimate_matrix(&matrix, self.datatype, &self.options)
                .with_context(|| format!("Failed estimating tree for {}", file.display()))?;
            let locus = locus_stem(file);
            res.tree.metadata.name = Some(locus.clone());
            estimated.push((file.clone(), locus, res));
        }

        let mut tree_files = Vec::new();
        let mut csv_rows = Vec::with_capacity(estimated.len());
        let mut stat_entries = Vec::with_capacity(estimated.len());

        let clean_prefix = prefix.filter(|p| !p.is_empty());
        if self.options.single_file {
            let base_name = clean_prefix.unwrap_or("trees");
            let mut filename = format!("{base_name}.{}", self.options.output_fmt.extension());
            if self.options.write_opts.gzip {
                filename.push_str(".gz");
            }
            let out_tree_path = output.join(&filename);
            let trees: Vec<Tree> = estimated.iter().map(|(_, _, r)| r.tree.clone()).collect();
            let writer = TreeWriter::new(
                &out_tree_path,
                self.options.output_fmt,
                self.options.write_opts.clone(),
            );
            writer.write_trees(&trees)?;
            tree_files.push(out_tree_path);

            for (_, locus, res) in &estimated {
                csv_rows.push(self.build_csv_row(locus, res, &filename));
                stat_entries.push(self.build_stat_entry(locus, res));
            }
        } else {
            for (file, locus, res) in &estimated {
                let out_tree_path = output_tree_path(
                    output,
                    file,
                    self.options.output_fmt,
                    self.options.write_opts.gzip,
                );
                let out_name = out_tree_path
                    .file_name()
                    .and_then(|n| n.to_str())
                    .unwrap_or("tree")
                    .to_string();
                let writer = TreeWriter::new(
                    &out_tree_path,
                    self.options.output_fmt,
                    self.options.write_opts.clone(),
                );
                writer.write_trees(std::slice::from_ref(&res.tree))?;
                tree_files.push(out_tree_path);
                csv_rows.push(self.build_csv_row(locus, res, &out_name));
                stat_entries.push(self.build_stat_entry(locus, res));
            }
        }

        let summary_prefix = clean_prefix.unwrap_or(DEFAULT_PREFIX);
        let summary_path = output.join(format!("{summary_prefix}_summary.csv"));
        let mut csv_writer = csv::Writer::from_path(&summary_path)
            .with_context(|| format!("Failed creating {}", summary_path.display()))?;
        for row in &csv_rows {
            csv_writer.serialize(row)?;
        }
        csv_writer.flush()?;

        spin.finish_with_message("Finished estimating phylogenetic trees!\n");
        self.print_output_info(output, estimated.len(), &summary_path);

        let report_out = if let Some(builder) = report_builder {
            let written = snapshot.written_files(output);
            let stats = TreeInferStats {
                files: sorted.len(),
                trees_written: estimated.len(),
                method: self.options.method.canonical_str().to_string(),
                entries: stat_entries,
            };
            let rpath = report_path(output, clean_prefix);
            let report = builder
                .outputs(
                    output,
                    written.into_iter().map(|p| {
                        let role = if p == summary_path { "summary" } else { "tree" };
                        (role, p)
                    }),
                )
                .statistics(&stats)?
                .write_to(&rpath)?;
            Some((rpath, report))
        } else {
            None
        };

        let (report_file, report_obj) = match report_out {
            Some((p, r)) => (Some(p), Some(r)),
            None => (None, None),
        };

        Ok((
            TreeEstimateFiles {
                summary: summary_path,
                tree_files,
                report: report_file,
            },
            report_obj,
        ))
    }

    fn build_csv_row(
        &self,
        locus: &str,
        res: &TreeEstimateResult,
        output_tree: &str,
    ) -> CsvSummaryRow {
        let model = if self.options.method.is_distance() {
            canonical::distance_model(&self.options.model).to_string()
        } else {
            String::new()
        };
        CsvSummaryRow {
            locus: locus.to_string(),
            ntax: res.ntax,
            nchar: res.nchar,
            method: self.options.method.canonical_str().to_string(),
            model,
            parsimony_score: res.parsimony_score,
            ci: res.consistency_index,
            ri: res.retention_index,
            tree_length: res.tree_length,
            bootstrap_replicates: self.options.bootstrap,
            rooted: res.tree.metadata.rootedness == Rootedness::Rooted,
            undefined_pairs: res.undefined_pairs,
            output_tree: output_tree.to_string(),
        }
    }

    fn build_stat_entry(&self, locus: &str, res: &TreeEstimateResult) -> LocusTreeInferStats {
        LocusTreeInferStats {
            locus: locus.to_string(),
            ntax: res.ntax,
            nchar: res.nchar,
            parsimony_score: res.parsimony_score,
            ci: res.consistency_index,
            ri: res.retention_index,
            tree_length: res.tree_length,
            bootstrap_replicates: self.options.bootstrap,
            rooted: res.tree.metadata.rootedness == Rootedness::Rooted,
            undefined_pairs: res.undefined_pairs,
        }
    }

    fn print_output_info(&self, output: &Path, total_trees: usize, summary_path: &Path) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Output dir", output.display());
        log::info!("{:18}: {}", "Method", self.options.method);
        if self.options.method.is_distance() {
            log::info!("{:18}: {}", "Distance model", self.options.model);
        }
        if self.options.bootstrap > 0 {
            log::info!(
                "{:18}: {} (seed={})",
                "Bootstrap",
                utils::fmt_num(&self.options.bootstrap),
                self.options.seed
            );
        }
        log::info!("{:18}: {}", "Output format", self.options.output_fmt);
        log::info!("{:18}: {}", "Trees written", utils::fmt_num(&total_trees));
        log::info!("{:18}: {}\n", "Summary CSV", summary_path.display());
    }
}

fn locus_stem(path: &Path) -> String {
    path.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("tree")
        .to_string()
}
