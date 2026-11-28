//! High-performance phylogenetic summary tree construction and annotation.
//!
//! Supports Maximum Clade Credibility (MCC), Majority-Rule Consensus (MRC),
//! Extended Majority-Rule Consensus (MRE), Strict Consensus, Maximum A Posteriori
//! (MAP), and Target-tree annotation with $O(N \log N)$ Common Ancestor heights,
//! exact Chen-Shao HPD intervals, BEAST/FigTree/NHX metadata summarization, and
//! optional partial-taxa support.

pub(crate) mod annotate;
pub(crate) mod bitset;
pub(crate) mod consensus;
pub(crate) mod engine;
pub(crate) mod stats_acc;

use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use anyhow::{Context, Result, bail};
use colored::Colorize;
use serde::{Deserialize, Serialize};

use crate::core::tree::estimate::RootingMethod;
use crate::helper::progress::{self, ProgressReporter};
use crate::helper::utils;
use crate::parser::tree::{
    ParseOptions, Tree, TreeFmt, TreeOutputFmt, count_trees_file, parse_tree_file,
    read_single_tree_at,
};
use crate::report::spec::{TaskKind, TreeSummaryParams, canonical};
use crate::report::stats::{TreeSummaryFileStats, TreeSummaryStats};
use crate::report::{DirSnapshot, ReportContext, RunReport, RunReportBuilder, report_path};
use crate::writer::tree::{TreeWriteOptions, TreeWriter};

use annotate::{
    annotate_and_finalize_tree, summary_method_label, tree_length_and_height, write_clades_csv,
    write_topologies_csv,
};
use consensus::{
    TreeLocation, build_consensus_tree, map_tree_nodes_to_clades, select_map_tree, select_mcc_tree,
};
use engine::{
    FilePlan, compute_burnin_skip, compute_ca_heights_files, compute_ca_heights_in_memory,
    needs_sample_storage, summarize_files_pass1, summarize_in_memory_pass1,
};

const DEFAULT_PREFIX: &str = "tree_summary";
const MAX_TOPOLOGY_CSV_ROWS: usize = 1000;

/// Summary tree construction method.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum SummaryMethod {
    /// Maximum Clade Credibility tree (default).
    #[default]
    Mcc,
    /// Majority-rule consensus tree.
    Mrc,
    /// Extended (greedy) majority-rule consensus tree.
    Mre,
    /// Strict consensus tree (100% clade agreement).
    Strict,
    /// Maximum A Posteriori topology tree.
    Map,
    /// Map support and node summaries onto a user-supplied target tree.
    Target,
}

impl SummaryMethod {
    /// Canonical lowercase identifier recorded in `report.json`.
    pub fn canonical_str(self) -> &'static str {
        match self {
            Self::Mcc => "mcc",
            Self::Mrc => "mrc",
            Self::Mre => "mre",
            Self::Strict => "strict",
            Self::Map => "map",
            Self::Target => "target",
        }
    }
}

impl fmt::Display for SummaryMethod {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.canonical_str())
    }
}

impl FromStr for SummaryMethod {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "mcc" => Ok(Self::Mcc),
            "mrc" => Ok(Self::Mrc),
            "mre" => Ok(Self::Mre),
            "strict" => Ok(Self::Strict),
            "map" => Ok(Self::Map),
            "target" => Ok(Self::Target),
            _ => Err(format!("'{}' is not a valid tree summary method", s)),
        }
    }
}

/// Strategy for summarizing branch lengths and node heights onto the summary tree.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum BranchLengthStrategy {
    /// Mean observed branch length (and mean node height when rooted).
    #[default]
    Mean,
    /// Median observed branch length (and median node height when rooted).
    Median,
    /// Common Ancestor (CA) node heights for rooted trees, guaranteeing non-negative branches.
    CommonAncestor,
    /// Keep original branch lengths from the selected MCC, MAP, or target tree.
    Keep,
    /// Strip branch lengths from the summary tree.
    None,
}

impl BranchLengthStrategy {
    pub fn canonical_str(self) -> &'static str {
        match self {
            Self::Mean => "mean",
            Self::Median => "median",
            Self::CommonAncestor => "ca",
            Self::Keep => "keep",
            Self::None => "none",
        }
    }
}

impl fmt::Display for BranchLengthStrategy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.canonical_str())
    }
}

impl FromStr for BranchLengthStrategy {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "mean" => Ok(Self::Mean),
            "median" => Ok(Self::Median),
            "ca" | "common-ancestor" | "common_ancestor" => Ok(Self::CommonAncestor),
            "keep" => Ok(Self::Keep),
            "none" => Ok(Self::None),
            _ => Err(format!("'{}' is not a valid branch length strategy", s)),
        }
    }
}

/// Format for node support values.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum SupportFormat {
    /// Proportions in `[0.0, 1.0]` (default).
    #[default]
    Proportion,
    /// Percentages in `[0.0, 100.0]`.
    Percent,
}

impl SupportFormat {
    pub fn canonical_str(self) -> &'static str {
        match self {
            Self::Proportion => "proportion",
            Self::Percent => "percent",
        }
    }
}

impl fmt::Display for SupportFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.canonical_str())
    }
}

impl FromStr for SupportFormat {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "proportion" | "prob" => Ok(Self::Proportion),
            "percent" | "pct" => Ok(Self::Percent),
            _ => Err(format!("'{}' is not a valid support format", s)),
        }
    }
}

/// Where support values are written on the summary tree.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum SupportPlacement {
    /// Annotation-only for NEXUS/PhyloXML, both for Newick/eNewick/NeXML (default).
    #[default]
    Auto,
    /// Write support inside `[&posterior=...]` / `[&support=...]` annotations only.
    Annotation,
    /// Write support as internal node labels only.
    Label,
    /// Write support both as internal node labels and inside annotations.
    Both,
}

impl SupportPlacement {
    pub fn canonical_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Annotation => "annotation",
            Self::Label => "label",
            Self::Both => "both",
        }
    }
}

impl fmt::Display for SupportPlacement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.canonical_str())
    }
}

impl FromStr for SupportPlacement {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "auto" => Ok(Self::Auto),
            "annotation" => Ok(Self::Annotation),
            "label" => Ok(Self::Label),
            "both" => Ok(Self::Both),
            _ => Err(format!("'{}' is not a valid support placement", s)),
        }
    }
}

/// Controls whether input trees are interpreted as rooted clades or unrooted bipartitions.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum RootedMode {
    /// Detect from the first retained tree (`[&R]`/`[&U]` or root degree).
    #[default]
    Auto,
    /// Treat all input trees as rooted.
    Rooted,
    /// Treat all input trees as unrooted bipartitions.
    Unrooted,
}

impl RootedMode {
    pub fn canonical_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Rooted => "rooted",
            Self::Unrooted => "unrooted",
        }
    }
}

impl fmt::Display for RootedMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.canonical_str())
    }
}

impl FromStr for RootedMode {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "auto" => Ok(Self::Auto),
            "rooted" => Ok(Self::Rooted),
            "unrooted" => Ok(Self::Unrooted),
            _ => Err(format!("'{}' is not a valid rooted mode", s)),
        }
    }
}

/// Scoring objective used to rank candidate trees in MCC selection.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum MccScoring {
    /// Maximize sum of log clade probabilities (default, matches TreeAnnotator).
    #[default]
    LogProduct,
    /// Maximize sum of clade probabilities.
    Sum,
}

impl MccScoring {
    pub fn canonical_str(self) -> &'static str {
        match self {
            Self::LogProduct => "log-product",
            Self::Sum => "sum",
        }
    }
}

impl fmt::Display for MccScoring {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.canonical_str())
    }
}

impl FromStr for MccScoring {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "log-product" | "log_product" | "log" | "product" => Ok(Self::LogProduct),
            "sum" => Ok(Self::Sum),
            _ => Err(format!("'{}' is not a valid MCC scoring mode", s)),
        }
    }
}

/// Per-file burn-in specification (fixed tree count or percentage/fraction).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BurninSpec {
    Count(usize),
    Fraction(f64),
}

impl Default for BurninSpec {
    fn default() -> Self {
        Self::Count(0)
    }
}

impl fmt::Display for BurninSpec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Count(n) => write!(f, "{n}"),
            Self::Fraction(frac) => write!(f, "{}%", frac * 100.0),
        }
    }
}

impl FromStr for BurninSpec {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        let trimmed = s.trim();
        if let Some(pct_str) = trimmed.strip_suffix('%') {
            let pct: f64 = pct_str
                .trim()
                .parse()
                .map_err(|_| format!("Invalid burn-in percentage '{}'", s))?;
            if !(0.0..100.0).contains(&pct) {
                return Err(format!(
                    "Burn-in percentage must be in [0%, 100%), got {}",
                    s
                ));
            }
            return Ok(Self::Fraction(pct / 100.0));
        }
        if let Ok(n) = trimmed.parse::<usize>() {
            return Ok(Self::Count(n));
        }
        if let Ok(frac) = trimmed.parse::<f64>()
            && (0.0..1.0).contains(&frac)
        {
            return Ok(Self::Fraction(frac));
        }
        Err(format!(
            "Invalid burn-in '{}': use a non-negative integer count, a percentage like '10%', or a fraction in [0, 1)",
            s
        ))
    }
}

/// Configuration options for [`summarize_trees`] and [`TreeSummarization`].
#[derive(Debug, Clone)]
pub struct TreeSummaryOptions {
    pub method: SummaryMethod,
    pub target_tree: Option<PathBuf>,
    pub burnin: BurninSpec,
    /// Minimum clade frequency threshold; `None` defaults to `0.5` for MRC, `1.0` for Strict, and `0.0` for MRE.
    pub min_clade_prob: Option<f64>,
    pub rooted_mode: RootedMode,
    pub mcc_scoring: MccScoring,
    pub branch_lengths: BranchLengthStrategy,
    pub support_format: SupportFormat,
    pub support_placement: SupportPlacement,
    pub hpd_level: f64,
    pub summarize_annotations: bool,
    pub allow_partial_taxa: bool,
    pub rooting: RootingMethod,
    pub output_fmt: TreeOutputFmt,
    pub write_opts: TreeWriteOptions,
}

impl TreeSummaryOptions {
    /// Return the resolved minimum clade probability for the configured `method`.
    pub fn effective_min_clade_prob(&self) -> f64 {
        self.min_clade_prob.unwrap_or(match self.method {
            SummaryMethod::Strict => 1.0,
            SummaryMethod::Mrc => 0.5,
            SummaryMethod::Mre
            | SummaryMethod::Mcc
            | SummaryMethod::Map
            | SummaryMethod::Target => 0.0,
        })
    }
}

impl Default for TreeSummaryOptions {
    fn default() -> Self {
        Self {
            method: SummaryMethod::Mcc,
            target_tree: None,
            burnin: BurninSpec::Count(0),
            min_clade_prob: None,
            rooted_mode: RootedMode::Auto,
            mcc_scoring: MccScoring::LogProduct,
            branch_lengths: BranchLengthStrategy::Mean,
            support_format: SupportFormat::Proportion,
            support_placement: SupportPlacement::Auto,
            hpd_level: 0.95,
            summarize_annotations: true,
            allow_partial_taxa: false,
            rooting: RootingMethod::None,
            output_fmt: TreeOutputFmt::Newick,
            write_opts: TreeWriteOptions::default(),
        }
    }
}

/// Summary statistics for a single clade or bipartition, serialized to CSV and SEGUI DTOs.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CladeSummaryEntry {
    pub clade_id: usize,
    pub ntax_in_clade: usize,
    pub count: u64,
    pub decisive_trees: u64,
    pub support: f64,
    pub support_percent: f64,
    pub mean_length: Option<f64>,
    pub median_length: Option<f64>,
    pub length_hpd_low: Option<f64>,
    pub length_hpd_high: Option<f64>,
    pub mean_height: Option<f64>,
    pub median_height: Option<f64>,
    pub height_hpd_low: Option<f64>,
    pub height_hpd_high: Option<f64>,
    pub ca_height: Option<f64>,
    pub taxa: String,
}

/// Summary row for a unique topology in the posterior or bootstrap collection.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TopologySummaryEntry {
    pub rank: usize,
    pub count: u64,
    pub frequency: f64,
    pub cumulative_frequency: f64,
    pub n_clades: usize,
}

/// Complete in-memory result of summarizing a tree collection.
#[derive(Debug, Clone)]
pub struct TreeSummaryResult {
    pub tree: Tree,
    pub stats: TreeSummaryStats,
    pub clades: Vec<CladeSummaryEntry>,
    pub topologies: Vec<TopologySummaryEntry>,
}

/// Files written by [`TreeSummarization`].
#[derive(Debug, Clone)]
pub struct TreeSummaryFiles {
    /// Annotated summary tree file path.
    pub tree_file: PathBuf,
    /// Per-clade summary CSV path.
    pub clade_summary: PathBuf,
    /// Topology credible-set CSV path (omitted in partial-taxa target mode).
    pub topology_summary: Option<PathBuf>,
    /// Unified `report.json` path, if written.
    pub report: Option<PathBuf>,
}

/// Summarize an in-memory slice of [`Tree`]s (SEGUI previews, PySEGUL, and tests).
pub fn summarize_trees(
    trees: &[Tree],
    target: Option<&Tree>,
    options: &TreeSummaryOptions,
) -> Result<TreeSummaryResult> {
    if trees.is_empty() {
        bail!("Cannot summarize an empty tree slice");
    }
    validate_options(options, target.is_some())?;

    let mut data = summarize_in_memory_pass1(trees, target, options)?;
    if !data.rooted && options.branch_lengths == BranchLengthStrategy::CommonAncestor {
        bail!(
            "Common Ancestor (ca) heights require rooted trees (specify --rooted or use mean/median/keep)"
        );
    }

    let total_u64 = data.total_summarized as u64;
    let min_prob = options.effective_min_clade_prob();
    let (mut summary_tree, node_clades, node_taxa, mcc_log_score) = match options.method {
        SummaryMethod::Strict | SummaryMethod::Mrc | SummaryMethod::Mre => {
            let (t, nc) = build_consensus_tree(
                &data.table,
                total_u64,
                options.method,
                min_prob,
                data.rooted,
            )?;
            let (_, nt) = map_tree_nodes_to_clades(&t, &mut data.table, data.rooted, false)?;
            (t, nc, nt, None)
        }
        SummaryMethod::Mcc => {
            let (loc, score) = select_mcc_tree(
                &data.tree_records,
                &data.table,
                &data.topologies,
                total_u64,
                options.mcc_scoring,
            )?;
            let TreeLocation::InMemory(idx) = loc else {
                unreachable!("In-memory summary always returns InMemory location");
            };
            let winner = trees[idx].clone();
            let winner = if options.allow_partial_taxa {
                engine::prune_tree_to_taxa(&winner, &data.table.taxa)?
            } else {
                winner
            };
            let (nc, nt) = map_tree_nodes_to_clades(&winner, &mut data.table, data.rooted, false)?;
            (winner, nc, nt, Some(score))
        }
        SummaryMethod::Map => {
            let (loc, _, score) = select_map_tree(
                &data.topologies,
                &data.table,
                total_u64,
                options.mcc_scoring,
            )?;
            let TreeLocation::InMemory(idx) = loc else {
                unreachable!("In-memory summary always returns InMemory location");
            };
            let winner = trees[idx].clone();
            let winner = if options.allow_partial_taxa {
                engine::prune_tree_to_taxa(&winner, &data.table.taxa)?
            } else {
                winner
            };
            let (nc, nt) = map_tree_nodes_to_clades(&winner, &mut data.table, data.rooted, false)?;
            (winner, nc, nt, Some(score))
        }
        SummaryMethod::Target => {
            let tgt = target
                .context("Target method requires a target tree")?
                .clone();
            let (nc, nt) = map_tree_nodes_to_clades(&tgt, &mut data.table, data.rooted, true)?;
            (tgt, nc, nt, None)
        }
    };

    if options.branch_lengths == BranchLengthStrategy::CommonAncestor && data.allow_ca_on_slice() {
        let store_samples = needs_sample_storage(options);
        let pruned_holder;
        let ca_slice = if options.allow_partial_taxa && options.method != SummaryMethod::Target {
            let retained = &trees[data.total_burned_in..];
            pruned_holder = retained
                .iter()
                .map(|t| engine::prune_tree_to_taxa(t, &data.table.taxa))
                .collect::<Result<Vec<_>>>()?;
            &pruned_holder[..]
        } else {
            &trees[data.total_burned_in..]
        };
        compute_ca_heights_in_memory(
            ca_slice,
            0,
            &summary_tree,
            &node_clades,
            &mut data.table,
            store_samples,
        )?;
    }

    summary_tree.metadata.name = Some(DEFAULT_PREFIX.to_string());
    let (final_tree, clades, topologies) = annotate_and_finalize_tree(
        summary_tree,
        &node_clades,
        &node_taxa,
        &mut data.table,
        &data.topologies,
        total_u64,
        data.rooted,
        data.is_partial_target,
        options,
    )?;

    let stats = build_summary_stats(
        0,
        &data,
        &final_tree,
        &node_clades,
        &topologies,
        mcc_log_score,
        options.method,
        Vec::new(),
    );

    Ok(TreeSummaryResult {
        tree: final_tree,
        stats,
        clades,
        topologies,
    })
}

impl engine::SummaryRunData {
    fn allow_ca_on_slice(&self) -> bool {
        self.rooted
    }
}

/// Batch phylogenetic tree summarization task across one or more tree files.
pub struct TreeSummarization<'a> {
    input_fmt: &'a TreeFmt,
    options: TreeSummaryOptions,
    report_context: ReportContext,
    progress: Option<&'a dyn ProgressReporter>,
}

impl<'a> TreeSummarization<'a> {
    /// Create a new batch tree summary task.
    pub fn new(input_fmt: &'a TreeFmt, options: TreeSummaryOptions) -> Self {
        Self {
            input_fmt,
            options,
            report_context: ReportContext::default(),
            progress: None,
        }
    }

    /// Attach frontend provenance and SHA-256 checksum settings for `report.json`.
    pub fn with_report_context(mut self, context: ReportContext) -> Self {
        self.report_context = context;
        self
    }

    /// Attach a [`ProgressReporter`] for per-file progress updates and cancellation.
    pub fn with_progress(mut self, reporter: &'a dyn ProgressReporter) -> Self {
        self.progress = Some(reporter);
        self
    }

    /// Summarize `files` into `output` without writing `report.json`.
    pub fn summarize(
        &self,
        files: &[PathBuf],
        output: &Path,
        prefix: Option<&str>,
    ) -> Result<TreeSummaryFiles> {
        let (files_out, _) = self.run_inner(files, output, prefix, false)?;
        Ok(files_out)
    }

    /// Summarize `files` into `output` and write a unified Schema v2 `report.json`.
    pub fn summarize_with_report(
        &self,
        files: &[PathBuf],
        output: &Path,
        prefix: Option<&str>,
    ) -> Result<(TreeSummaryFiles, RunReport)> {
        let (files_out, report) = self.run_inner(files, output, prefix, true)?;
        Ok((files_out, report.expect("report requested")))
    }

    fn run_inner(
        &self,
        files: &[PathBuf],
        output: &Path,
        prefix: Option<&str>,
        write_report: bool,
    ) -> Result<(TreeSummaryFiles, Option<RunReport>)> {
        if files.is_empty() {
            bail!("No tree files provided for tree summary");
        }
        validate_options(&self.options, self.options.target_tree.is_some())?;

        let mut sorted = files.to_vec();
        alphanumeric_sort::sort_path_slice(&mut sorted);
        for f in &sorted {
            if !f.is_file() {
                bail!("Input tree file not found: {}", f.display());
            }
        }

        let target_tree = if let Some(ref tgt_path) = self.options.target_tree {
            if !tgt_path.is_file() {
                bail!("Target tree file not found: {}", tgt_path.display());
            }
            Some(parse_tree_file(tgt_path, TreeFmt::Auto)?)
        } else {
            None
        };

        fs::create_dir_all(output)
            .with_context(|| format!("Failed creating output directory {}", output.display()))?;
        let snapshot = DirSnapshot::capture(output);

        let clean_prefix = prefix.filter(|p| !p.is_empty());
        let out_prefix = clean_prefix.unwrap_or(DEFAULT_PREFIX);
        let rooting_str = match &self.options.rooting {
            RootingMethod::None => "none".to_string(),
            RootingMethod::Midpoint => "midpoint".to_string(),
            RootingMethod::Outgroup(_) => "outgroup".to_string(),
        };
        let outgroup_opt = match &self.options.rooting {
            RootingMethod::Outgroup(taxa) => Some(taxa.clone()),
            _ => None,
        };
        let min_prob = self.options.effective_min_clade_prob();
        let params = TreeSummaryParams {
            input_format: canonical::tree_fmt(self.input_fmt).to_string(),
            method: self.options.method.canonical_str().to_string(),
            burnin: self.options.burnin.to_string(),
            min_clade_prob: min_prob,
            rooted_mode: self.options.rooted_mode.canonical_str().to_string(),
            mcc_scoring: self.options.mcc_scoring.canonical_str().to_string(),
            branch_lengths: self.options.branch_lengths.canonical_str().to_string(),
            support_format: self.options.support_format.canonical_str().to_string(),
            support_placement: self.options.support_placement.canonical_str().to_string(),
            hpd_level: self.options.hpd_level,
            summarize_annotations: self.options.summarize_annotations,
            allow_partial_taxa: self.options.allow_partial_taxa,
            rooting: rooting_str,
            outgroup: outgroup_opt,
            output_format: canonical::tree_output_fmt(&self.options.output_fmt).to_string(),
            translate: self.options.write_opts.use_translate,
            strip_annotations: !self.options.write_opts.write_annotations,
            gzip: self.options.write_opts.gzip,
            prefix: clean_prefix.map(str::to_string),
        };

        let report_builder = if write_report {
            let mut b = RunReportBuilder::new(TaskKind::TreeSummary, &self.report_context)
                .primary_inputs(&sorted);
            if let Some(ref tgt_path) = self.options.target_tree {
                b = b.auxiliary_inputs([("target", tgt_path.clone())]);
            }
            Some(b.parameters(&params)?)
        } else {
            None
        };

        let spin = utils::set_spinner();
        spin.set_message(format!(
            "Summarizing phylogenetic trees ({})...",
            summary_method_label(self.options.method)
        ));

        progress::report_start(self.progress, sorted.len());
        progress::check_cancelled(self.progress)?;

        // Compute per-file tree counts and burn-in skips.
        let mut plans = Vec::with_capacity(sorted.len());
        for path in &sorted {
            progress::check_cancelled(self.progress)?;
            let total = count_trees_file(path, *self.input_fmt)?;
            if total == 0 {
                bail!("Tree file contains no trees: {}", path.display());
            }
            let skip = compute_burnin_skip(total, self.options.burnin)?;
            plans.push(FilePlan {
                path: path.clone(),
                total_trees: total,
                burnin_skip: skip,
            });
        }

        let (mut data, per_file) = summarize_files_pass1(
            &plans,
            *self.input_fmt,
            target_tree.as_ref(),
            &self.options,
            self.progress,
        )?;

        if !data.rooted && self.options.branch_lengths == BranchLengthStrategy::CommonAncestor {
            bail!(
                "Common Ancestor (ca) heights require rooted trees (specify --rooted or use mean/median/keep)"
            );
        }

        progress::report_stage(self.progress, "Building summary tree");
        let total_u64 = data.total_summarized as u64;
        let parse_opts = ParseOptions {
            parse_annotations: self.options.summarize_annotations,
            ..Default::default()
        };

        let (mut summary_tree, node_clades, node_taxa, mcc_log_score) = match self.options.method {
            SummaryMethod::Strict | SummaryMethod::Mrc | SummaryMethod::Mre => {
                let (t, nc) = build_consensus_tree(
                    &data.table,
                    total_u64,
                    self.options.method,
                    min_prob,
                    data.rooted,
                )?;
                let (_, nt) = map_tree_nodes_to_clades(&t, &mut data.table, data.rooted, false)?;
                (t, nc, nt, None)
            }
            SummaryMethod::Mcc => {
                let (loc, score) = select_mcc_tree(
                    &data.tree_records,
                    &data.table,
                    &data.topologies,
                    total_u64,
                    self.options.mcc_scoring,
                )?;
                let winner = load_tree_at_location(
                    loc,
                    &plans,
                    *self.input_fmt,
                    &parse_opts,
                    &data.table.taxa,
                    self.options.allow_partial_taxa,
                )?;
                let (nc, nt) =
                    map_tree_nodes_to_clades(&winner, &mut data.table, data.rooted, false)?;
                (winner, nc, nt, Some(score))
            }
            SummaryMethod::Map => {
                let (loc, _, score) = select_map_tree(
                    &data.topologies,
                    &data.table,
                    total_u64,
                    self.options.mcc_scoring,
                )?;
                let winner = load_tree_at_location(
                    loc,
                    &plans,
                    *self.input_fmt,
                    &parse_opts,
                    &data.table.taxa,
                    self.options.allow_partial_taxa,
                )?;
                let (nc, nt) =
                    map_tree_nodes_to_clades(&winner, &mut data.table, data.rooted, false)?;
                (winner, nc, nt, Some(score))
            }
            SummaryMethod::Target => {
                let tgt = target_tree
                    .as_ref()
                    .context("Target method requires a target tree")?
                    .clone();
                let (nc, nt) = map_tree_nodes_to_clades(&tgt, &mut data.table, data.rooted, true)?;
                (tgt, nc, nt, None)
            }
        };

        if self.options.branch_lengths == BranchLengthStrategy::CommonAncestor {
            progress::report_stage(self.progress, "Computing Common Ancestor heights");
            let store_samples = needs_sample_storage(&self.options);
            let prune_to_shared =
                self.options.allow_partial_taxa && self.options.method != SummaryMethod::Target;
            compute_ca_heights_files(
                &plans,
                *self.input_fmt,
                &parse_opts,
                &summary_tree,
                &node_clades,
                &mut data.table,
                store_samples,
                prune_to_shared,
                self.progress,
            )?;
        }

        progress::report_stage(self.progress, "Writing summary outputs");
        summary_tree.metadata.name = Some(out_prefix.to_string());
        let (final_tree, clades, topologies) = annotate_and_finalize_tree(
            summary_tree,
            &node_clades,
            &node_taxa,
            &mut data.table,
            &data.topologies,
            total_u64,
            data.rooted,
            data.is_partial_target,
            &self.options,
        )?;

        let mut tree_filename = format!("{out_prefix}.{}", self.options.output_fmt.extension());
        if self.options.write_opts.gzip {
            tree_filename.push_str(".gz");
        }
        let tree_path = output.join(tree_filename);
        let writer = TreeWriter::new(
            &tree_path,
            self.options.output_fmt,
            self.options.write_opts.clone(),
        );
        writer.write_trees(std::slice::from_ref(&final_tree))?;

        let clade_csv_path = output.join(format!("{out_prefix}_clades.csv"));
        write_clades_csv(&clade_csv_path, &clades)?;

        let topo_csv_path = if !data.is_partial_target {
            let p = output.join(format!("{out_prefix}_topologies.csv"));
            write_topologies_csv(&p, &topologies, MAX_TOPOLOGY_CSV_ROWS)?;
            Some(p)
        } else {
            None
        };

        let stats = build_summary_stats(
            sorted.len(),
            &data,
            &final_tree,
            &node_clades,
            &topologies,
            mcc_log_score,
            self.options.method,
            per_file,
        );

        spin.finish_with_message("Finished summarizing phylogenetic trees!\n");
        self.print_output_info(
            output,
            &tree_path,
            &clade_csv_path,
            topo_csv_path.as_deref(),
            &stats,
        );

        let report_out = if let Some(builder) = report_builder {
            let written = snapshot.written_files(output);
            let rpath = report_path(output, clean_prefix);
            let report = builder
                .outputs(
                    output,
                    written.into_iter().map(|p| {
                        let role = if p == tree_path {
                            "tree"
                        } else if p == clade_csv_path {
                            "clades"
                        } else {
                            "topologies"
                        };
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
            TreeSummaryFiles {
                tree_file: tree_path,
                clade_summary: clade_csv_path,
                topology_summary: topo_csv_path,
                report: report_file,
            },
            report_obj,
        ))
    }

    fn print_output_info(
        &self,
        output: &Path,
        tree_path: &Path,
        clade_csv: &Path,
        topo_csv: Option<&Path>,
        stats: &TreeSummaryStats,
    ) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Output dir", output.display());
        log::info!(
            "{:18}: {}",
            "Method",
            summary_method_label(self.options.method)
        );
        log::info!("{:18}: {}", "Trees read", utils::fmt_num(&stats.trees_read));
        if stats.trees_burned_in > 0 {
            log::info!(
                "{:18}: {}",
                "Burned in",
                utils::fmt_num(&stats.trees_burned_in)
            );
        }
        log::info!(
            "{:18}: {}",
            "Trees summarized",
            utils::fmt_num(&stats.trees_summarized)
        );
        log::info!("{:18}: {}", "Taxa", utils::fmt_num(&stats.ntax));
        if stats.union_ntax > stats.ntax {
            log::info!("{:18}: {}", "Union taxa", utils::fmt_num(&stats.union_ntax));
        }
        if stats.unique_topologies > 0 {
            log::info!(
                "{:18}: {}",
                "Unique topologies",
                utils::fmt_num(&stats.unique_topologies)
            );
        }
        if let Some(cs95) = stats.credible_set_95_topologies {
            log::info!("{:18}: {}", "95% credible set", utils::fmt_num(&cs95));
        }
        log::info!("{:18}: {}", "Summary tree", tree_path.display());
        log::info!("{:18}: {}", "Clade CSV", clade_csv.display());
        if let Some(tp) = topo_csv {
            log::info!("{:18}: {}", "Topologies CSV", tp.display());
        }
        println!();
    }
}

fn load_tree_at_location(
    loc: TreeLocation,
    plans: &[FilePlan],
    input_fmt: TreeFmt,
    parse_opts: &ParseOptions,
    shared_taxa: &crate::parser::tree::TaxonTable,
    prune: bool,
) -> Result<Tree> {
    let TreeLocation::File { file_idx, tree_idx } = loc else {
        unreachable!("File summary always returns File location");
    };
    let tree = read_single_tree_at(&plans[file_idx].path, input_fmt, parse_opts, tree_idx)?;
    if prune {
        engine::prune_tree_to_taxa(&tree, shared_taxa)
    } else {
        Ok(tree)
    }
}

fn validate_options(options: &TreeSummaryOptions, has_target: bool) -> Result<()> {
    if options.method == SummaryMethod::Target && !has_target {
        bail!("Summary method 'target' requires a target tree (--target <FILE>)");
    }
    if options.method != SummaryMethod::Target && options.target_tree.is_some() {
        bail!("--target can only be used with `--method target`");
    }
    if !(0.0..1.0).contains(&options.hpd_level) {
        bail!(
            "HPD level must be in [0.0, 1.0) (0 disables HPD), got {}",
            options.hpd_level
        );
    }
    let min_prob = options.effective_min_clade_prob();
    if !(0.0..=1.0).contains(&min_prob) {
        bail!(
            "Minimum clade probability must be in [0.0, 1.0], got {}",
            min_prob
        );
    }
    if options.method == SummaryMethod::Mrc && min_prob < 0.5 {
        bail!(
            "Majority-rule consensus (--method mrc) requires threshold >= 0.5, got {}",
            min_prob
        );
    }
    if let BurninSpec::Fraction(f) = options.burnin
        && !(0.0..1.0).contains(&f)
    {
        bail!("Burn-in fraction must be in [0.0, 1.0), got {}", f);
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn build_summary_stats(
    files: usize,
    data: &engine::SummaryRunData,
    tree: &Tree,
    node_clades: &[Option<bitset::CladeId>],
    topologies: &[TopologySummaryEntry],
    mcc_log_score: Option<f64>,
    method: SummaryMethod,
    per_file: Vec<TreeSummaryFileStats>,
) -> TreeSummaryStats {
    let (tree_length, tree_height) = tree_length_and_height(tree, data.rooted);
    let summary_internal_nodes = tree
        .preorder()
        .filter(|&id| !tree.is_leaf(id) && !tree.is_root(id))
        .count();

    let total_u64 = data.total_summarized as u64;
    let mut supports: Vec<f64> = node_clades
        .iter()
        .flatten()
        .map(|cid| {
            let st = &data.table.clades[cid.index()];
            let denom = if data.is_partial_target {
                st.decisive_count
            } else {
                total_u64
            };
            if denom > 0 {
                st.count as f64 / denom as f64
            } else {
                0.0
            }
        })
        .collect();
    let mean_support = if supports.is_empty() {
        None
    } else {
        Some(supports.iter().sum::<f64>() / (supports.len() as f64))
    };
    let median_support = if supports.is_empty() {
        None
    } else {
        supports.sort_unstable_by(f64::total_cmp);
        stats_acc::compute_median(&supports)
    };

    let credible_set_95_topologies = if topologies.is_empty() {
        None
    } else {
        let idx = topologies
            .iter()
            .position(|t| t.cumulative_frequency >= 0.95 - 1e-12)
            .map(|i| i + 1)
            .unwrap_or(topologies.len());
        Some(idx)
    };
    let map_topology_frequency = topologies.first().map(|t| t.frequency);

    TreeSummaryStats {
        files,
        trees_read: data.total_read,
        trees_burned_in: data.total_burned_in,
        trees_summarized: data.total_summarized,
        unique_topologies: topologies.len(),
        credible_set_95_topologies,
        map_topology_frequency,
        ntax: data.table.taxa.len(),
        union_ntax: data.union_ntax,
        rooted: data.rooted,
        method: method.canonical_str().to_string(),
        summary_internal_nodes,
        mean_support,
        median_support,
        tree_length,
        tree_height,
        mcc_log_score,
        per_file,
    }
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::{AtomicUsize, Ordering};

    use assert_approx_eq::assert_approx_eq;
    use tempfile::TempDir;

    use crate::helper::progress::CancelAfter;
    use crate::parser::tree::{parse_tree_str, parse_trees_file, parse_trees_str};

    use super::*;

    fn unrooted_bootstrap_trees() -> Vec<Tree> {
        // 10 unrooted 5-taxon trees:
        // {A,B}: 10/10 = 1.0
        // {C,D}: 7/10 = 0.7 (conflicts with {C,E} in 3 trees)
        // Wait: on 5 taxa, non-trivial splits are size 2 (or 3).
        // Let's use 6 taxa so we can test 100%, 70%, 40%, 30%:
        // Trees 1..=7: ((A,B),(C,D),(E,F)) -> {A,B}=10, {C,D}=7, {E,F}=7
        // Trees 8..=10: ((A,B),(C,E),(D,F)) -> {A,B}=10, {C,E}=3, {D,F}=3
        let mut out = Vec::new();
        for _ in 0..7 {
            out.push(
                parse_tree_str("((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2,(E:0.1,F:0.1):0.2);").unwrap(),
            );
        }
        for _ in 0..3 {
            out.push(
                parse_tree_str("((A:0.1,B:0.1):0.4,(C:0.1,E:0.1):0.1,(D:0.1,F:0.1):0.1);").unwrap(),
            );
        }
        out
    }

    #[test]
    fn test_all_methods_analytical_unrooted() {
        let trees = unrooted_bootstrap_trees();

        // Strict: only {A,B} (100%)
        let strict = summarize_trees(
            &trees,
            None,
            &TreeSummaryOptions {
                method: SummaryMethod::Strict,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(strict.stats.summary_internal_nodes, 1);

        // MRC: {A,B} (1.0), {C,D} (0.7), {E,F} (0.7) -> 3 internal nodes
        let mrc = summarize_trees(
            &trees,
            None,
            &TreeSummaryOptions {
                method: SummaryMethod::Mrc,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(mrc.stats.summary_internal_nodes, 3);
        assert_eq!(mrc.stats.unique_topologies, 2);
        assert_approx_eq!(mrc.stats.map_topology_frequency.unwrap(), 0.7, 1e-9);

        // MRE, MCC, MAP all pick the 70% topology with 3 internal nodes
        for method in [SummaryMethod::Mre, SummaryMethod::Mcc, SummaryMethod::Map] {
            let res = summarize_trees(
                &trees,
                None,
                &TreeSummaryOptions {
                    method,
                    ..Default::default()
                },
            )
            .unwrap();
            assert_eq!(res.stats.summary_internal_nodes, 3);
            assert_approx_eq!(
                res.stats.mean_support.unwrap(),
                (1.0 + 0.7 + 0.7) / 3.0,
                1e-6
            );
        }

        // Target with alternative topology ((A,B),(C,E),(D,F)) -> supports 1.0, 0.3, 0.3
        let target = parse_tree_str("((A:1,B:1):1,(C:1,E:1):1,(D:1,F:1):1);").unwrap();
        let tgt_res = summarize_trees(
            &trees,
            Some(&target),
            &TreeSummaryOptions {
                method: SummaryMethod::Target,
                ..Default::default()
            },
        )
        .unwrap();
        assert_approx_eq!(
            tgt_res.stats.mean_support.unwrap(),
            (1.0 + 0.3 + 0.3) / 3.0,
            1e-6
        );
    }

    #[test]
    fn test_ca_heights_guarantee_nonnegative_branches_and_beast_annotations() {
        let input = "\
            [&R] ((A[&rate=0.1,loc=north]:1.0,B[&rate=0.2,loc=north]:1.0)[&rate=0.15]:2.0,C:3.0);\n\
            [&R] ((A[&rate=0.3,loc=north]:1.5,B[&rate=0.4,loc=south]:1.5)[&rate=0.25]:1.5,C:3.0);\n\
            [&R] ((A[&rate=0.2,loc=north]:0.5,C:0.5):2.5,B:3.0);\n";
        let trees = parse_trees_str(input).unwrap();
        let res = summarize_trees(
            &trees,
            None,
            &TreeSummaryOptions {
                method: SummaryMethod::Mcc,
                branch_lengths: BranchLengthStrategy::CommonAncestor,
                hpd_level: 0.95,
                ..Default::default()
            },
        )
        .unwrap();

        for id in res.tree.preorder() {
            if !res.tree.is_root(id) {
                let bl = res.tree.branch_length(id).unwrap();
                assert!(bl >= 0.0, "CA branch length must be non-negative, got {bl}");
            }
        }
        let nwk = crate::writer::tree::to_newick_string(&res.tree, &Default::default());
        assert!(nwk.contains("height_95%_HPD="));
        assert!(nwk.contains("rate="));
        assert!(nwk.contains("loc=north"));
    }

    #[test]
    fn test_partial_taxa_target_gcf_and_shared_pruning() {
        let target = parse_tree_str("((A:1,B:1):1,(C:1,D:1):1);").unwrap();
        let trees = parse_trees_str(
            "\
            ((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\n\
            ((A:0.1,B:0.1):0.2,(C:0.1,E:0.1):0.2);\n\
            ((A:0.1,C:0.1):0.2,(B:0.1,D:0.1):0.2);\n",
        )
        .unwrap();
        // Unrooted split {A,B}|{C,D}:
        // Tree 1: {A,B}|{C,D} -> decisive & supporting
        // Tree 2: has {A,B,C} only (1 taxon on right) -> NOT decisive!
        // Tree 3: {A,C}|{B,D} -> decisive & NOT supporting
        // So {A,B}|{C,D} has decisive=2, supporting=1 -> support = 0.5!
        let res = summarize_trees(
            &trees,
            Some(&target),
            &TreeSummaryOptions {
                method: SummaryMethod::Target,
                rooted_mode: RootedMode::Unrooted,
                allow_partial_taxa: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(res.stats.ntax, 4);
        assert_eq!(res.stats.union_ntax, 5);
        assert_eq!(res.clades.len(), 1);
        assert_eq!(res.clades[0].decisive_trees, 2);
        assert_eq!(res.clades[0].count, 1);
        assert_approx_eq!(res.clades[0].support, 0.5, 1e-9);

        // Non-target Mrc with allow_partial_taxa prunes to shared taxa {A,B,C}
        let mrc = summarize_trees(
            &trees,
            None,
            &TreeSummaryOptions {
                method: SummaryMethod::Mrc,
                rooted_mode: RootedMode::Unrooted,
                allow_partial_taxa: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(mrc.stats.ntax, 3);
        assert_eq!(mrc.stats.union_ntax, 5);
    }

    #[test]
    fn test_progress_and_cancellation_and_format_roundtrip() {
        struct CountProg(AtomicUsize);
        impl ProgressReporter for CountProg {
            fn on_file_finished(&self, _path: &Path) {
                self.0.fetch_add(1, Ordering::Relaxed);
            }
        }

        let dir = TempDir::with_prefix("segul_tree_sum_test").unwrap();
        let f1 = dir.path().join("chain1.nwk");
        let f2 = dir.path().join("chain2.nwk");
        let content = "[&R] ((A:0.1,B:0.1):0.2,C:0.3);\n[&R] ((A:0.2,B:0.2):0.1,C:0.3);\n";
        fs::write(&f1, content).unwrap();
        fs::write(&f2, content).unwrap();

        let out_cancel = dir.path().join("out_cancel");
        let cancel = CancelAfter::new(1);
        let err = TreeSummarization::new(&TreeFmt::Newick, TreeSummaryOptions::default())
            .with_progress(&cancel)
            .summarize(&[f1.clone(), f2.clone()], &out_cancel, None)
            .unwrap_err();
        assert!(
            err.downcast_ref::<crate::helper::progress::Cancelled>()
                .is_some()
        );

        for fmt in [
            TreeOutputFmt::Newick,
            TreeOutputFmt::Nexus,
            TreeOutputFmt::PhyloXml,
            TreeOutputFmt::NeXml,
        ] {
            let out = dir.path().join(format!("out_{fmt}"));
            let prog = CountProg(AtomicUsize::new(0));
            let files = TreeSummarization::new(
                &TreeFmt::Newick,
                TreeSummaryOptions {
                    output_fmt: fmt,
                    ..Default::default()
                },
            )
            .with_progress(&prog)
            .summarize(&[f1.clone(), f2.clone()], &out, Some("run1"))
            .unwrap();
            assert_eq!(prog.0.load(Ordering::Relaxed), 2);
            assert!(files.tree_file.is_file());
            assert!(files.clade_summary.is_file());
            let parsed = parse_trees_file(&files.tree_file, TreeFmt::Auto).unwrap();
            assert_eq!(parsed.len(), 1);
            assert_eq!(parsed[0].tip_count(), 3);
        }
    }

    #[test]
    fn test_high_throughput_determinism_1000_trees_100_taxa() {
        // Build two 100-taxon pectinate/balanced trees and summarize 1,000 trees.
        let mut s1 = String::new();
        let mut s2 = String::new();
        for i in 0..99 {
            s1.push_str(&format!("(T{i}:0.01,"));
            s2.push_str(&format!("(T{i}:0.02,"));
        }
        s1.push_str("T99:0.01");
        s2.push_str("T99:0.02");
        for _ in 0..99 {
            s1.push_str("):0.01");
            s2.push_str("):0.02");
        }
        s1.push(';');
        s2.push(';');

        let t1 = parse_tree_str(&s1).unwrap();
        let t2 = parse_tree_str(&s2).unwrap();
        let mut trees = Vec::with_capacity(1000);
        for i in 0..1000 {
            trees.push(if i % 2 == 0 { t1.clone() } else { t2.clone() });
        }

        let opts = TreeSummaryOptions {
            method: SummaryMethod::Mcc,
            branch_lengths: BranchLengthStrategy::CommonAncestor,
            hpd_level: 0.95,
            ..Default::default()
        };
        let r1 = summarize_trees(&trees, None, &opts).unwrap();
        let r2 = summarize_trees(&trees, None, &opts).unwrap();
        assert_eq!(
            crate::writer::tree::to_newick_string(&r1.tree, &Default::default()),
            crate::writer::tree::to_newick_string(&r2.tree, &Default::default())
        );
        assert_eq!(r1.stats.ntax, 100);
        assert_eq!(r1.stats.trees_summarized, 1000);
    }

    #[test]
    fn test_mre_resolves_below_half_by_default() {
        // 5 rooted 4-taxon trees:
        // 2x ((A,B),(C,D)) -> {A,B}=2 (40%), {C,D}=2 (40%)
        // 1x ((A,C),(B,D)) -> {A,C}=1 (20%), {B,D}=1 (20%)
        // 1x ((A,D),(B,C)) -> {A,D}=1 (20%), {B,C}=1 (20%)
        // 1x (A,B,C,D)
        // MRC has 0 clades > 50%; MRE with default min_clade_prob=None resolves {A,B} and {C,D} (40%).
        let input = "\
            [&R] ((A:1,B:1):1,(C:1,D:1):1);\n\
            [&R] ((A:1,B:1):1,(C:1,D:1):1);\n\
            [&R] ((A:1,C:1):1,(B:1,D:1):1);\n\
            [&R] ((A:1,D:1):1,(B:1,C:1):1);\n\
            [&R] (A:1,B:1,C:1,D:1);\n";
        let trees = parse_trees_str(input).unwrap();
        let mrc = summarize_trees(
            &trees,
            None,
            &TreeSummaryOptions {
                method: SummaryMethod::Mrc,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(mrc.stats.summary_internal_nodes, 0);

        let mre = summarize_trees(
            &trees,
            None,
            &TreeSummaryOptions {
                method: SummaryMethod::Mre,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(mre.stats.summary_internal_nodes, 2);
        assert_approx_eq!(mre.stats.mean_support.unwrap(), 0.4, 1e-9);
    }
}
