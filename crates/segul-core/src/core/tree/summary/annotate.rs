//! Summary tree annotation, branch length / height assignment, and CSV output writers.

use std::fmt::Write as _;
use std::path::Path;

use anyhow::{Context, Result};
use serde::Serialize;

use crate::core::tree::estimate::{RootingMethod, reroot_midpoint, reroot_outgroup};
use crate::parser::tree::{
    Annotation, AnnotationFormat, AnnotationValue, NodeId, Tree, TreeOutputFmt,
};

use super::bitset::{CladeId, popcount};
use super::engine::{CladeStats, GlobalCladeTable, TopologyRecord};
use super::stats_acc::{AnnotationKeyTable, AnnotationValueAcc, MetricSummary};
use super::{
    BranchLengthStrategy, CladeSummaryEntry, SummaryMethod, SupportFormat, SupportPlacement,
    TopologySummaryEntry, TreeSummaryOptions,
};

#[derive(Debug, Serialize)]
struct CladeCsvRow<'a> {
    clade_id: usize,
    ntax_in_clade: usize,
    count: u64,
    decisive_trees: u64,
    support: f64,
    support_percent: f64,
    mean_length: Option<f64>,
    median_length: Option<f64>,
    length_hpd_low: Option<f64>,
    length_hpd_high: Option<f64>,
    mean_height: Option<f64>,
    median_height: Option<f64>,
    height_hpd_low: Option<f64>,
    height_hpd_high: Option<f64>,
    ca_height: Option<f64>,
    taxa: &'a str,
}

#[derive(Debug, Serialize)]
struct TopologyCsvRow {
    rank: usize,
    count: u64,
    frequency: f64,
    cumulative_frequency: f64,
    n_clades: usize,
}

/// Apply branch lengths, heights, support values, and BEAST/FigTree annotations to `tree`.
///
/// Also returns finalized clade and topology summary entries for reports and CSV files.
#[allow(clippy::too_many_arguments)]
pub(crate) fn annotate_and_finalize_tree(
    mut tree: Tree,
    node_clades: &[Option<CladeId>],
    node_taxa: &[Option<usize>],
    table: &mut GlobalCladeTable,
    topologies: &[TopologyRecord],
    total_trees: u64,
    rooted: bool,
    is_partial_target: bool,
    options: &TreeSummaryOptions,
) -> Result<(Tree, Vec<CladeSummaryEntry>, Vec<TopologySummaryEntry>)> {
    let hpd_level = options.hpd_level;
    let write_annots = options.summarize_annotations && options.write_opts.write_annotations;

    // Strip existing annotations on winning/target tree so synthesized summary tags are clean.
    tree.clear_annotations();

    // Finalize all clade metrics once so CSV and tree annotations share exact values.
    let mut clade_summaries = Vec::with_capacity(table.clades.len());
    let mut clade_entries = Vec::with_capacity(table.clades.len());
    for (idx, stats) in table.clades.iter_mut().enumerate() {
        let len_sum = stats.branch_lengths.finalize(hpd_level);
        let h_sum = stats.heights.finalize(hpd_level);
        let ca_sum = stats.ca_heights.finalize(hpd_level);

        let denom = if is_partial_target {
            stats.decisive_count
        } else {
            total_trees
        };
        let support = if denom > 0 {
            stats.count as f64 / denom as f64
        } else {
            0.0
        };
        let bits = table.arena.get(CladeId(idx as u32));
        let ntax_in_clade = popcount(bits);
        let taxa_label = format_clade_taxa(bits, &table.taxa);

        clade_entries.push(CladeSummaryEntry {
            clade_id: idx + 1,
            ntax_in_clade,
            count: stats.count,
            decisive_trees: denom,
            support,
            support_percent: support * 100.0,
            mean_length: len_sum.mean,
            median_length: len_sum.median,
            length_hpd_low: len_sum.hpd_low,
            length_hpd_high: len_sum.hpd_high,
            mean_height: h_sum.mean,
            median_height: h_sum.median,
            height_hpd_low: h_sum.hpd_low,
            height_hpd_high: h_sum.hpd_high,
            ca_height: ca_sum.mean,
            taxa: taxa_label,
        });
        clade_summaries.push((len_sum, h_sum, ca_sum, support, denom));
    }

    // Finalize tip and root metric summaries.
    let mut tip_summaries = Vec::with_capacity(table.tip_stats.len());
    for tip in &mut table.tip_stats {
        tip_summaries.push((
            tip.branch_lengths.finalize(hpd_level),
            tip.heights.finalize(hpd_level),
            tip.ca_heights.finalize(hpd_level),
        ));
    }
    let root_len_sum = table.root_stats.branch_lengths.finalize(hpd_level);
    let root_h_sum = table.root_stats.heights.finalize(hpd_level);
    let root_ca_sum = table.root_stats.ca_heights.finalize(hpd_level);

    // Precompute chosen node heights when CommonAncestor is active so branch lengths are >= 0.
    let mut chosen_heights = vec![0.0f64; tree.len()];
    if rooted && options.branch_lengths == BranchLengthStrategy::CommonAncestor {
        for id in tree.preorder() {
            let idx = id.index();
            let h = if tree.is_root(id) {
                root_ca_sum.mean.or(root_h_sum.mean).unwrap_or(0.0)
            } else if tree.is_leaf(id) {
                if let Some(tid) = node_taxa[idx] {
                    tip_summaries[tid]
                        .2
                        .mean
                        .or(tip_summaries[tid].1.mean)
                        .unwrap_or(0.0)
                } else {
                    0.0
                }
            } else if let Some(cid) = node_clades[idx] {
                let (_, ref h_s, ref ca_s, _, _) = clade_summaries[cid.index()];
                ca_s.mean.or(h_s.mean).unwrap_or(0.0)
            } else {
                0.0
            };
            chosen_heights[idx] = h;
        }
    }

    let placement =
        resolve_support_placement(options.support_placement, options.output_fmt, write_annots);

    let keys = &table.keys;
    let node_ids: Vec<NodeId> = tree.preorder().collect();

    for id in node_ids {
        let idx = id.index();
        if tree.is_root(id) {
            if options.branch_lengths == BranchLengthStrategy::None {
                tree.clear_branch_length(id);
            }
            if write_annots
                && let Some(annot) = build_beast_annotation(
                    None,
                    None,
                    &root_len_sum,
                    &root_h_sum,
                    &root_ca_sum,
                    &mut table.root_stats,
                    keys,
                    rooted,
                    false,
                    options,
                )
            {
                tree.set_node_annotations(id, vec![annot]);
            }
            continue;
        }

        if tree.is_leaf(id) {
            let Some(tid) = node_taxa[idx] else {
                continue;
            };
            let (ref len_s, ref h_s, ref ca_s) = tip_summaries[tid];
            apply_branch_length(
                &mut tree,
                id,
                options.branch_lengths,
                len_s,
                &chosen_heights,
            );
            if write_annots
                && let Some(annot) = build_beast_annotation(
                    None,
                    None,
                    len_s,
                    h_s,
                    ca_s,
                    &mut table.tip_stats[tid],
                    keys,
                    rooted,
                    false,
                    options,
                )
            {
                tree.set_node_annotations(id, vec![annot]);
            }
            continue;
        }

        // Internal non-root node.
        let Some(cid) = node_clades[idx] else {
            // Unmapped internal node (e.g., second child of degree-2 root in unrooted mode)
            if options.branch_lengths == BranchLengthStrategy::None {
                tree.clear_branch_length(id);
            }
            continue;
        };

        let (ref len_s, ref h_s, ref ca_s, support, denom) = clade_summaries[cid.index()];
        apply_branch_length(
            &mut tree,
            id,
            options.branch_lengths,
            len_s,
            &chosen_heights,
        );

        let emit_label = matches!(
            placement,
            ResolvedPlacement::Label | ResolvedPlacement::Both
        );
        let emit_annot = matches!(
            placement,
            ResolvedPlacement::Annotation | ResolvedPlacement::Both
        );

        if emit_label {
            let formatted = format_support_value(support, options.support_format);
            tree.set_label(id, &formatted, false);
        } else {
            tree.clear_label(id);
        }

        if write_annots {
            let sup_opt = if emit_annot { Some(support) } else { None };
            let partial_counts = if is_partial_target {
                Some((table.clades[cid.index()].count, denom))
            } else {
                None
            };
            if let Some(annot) = build_beast_annotation(
                sup_opt,
                partial_counts,
                len_s,
                h_s,
                ca_s,
                &mut table.clades[cid.index()],
                keys,
                rooted,
                true,
                options,
            ) {
                tree.set_node_annotations(id, vec![annot]);
            }
        }
    }

    // Optional post-summary rooting (`--outgroup` or `--midpoint`).
    let final_tree = match &options.rooting {
        RootingMethod::None => tree,
        RootingMethod::Midpoint => reroot_midpoint(&tree)?,
        RootingMethod::Outgroup(outgroup) => reroot_outgroup(&tree, outgroup)?,
    };

    // Sort clades descending by support, then size, then taxa label.
    clade_entries.sort_by(|a, b| {
        b.support
            .total_cmp(&a.support)
            .then_with(|| b.ntax_in_clade.cmp(&a.ntax_in_clade))
            .then_with(|| a.taxa.cmp(&b.taxa))
    });
    for (rank, entry) in clade_entries.iter_mut().enumerate() {
        entry.clade_id = rank + 1;
    }

    let topo_entries = build_topology_entries(topologies, total_trees);

    Ok((final_tree, clade_entries, topo_entries))
}

#[derive(Clone, Copy)]
enum ResolvedPlacement {
    Annotation,
    Label,
    Both,
}

fn resolve_support_placement(
    placement: SupportPlacement,
    output_fmt: TreeOutputFmt,
    write_annots: bool,
) -> ResolvedPlacement {
    match placement {
        SupportPlacement::Annotation => {
            if write_annots {
                ResolvedPlacement::Annotation
            } else {
                ResolvedPlacement::Label
            }
        }
        SupportPlacement::Label => ResolvedPlacement::Label,
        SupportPlacement::Both => {
            if write_annots {
                ResolvedPlacement::Both
            } else {
                ResolvedPlacement::Label
            }
        }
        SupportPlacement::Auto => {
            if !write_annots {
                ResolvedPlacement::Label
            } else if matches!(output_fmt, TreeOutputFmt::Nexus | TreeOutputFmt::PhyloXml) {
                ResolvedPlacement::Annotation
            } else {
                ResolvedPlacement::Both
            }
        }
    }
}

fn apply_branch_length(
    tree: &mut Tree,
    id: NodeId,
    strategy: BranchLengthStrategy,
    len_summary: &MetricSummary,
    chosen_heights: &[f64],
) {
    match strategy {
        BranchLengthStrategy::Keep => {}
        BranchLengthStrategy::None => {
            tree.clear_branch_length(id);
        }
        BranchLengthStrategy::Mean => {
            if let Some(m) = len_summary.mean {
                tree.set_branch_length(id, clean_zero(m));
            }
        }
        BranchLengthStrategy::Median => {
            if let Some(m) = len_summary.median.or(len_summary.mean) {
                tree.set_branch_length(id, clean_zero(m));
            }
        }
        BranchLengthStrategy::CommonAncestor => {
            if let Some(parent) = tree.parent(id) {
                let bl = (chosen_heights[parent.index()] - chosen_heights[id.index()]).max(0.0);
                tree.set_branch_length(id, clean_zero(bl));
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn build_beast_annotation(
    support: Option<f64>,
    partial_counts: Option<(u64, u64)>,
    len_sum: &MetricSummary,
    h_sum: &MetricSummary,
    ca_sum: &MetricSummary,
    stats: &mut CladeStats,
    keys: &AnnotationKeyTable,
    rooted: bool,
    is_internal: bool,
    options: &TreeSummaryOptions,
) -> Option<Annotation> {
    let mut entries: Vec<(String, AnnotationValue)> = Vec::new();
    let hpd_pct = format_hpd_pct(options.hpd_level);

    if let Some(sup) = support {
        let key = if rooted { "posterior" } else { "support" };
        let val = match options.support_format {
            SupportFormat::Proportion => clean_zero(round_sig(sup, 6)),
            SupportFormat::Percent => clean_zero(round_sig(sup * 100.0, 4)),
        };
        entries.push((key.to_string(), AnnotationValue::Float(val)));
    }
    if let Some((supporting, decisive)) = partial_counts {
        entries.push((
            "supporting".to_string(),
            AnnotationValue::Int(supporting as i64),
        ));
        entries.push((
            "decisive".to_string(),
            AnnotationValue::Int(decisive as i64),
        ));
    }

    if is_internal || len_sum.count > 0 {
        push_metric_tags(
            &mut entries,
            "length",
            len_sum,
            &hpd_pct,
            options.branch_lengths,
            false,
        );
    }

    if rooted && (h_sum.count > 0 || ca_sum.count > 0) {
        let primary_h = if options.branch_lengths == BranchLengthStrategy::CommonAncestor {
            ca_sum.mean.or(h_sum.mean)
        } else {
            h_sum.mean
        };
        let active_h_sum =
            if options.branch_lengths == BranchLengthStrategy::CommonAncestor && ca_sum.count > 0 {
                ca_sum
            } else {
                h_sum
            };
        if let Some(h) = primary_h {
            entries.push((
                "height".to_string(),
                AnnotationValue::Float(clean_zero(round_sig(h, 6))),
            ));
        }
        if let Some(med) = active_h_sum.median {
            entries.push((
                "height_median".to_string(),
                AnnotationValue::Float(clean_zero(round_sig(med, 6))),
            ));
        }
        if let (Some(lo), Some(hi)) = (active_h_sum.hpd_low, active_h_sum.hpd_high) {
            entries.push((
                format!("height_{hpd_pct}%_HPD"),
                AnnotationValue::Array(vec![
                    AnnotationValue::Float(clean_zero(round_sig(lo, 6))),
                    AnnotationValue::Float(clean_zero(round_sig(hi, 6))),
                ]),
            ));
        }
        if let (Some(min), Some(max)) = (active_h_sum.min, active_h_sum.max) {
            entries.push((
                "height_range".to_string(),
                AnnotationValue::Array(vec![
                    AnnotationValue::Float(clean_zero(round_sig(min, 6))),
                    AnnotationValue::Float(clean_zero(round_sig(max, 6))),
                ]),
            ));
        }
    }

    // Generic BEAST/FigTree/NHX numeric & categorical tags.
    if !stats.annotations.is_empty() {
        let mut collected: Vec<(String, &mut AnnotationValueAcc)> = stats
            .annotations
            .iter_mut()
            .filter_map(|(kid, acc)| keys.get(kid).map(|name| (name.to_string(), acc)))
            .collect();
        collected.sort_by(|a, b| a.0.cmp(&b.0));

        for (name, acc) in collected {
            match acc {
                AnnotationValueAcc::Numeric(m_acc) => {
                    let sum = m_acc.finalize(options.hpd_level);
                    if let Some(mean) = sum.mean {
                        entries.push((
                            name.clone(),
                            AnnotationValue::Float(clean_zero(round_sig(mean, 6))),
                        ));
                    }
                    if let Some(med) = sum.median {
                        entries.push((
                            format!("{name}_median"),
                            AnnotationValue::Float(clean_zero(round_sig(med, 6))),
                        ));
                    }
                    if let (Some(lo), Some(hi)) = (sum.hpd_low, sum.hpd_high) {
                        entries.push((
                            format!("{name}_{hpd_pct}%_HPD"),
                            AnnotationValue::Array(vec![
                                AnnotationValue::Float(clean_zero(round_sig(lo, 6))),
                                AnnotationValue::Float(clean_zero(round_sig(hi, 6))),
                            ]),
                        ));
                    }
                }
                AnnotationValueAcc::Categorical(counts) => {
                    let total: u64 = counts.values().sum();
                    if total > 0
                        && let Some((best_val, best_cnt)) = counts
                            .iter()
                            .max_by(|(ka, ca), (kb, cb)| ca.cmp(cb).then_with(|| kb.cmp(ka)))
                    {
                        entries.push((name.clone(), AnnotationValue::Str(best_val.clone())));
                        entries.push((
                            format!("{name}.prob"),
                            AnnotationValue::Float(clean_zero(round_sig(
                                *best_cnt as f64 / total as f64,
                                4,
                            ))),
                        ));
                    }
                }
            }
        }
    }

    if entries.is_empty() {
        return None;
    }

    let raw = format_beast_raw(&entries);
    Some(Annotation {
        format: AnnotationFormat::Beast,
        entries,
        raw,
    })
}

fn push_metric_tags(
    entries: &mut Vec<(String, AnnotationValue)>,
    prefix: &str,
    sum: &MetricSummary,
    hpd_pct: &str,
    _strategy: BranchLengthStrategy,
    include_range: bool,
) {
    if let Some(mean) = sum.mean {
        entries.push((
            prefix.to_string(),
            AnnotationValue::Float(clean_zero(round_sig(mean, 6))),
        ));
    }
    if let Some(med) = sum.median {
        entries.push((
            format!("{prefix}_median"),
            AnnotationValue::Float(clean_zero(round_sig(med, 6))),
        ));
    }
    if let (Some(lo), Some(hi)) = (sum.hpd_low, sum.hpd_high) {
        entries.push((
            format!("{prefix}_{hpd_pct}%_HPD"),
            AnnotationValue::Array(vec![
                AnnotationValue::Float(clean_zero(round_sig(lo, 6))),
                AnnotationValue::Float(clean_zero(round_sig(hi, 6))),
            ]),
        ));
    }
    if include_range && let (Some(min), Some(max)) = (sum.min, sum.max) {
        entries.push((
            format!("{prefix}_range"),
            AnnotationValue::Array(vec![
                AnnotationValue::Float(clean_zero(round_sig(min, 6))),
                AnnotationValue::Float(clean_zero(round_sig(max, 6))),
            ]),
        ));
    }
}

fn format_beast_raw(entries: &[(String, AnnotationValue)]) -> String {
    let mut out = String::from("&");
    for (i, (k, v)) in entries.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(k);
        out.push('=');
        append_annot_value(&mut out, v);
    }
    out
}

fn append_annot_value(out: &mut String, val: &AnnotationValue) {
    match val {
        AnnotationValue::Str(s) => {
            if s.contains([',', '=', ' ', '{', '}', '[', ']']) {
                out.push('"');
                out.push_str(&s.replace('"', "'"));
                out.push('"');
            } else {
                out.push_str(s);
            }
        }
        AnnotationValue::Float(v) => {
            let _ = write!(out, "{}", format_float(*v));
        }
        AnnotationValue::Int(v) => {
            let _ = write!(out, "{v}");
        }
        AnnotationValue::Bool(v) => {
            let _ = write!(out, "{v}");
        }
        AnnotationValue::Array(items) => {
            out.push('{');
            for (i, item) in items.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                append_annot_value(out, item);
            }
            out.push('}');
        }
    }
}

fn format_support_value(support: f64, format: SupportFormat) -> String {
    match format {
        SupportFormat::Proportion => format_float(clean_zero(round_sig(support, 4))),
        SupportFormat::Percent => format_float(clean_zero(round_sig(support * 100.0, 2))),
    }
}

fn format_float(v: f64) -> String {
    let v = clean_zero(v);
    if v == 0.0 {
        return "0".to_string();
    }
    let abs = v.abs();
    if !(1e-5..1e7).contains(&abs) {
        return format!("{v:.6e}");
    }
    let s = format!("{v:.6}");
    let trimmed = s.trim_end_matches('0').trim_end_matches('.');
    trimmed.to_string()
}

fn round_sig(v: f64, decimals: i32) -> f64 {
    if !v.is_finite() {
        return v;
    }
    let factor = 10f64.powi(decimals);
    (v * factor).round() / factor
}

#[inline]
fn clean_zero(v: f64) -> f64 {
    if v == 0.0 { 0.0 } else { v }
}

fn format_hpd_pct(level: f64) -> String {
    let pct = level * 100.0;
    if (pct - pct.round()).abs() < 1e-6 {
        format!("{}", pct.round() as u32)
    } else {
        format_float(round_sig(pct, 2))
    }
}

fn format_clade_taxa(bits: &[u64], taxa: &crate::parser::tree::TaxonTable) -> String {
    let mut names = Vec::new();
    for (w_idx, &word) in bits.iter().enumerate() {
        let mut w = word;
        let base = w_idx * 64;
        while w != 0 {
            let tz = w.trailing_zeros() as usize;
            let tid = base + tz;
            if let Some(name) = taxa.get(crate::parser::tree::TaxonId::from_raw(tid as u32)) {
                names.push(name);
            }
            w &= w - 1;
        }
    }
    names.join(",")
}

fn build_topology_entries(
    topologies: &[TopologyRecord],
    total_trees: u64,
) -> Vec<TopologySummaryEntry> {
    if topologies.is_empty() || total_trees == 0 {
        return Vec::new();
    }
    let mut sorted: Vec<&TopologyRecord> = topologies.iter().collect();
    sorted.sort_by(|a, b| {
        b.count
            .cmp(&a.count)
            .then_with(|| a.first_location.cmp(&b.first_location))
    });
    let total_f = total_trees as f64;
    let mut cum = 0.0f64;
    sorted
        .into_iter()
        .enumerate()
        .map(|(i, rec)| {
            let freq = rec.count as f64 / total_f;
            cum += freq;
            TopologySummaryEntry {
                rank: i + 1,
                count: rec.count,
                frequency: freq,
                cumulative_frequency: cum.min(1.0),
                n_clades: rec.clade_ids.len(),
            }
        })
        .collect()
}

pub(crate) fn write_clades_csv(path: &Path, clades: &[CladeSummaryEntry]) -> Result<()> {
    let mut w = csv::Writer::from_path(path)
        .with_context(|| format!("Failed creating {}", path.display()))?;
    for c in clades {
        w.serialize(CladeCsvRow {
            clade_id: c.clade_id,
            ntax_in_clade: c.ntax_in_clade,
            count: c.count,
            decisive_trees: c.decisive_trees,
            support: c.support,
            support_percent: c.support_percent,
            mean_length: c.mean_length,
            median_length: c.median_length,
            length_hpd_low: c.length_hpd_low,
            length_hpd_high: c.length_hpd_high,
            mean_height: c.mean_height,
            median_height: c.median_height,
            height_hpd_low: c.height_hpd_low,
            height_hpd_high: c.height_hpd_high,
            ca_height: c.ca_height,
            taxa: &c.taxa,
        })?;
    }
    w.flush()?;
    Ok(())
}

pub(crate) fn write_topologies_csv(
    path: &Path,
    topologies: &[TopologySummaryEntry],
    max_rows: usize,
) -> Result<()> {
    let mut w = csv::Writer::from_path(path)
        .with_context(|| format!("Failed creating {}", path.display()))?;
    for t in topologies.iter().take(max_rows) {
        w.serialize(TopologyCsvRow {
            rank: t.rank,
            count: t.count,
            frequency: t.frequency,
            cumulative_frequency: t.cumulative_frequency,
            n_clades: t.n_clades,
        })?;
    }
    w.flush()?;
    Ok(())
}

pub(crate) fn tree_length_and_height(tree: &Tree, rooted: bool) -> (Option<f64>, Option<f64>) {
    let mut has_any = false;
    let mut sum = 0.0f64;
    let mut depths = vec![0.0f64; tree.len()];
    let mut max_depth = 0.0f64;

    for id in tree.preorder() {
        let d = if let Some(parent) = tree.parent(id) {
            let bl = tree.branch_length(id);
            if let Some(len) = bl
                && len.is_finite()
            {
                has_any = true;
                sum += len;
                depths[parent.index()] + len.max(0.0)
            } else {
                depths[parent.index()]
            }
        } else {
            0.0
        };
        depths[id.index()] = d;
        if tree.is_leaf(id) && d > max_depth {
            max_depth = d;
        }
    }

    let len = if has_any { Some(sum) } else { None };
    let height = if rooted && has_any {
        Some(max_depth)
    } else {
        None
    };
    (len, height)
}

pub(crate) fn summary_method_label(method: SummaryMethod) -> &'static str {
    match method {
        SummaryMethod::Mcc => "Maximum Clade Credibility (MCC)",
        SummaryMethod::Mrc => "Majority-Rule Consensus (MRC)",
        SummaryMethod::Mre => "Extended Majority-Rule Consensus (MRE)",
        SummaryMethod::Strict => "Strict Consensus",
        SummaryMethod::Map => "Maximum A Posteriori (MAP)",
        SummaryMethod::Target => "Target Tree Annotation",
    }
}
