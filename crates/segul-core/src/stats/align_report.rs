//! Alignment summary statistics as one JSON report.
//!
//! The report is written to `alignment-report.json` next to the summary
//! CSV files, so applications can draw their own charts without parsing
//! the CSV and text outputs.
use std::collections::BTreeMap;
use std::ffi::OsStr;
use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use crate::helper::types::{DataType, TaxonRecords};
use crate::stats::sequence::{CharMatrix, CharSummary, Completeness, SiteSummary, Sites, Taxa};

/// Every statistic of an alignment summary.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AlignReport {
    /// `dna` or `aa`.
    pub datatype: String,
    pub summary: AlignOverview,
    /// Loci with at least each share of taxa, from 100% down.
    pub completeness: Vec<CompletenessLevel>,
    /// Count of each character across all alignments.
    pub chars: BTreeMap<String, usize>,
    /// One entry per alignment, sorted by path.
    pub loci: Vec<LocusReport>,
    /// One entry per taxon, sorted by name.
    pub taxa: Vec<TaxonReport>,
}

/// Totals across all alignments.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AlignOverview {
    pub taxa: usize,
    pub loci: usize,
    pub sites: usize,
    pub chars: usize,
    pub missing: usize,
    pub prop_missing: f64,
    /// DNA only.
    pub gc: Option<f64>,
    /// DNA only.
    pub at: Option<f64>,
    /// DNA only.
    pub nucleotides: Option<usize>,
    pub min_sites: usize,
    pub max_sites: usize,
    pub mean_sites: f64,
    pub min_taxa: usize,
    pub max_taxa: usize,
    pub mean_taxa: f64,
    pub conserved: SiteTypeSummary,
    pub variable: SiteTypeSummary,
    pub informative: SiteTypeSummary,
}

/// Conserved, variable, or parsimony informative sites across alignments.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SiteTypeSummary {
    /// Loci with at least one site of this type.
    pub loci: usize,
    pub prop_loci: f64,
    pub sites: usize,
    pub prop_sites: f64,
    /// Fewest, most, and mean sites of this type in one locus.
    pub min: usize,
    pub max: usize,
    pub mean: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CompletenessLevel {
    /// Share of taxa, in percent.
    pub percent: usize,
    pub loci: usize,
}

/// Statistics of one alignment.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LocusReport {
    pub path: String,
    pub name: String,
    pub ntax: usize,
    pub chars: usize,
    pub sites: usize,
    pub conserved: usize,
    pub variable: usize,
    pub informative: usize,
    pub missing: usize,
    pub prop_missing: f64,
    /// DNA only.
    pub gc: Option<f64>,
    /// DNA only.
    pub at: Option<f64>,
}

/// Statistics of one taxon across alignments.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TaxonReport {
    pub name: String,
    /// Loci the taxon appears in.
    pub loci: usize,
    pub chars: usize,
    pub missing: usize,
    pub prop_missing: f64,
    /// DNA only.
    pub gc: Option<f64>,
    /// DNA only.
    pub at: Option<f64>,
}

impl AlignReport {
    /// Build the report from the statistics of an alignment summary.
    pub fn new(
        sites: &SiteSummary,
        chars: &CharSummary,
        complete: &Completeness,
        stats: &[(Sites, CharMatrix, Taxa)],
        taxa: &BTreeMap<String, TaxonRecords>,
        datatype: &DataType,
    ) -> Self {
        let is_dna = DataType::Dna == *datatype;
        let dna = |value: f64| if is_dna { Some(value) } else { None };
        let summary = AlignOverview {
            taxa: complete.total_tax,
            loci: sites.total_loci,
            sites: sites.total_sites,
            chars: chars.total_chars,
            missing: chars.missing_data,
            prop_missing: chars.prop_missing_data,
            gc: dna(chars.gc_content),
            at: dna(chars.at_content),
            nucleotides: is_dna.then_some(chars.total_nucleotides),
            min_sites: sites.min_sites,
            max_sites: sites.max_sites,
            mean_sites: sites.mean_sites,
            min_taxa: chars.min_tax,
            max_taxa: chars.max_tax,
            mean_taxa: chars.mean_tax,
            conserved: SiteTypeSummary {
                loci: sites.cons_loci,
                prop_loci: sites.prop_cons_loci,
                sites: sites.total_cons_site,
                prop_sites: sites.prop_cons_site,
                min: sites.min_cons_site,
                max: sites.max_cons_site,
                mean: sites.mean_cons_site,
            },
            variable: SiteTypeSummary {
                loci: sites.var_loci,
                prop_loci: sites.prop_var_loci,
                sites: sites.total_var_site,
                prop_sites: sites.prop_var_site,
                min: sites.min_var_site,
                max: sites.max_var_site,
                mean: sites.mean_var_site,
            },
            informative: SiteTypeSummary {
                loci: sites.inf_loci,
                prop_loci: sites.prop_inf_loci,
                sites: sites.total_inf_site,
                prop_sites: sites.prop_inf_site,
                min: sites.min_inf_site,
                max: sites.max_inf_site,
                mean: sites.mean_inf_site,
            },
        };
        let loci = stats
            .iter()
            .map(|(site, matrix, _)| {
                let total = matrix.chars.total_chars;
                LocusReport {
                    path: site.path.display().to_string(),
                    name: site
                        .path
                        .file_stem()
                        .and_then(OsStr::to_str)
                        .unwrap_or_default()
                        .to_string(),
                    ntax: matrix.ntax,
                    chars: total,
                    sites: site.counts,
                    conserved: site.conserved,
                    variable: site.variable,
                    informative: site.pars_inf,
                    missing: matrix.chars.missing_data,
                    prop_missing: matrix.chars.prop_missing_data,
                    gc: dna(proportion(matrix.chars.gc_count, total)),
                    at: dna(proportion(matrix.chars.at_count, total)),
                }
            })
            .collect();
        let taxa = taxa
            .iter()
            .map(|(name, taxon)| TaxonReport {
                name: name.to_string(),
                loci: taxon.locus_counts,
                chars: taxon.total_chars,
                missing: taxon.missing_data,
                prop_missing: proportion(taxon.missing_data, taxon.total_chars),
                gc: dna(proportion(taxon.gc_count, taxon.total_chars)),
                at: dna(proportion(taxon.at_count, taxon.total_chars)),
            })
            .collect();
        Self {
            datatype: if is_dna { "dna" } else { "aa" }.to_string(),
            summary,
            completeness: complete
                .completeness
                .iter()
                .map(|(percent, loci)| CompletenessLevel {
                    percent: *percent,
                    loci: *loci,
                })
                .collect(),
            chars: chars
                .chars
                .iter()
                .map(|(ch, count)| (ch.to_string(), *count))
                .collect(),
            loci,
            taxa,
        }
    }

    /// Read a report written by `segul align summary`.
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

fn proportion(count: usize, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        count as f64 / total as f64
    }
}
