//! Structured statistics payloads embedded inside `RunReport::statistics`.

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use serde_json::Value;

#[cfg(feature = "align")]
use crate::stats::align_report::AlignReport;
#[cfg(feature = "contig")]
use crate::stats::contig_report::ContigReport;
#[cfg(feature = "align")]
use crate::stats::distance_report::DistanceReport;
#[cfg(feature = "read")]
use crate::stats::fastq::report::ReadReport;

use super::spec::TaskKind;

/// Typed `statistics` payload for any SEGUL task, dispatched by [`TaskKind`].
#[derive(Debug, Clone, PartialEq)]
pub enum ReportStatistics {
    #[cfg(feature = "align")]
    Align(AlignReport),
    #[cfg(feature = "align")]
    Distance(DistanceReport),
    #[cfg(feature = "contig")]
    Contig(ContigReport),
    #[cfg(feature = "read")]
    Read(ReadReport),
    MafSummary(MafSummaryReport),
    MafConvert(MafConvertStats),
    Concat(ConcatStats),
    Filter(FilterStats),
    Trim(TrimStats),
    Counts(TaskCounts),
    Raw(Value),
}

impl ReportStatistics {
    /// Parse `statistics` JSON using the report's `kind`.
    pub fn from_json_value(kind: TaskKind, value: &Value) -> Result<Self> {
        let v = value.clone();
        let parsed = match kind {
            #[cfg(feature = "align")]
            TaskKind::AlignSummary => Self::Align(serde_json::from_value(v)?),
            #[cfg(feature = "align")]
            TaskKind::AlignDistance => Self::Distance(serde_json::from_value(v)?),
            #[cfg(feature = "contig")]
            TaskKind::ContigSummary => Self::Contig(serde_json::from_value(v)?),
            #[cfg(feature = "read")]
            TaskKind::ReadSummary => Self::Read(serde_json::from_value(v)?),
            TaskKind::MafSummary => Self::MafSummary(serde_json::from_value(v)?),
            TaskKind::MafConvert => Self::MafConvert(serde_json::from_value(v)?),
            TaskKind::AlignConcat => Self::Concat(serde_json::from_value(v)?),
            TaskKind::AlignFilter => Self::Filter(serde_json::from_value(v)?),
            TaskKind::AlignTrim => Self::Trim(serde_json::from_value(v)?),
            _ => match serde_json::from_value::<TaskCounts>(v.clone()) {
                Ok(counts) => Self::Counts(counts),
                Err(_) => Self::Raw(v),
            },
        };
        Ok(parsed)
    }

    /// Serialize `statistics` to a JSON value.
    pub fn to_json_value(&self) -> Result<Value> {
        let res = match self {
            #[cfg(feature = "align")]
            Self::Align(r) => serde_json::to_value(r),
            #[cfg(feature = "align")]
            Self::Distance(r) => serde_json::to_value(r),
            #[cfg(feature = "contig")]
            Self::Contig(r) => serde_json::to_value(r),
            #[cfg(feature = "read")]
            Self::Read(r) => serde_json::to_value(r),
            Self::MafSummary(r) => serde_json::to_value(r),
            Self::MafConvert(r) => serde_json::to_value(r),
            Self::Concat(r) => serde_json::to_value(r),
            Self::Filter(r) => serde_json::to_value(r),
            Self::Trim(r) => serde_json::to_value(r),
            Self::Counts(r) => serde_json::to_value(r),
            Self::Raw(v) => Ok(v.clone()),
        };
        res.context("Failed serializing report statistics")
    }
}

/// Concatenation summary statistics.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ConcatStats {
    pub taxa: usize,
    pub loci: usize,
    pub sites: usize,
    pub partitions: usize,
}

/// Alignment filtering statistics, including resolved auto-thresholds.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FilterStats {
    pub input_loci: usize,
    pub retained_loci: usize,
    pub dropped_loci: usize,
    pub resolved_min_taxa: Option<usize>,
    pub resolved_min_pinf: Option<usize>,
}

/// Trimming result for a single alignment file.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TrimEntry {
    pub path: String,
    pub initial_sites: usize,
    pub final_sites: usize,
    pub removed_sites: usize,
}

/// Trimming summary across all processed alignments.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TrimStats {
    pub files: usize,
    pub entries: Vec<TrimEntry>,
}

/// MAF conversion counts (mirrors `MafConvertReport`).
#[derive(Debug, Default, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct MafConvertStats {
    pub files: usize,
    pub blocks_read: usize,
    pub blocks_no_ref: usize,
    pub blocks_filtered: usize,
    pub blocks_unmatched: usize,
    pub duplicate_rows: usize,
    pub overlapping_pieces: usize,
    pub loci_written: usize,
    pub loci_empty: usize,
}

/// MAF summary statistics across files and species.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MafSummaryReport {
    pub ref_species: String,
    pub files: Vec<MafFileSummaryReport>,
    pub species: Vec<MafSpeciesSummaryReport>,
}

/// Per-MAF-file summary statistics.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MafFileSummaryReport {
    pub path: String,
    pub blocks: usize,
    pub blocks_no_ref: usize,
    pub species: usize,
    pub ref_chromosomes: usize,
    pub ref_chrom_size_total: usize,
    pub ref_bases: usize,
    pub ref_genome_coverage: Option<f64>,
    pub aligned_columns: usize,
    pub mean_block_columns: Option<f64>,
    pub min_block_columns: usize,
    pub max_block_columns: usize,
    pub mean_block_ref_len: Option<f64>,
    pub mean_species_per_block: Option<f64>,
    pub mean_score: Option<f64>,
    pub gap_proportion: Option<f64>,
    pub duplicate_rows: usize,
}

/// Per-species MAF summary statistics (per file or `"all"` combined).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MafSpeciesSummaryReport {
    pub path: String,
    pub name: String,
    pub blocks_present: usize,
    pub block_proportion: Option<f64>,
    pub aligned_columns: usize,
    pub bases: usize,
    pub gaps: usize,
    pub gap_proportion: Option<f64>,
    pub n_count: usize,
    pub ref_bases_covered: usize,
    pub ref_coverage: Option<f64>,
}

/// Concise counts for transform and utility tasks.
#[derive(Debug, Default, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TaskCounts {
    pub files_read: usize,
    pub files_written: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub items_processed: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub items_skipped: Option<usize>,
}
