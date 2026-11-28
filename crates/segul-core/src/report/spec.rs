//! Canonical task kinds and strongly typed parameter specifications shared
//! across SEGUL CLI, SEGUI, and PySEGUL.

use std::str::FromStr;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use serde_json::Value;

/// Canonical identifier for every SEGUL operation that emits `report.json`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum TaskKind {
    #[serde(rename = "align.summary")]
    AlignSummary,
    #[serde(rename = "align.distance")]
    AlignDistance,
    #[serde(rename = "align.concat")]
    AlignConcat,
    #[serde(rename = "align.convert")]
    AlignConvert,
    #[serde(rename = "align.filter")]
    AlignFilter,
    #[serde(rename = "align.split")]
    AlignSplit,
    #[serde(rename = "align.trim")]
    AlignTrim,
    #[serde(rename = "align.unalign")]
    AlignUnalign,
    #[serde(rename = "partition.convert")]
    PartitionConvert,
    #[serde(rename = "contig.summary")]
    ContigSummary,
    #[serde(rename = "read.summary")]
    ReadSummary,
    #[serde(rename = "maf.summary")]
    MafSummary,
    #[serde(rename = "maf.convert")]
    MafConvert,
    #[serde(rename = "sequence.add")]
    SequenceAdd,
    #[serde(rename = "sequence.extract")]
    SequenceExtract,
    #[serde(rename = "sequence.filter")]
    SequenceFilter,
    #[serde(rename = "sequence.id")]
    SequenceId,
    #[serde(rename = "sequence.remove")]
    SequenceRemove,
    #[serde(rename = "sequence.rename")]
    SequenceRename,
    #[serde(rename = "sequence.translate")]
    SequenceTranslate,
    #[serde(rename = "sequence.transpose")]
    SequenceTranspose,
    #[serde(rename = "tree.convert")]
    TreeConvert,
}

impl TaskKind {
    /// Dotted canonical string name (e.g., `"align.summary"`).
    pub fn as_str(self) -> &'static str {
        match self {
            Self::AlignSummary => "align.summary",
            Self::AlignDistance => "align.distance",
            Self::AlignConcat => "align.concat",
            Self::AlignConvert => "align.convert",
            Self::AlignFilter => "align.filter",
            Self::AlignSplit => "align.split",
            Self::AlignTrim => "align.trim",
            Self::AlignUnalign => "align.unalign",
            Self::PartitionConvert => "partition.convert",
            Self::ContigSummary => "contig.summary",
            Self::ReadSummary => "read.summary",
            Self::MafSummary => "maf.summary",
            Self::MafConvert => "maf.convert",
            Self::SequenceAdd => "sequence.add",
            Self::SequenceExtract => "sequence.extract",
            Self::SequenceFilter => "sequence.filter",
            Self::SequenceId => "sequence.id",
            Self::SequenceRemove => "sequence.remove",
            Self::SequenceRename => "sequence.rename",
            Self::SequenceTranslate => "sequence.translate",
            Self::SequenceTranspose => "sequence.transpose",
            Self::TreeConvert => "tree.convert",
        }
    }

    /// Concise human-readable title for UI and Methods text.
    pub fn label(self) -> &'static str {
        match self {
            Self::AlignSummary => "Alignment summary",
            Self::AlignDistance => "Pairwise genetic distance",
            Self::AlignConcat => "Alignment concatenation",
            Self::AlignConvert => "Alignment format conversion",
            Self::AlignFilter => "Alignment filtering",
            Self::AlignSplit => "Alignment splitting",
            Self::AlignTrim => "Alignment trimming",
            Self::AlignUnalign => "Alignment unalignment",
            Self::PartitionConvert => "Partition conversion",
            Self::ContigSummary => "Contig summary",
            Self::ReadSummary => "Read summary",
            Self::MafSummary => "MAF summary",
            Self::MafConvert => "MAF conversion",
            Self::SequenceAdd => "Sequence addition",
            Self::SequenceExtract => "Sequence extraction",
            Self::SequenceFilter => "Sequence filtering",
            Self::SequenceId => "Sequence ID extraction",
            Self::SequenceRemove => "Sequence removal",
            Self::SequenceRename => "Sequence ID renaming",
            Self::SequenceTranslate => "Sequence translation",
            Self::SequenceTranspose => "Sequence transposition",
            Self::TreeConvert => "Tree format conversion",
        }
    }

    /// Whether this task has an interactive chart/table viewer in CLI/SEGUI.
    pub fn has_interactive_viewer(self) -> bool {
        matches!(
            self,
            Self::AlignSummary | Self::AlignDistance | Self::ContigSummary | Self::ReadSummary
        )
    }
}

impl std::fmt::Display for TaskKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FromStr for TaskKind {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.trim() {
            "align.summary" => Ok(Self::AlignSummary),
            "align.distance" => Ok(Self::AlignDistance),
            "align.concat" => Ok(Self::AlignConcat),
            "align.convert" => Ok(Self::AlignConvert),
            "align.filter" => Ok(Self::AlignFilter),
            "align.split" => Ok(Self::AlignSplit),
            "align.trim" => Ok(Self::AlignTrim),
            "align.unalign" => Ok(Self::AlignUnalign),
            "partition.convert" => Ok(Self::PartitionConvert),
            "contig.summary" => Ok(Self::ContigSummary),
            "read.summary" => Ok(Self::ReadSummary),
            "maf.summary" => Ok(Self::MafSummary),
            "maf.convert" => Ok(Self::MafConvert),
            "sequence.add" => Ok(Self::SequenceAdd),
            "sequence.extract" => Ok(Self::SequenceExtract),
            "sequence.filter" => Ok(Self::SequenceFilter),
            "sequence.id" => Ok(Self::SequenceId),
            "sequence.remove" => Ok(Self::SequenceRemove),
            "sequence.rename" => Ok(Self::SequenceRename),
            "sequence.translate" => Ok(Self::SequenceTranslate),
            "sequence.transpose" => Ok(Self::SequenceTranspose),
            "tree.convert" => Ok(Self::TreeConvert),
            _ => Err(format!("Unknown SEGUL task kind '{}'", s)),
        }
    }
}

/// Typed parameters for all 21 SEGUL tasks, dispatched explicitly by [`TaskKind`].
#[derive(Debug, Clone, PartialEq)]
pub enum TaskParameters {
    AlignSummary(AlignSummaryParams),
    AlignDistance(AlignDistanceParams),
    AlignConcat(AlignConcatParams),
    AlignConvert(AlignConvertParams),
    AlignFilter(AlignFilterParams),
    AlignSplit(AlignSplitParams),
    AlignTrim(AlignTrimParams),
    AlignUnalign(AlignUnalignParams),
    PartitionConvert(PartitionConvertParams),
    ContigSummary(ContigSummaryParams),
    ReadSummary(ReadSummaryParams),
    MafSummary(MafSummaryParams),
    MafConvert(MafConvertParams),
    SequenceAdd(SequenceAddParams),
    SequenceExtract(SequenceExtractParams),
    SequenceFilter(SequenceFilterParams),
    SequenceId(SequenceIdParams),
    SequenceRemove(SequenceRemoveParams),
    SequenceRename(SequenceRenameParams),
    SequenceTranslate(SequenceTranslateParams),
    SequenceTranspose(SequenceTransposeParams),
    TreeConvert(TreeConvertParams),
}

impl TaskParameters {
    /// Parse raw JSON parameters according to `kind`.
    pub fn from_json_value(kind: TaskKind, value: &Value) -> Result<Self> {
        let v = value.clone();
        let parsed = match kind {
            TaskKind::AlignSummary => Self::AlignSummary(serde_json::from_value(v)?),
            TaskKind::AlignDistance => Self::AlignDistance(serde_json::from_value(v)?),
            TaskKind::AlignConcat => Self::AlignConcat(serde_json::from_value(v)?),
            TaskKind::AlignConvert => Self::AlignConvert(serde_json::from_value(v)?),
            TaskKind::AlignFilter => Self::AlignFilter(serde_json::from_value(v)?),
            TaskKind::AlignSplit => Self::AlignSplit(serde_json::from_value(v)?),
            TaskKind::AlignTrim => Self::AlignTrim(serde_json::from_value(v)?),
            TaskKind::AlignUnalign => Self::AlignUnalign(serde_json::from_value(v)?),
            TaskKind::PartitionConvert => Self::PartitionConvert(serde_json::from_value(v)?),
            TaskKind::ContigSummary => Self::ContigSummary(serde_json::from_value(v)?),
            TaskKind::ReadSummary => Self::ReadSummary(serde_json::from_value(v)?),
            TaskKind::MafSummary => Self::MafSummary(serde_json::from_value(v)?),
            TaskKind::MafConvert => Self::MafConvert(serde_json::from_value(v)?),
            TaskKind::SequenceAdd => Self::SequenceAdd(serde_json::from_value(v)?),
            TaskKind::SequenceExtract => Self::SequenceExtract(serde_json::from_value(v)?),
            TaskKind::SequenceFilter => Self::SequenceFilter(serde_json::from_value(v)?),
            TaskKind::SequenceId => Self::SequenceId(serde_json::from_value(v)?),
            TaskKind::SequenceRemove => Self::SequenceRemove(serde_json::from_value(v)?),
            TaskKind::SequenceRename => Self::SequenceRename(serde_json::from_value(v)?),
            TaskKind::SequenceTranslate => Self::SequenceTranslate(serde_json::from_value(v)?),
            TaskKind::SequenceTranspose => Self::SequenceTranspose(serde_json::from_value(v)?),
            TaskKind::TreeConvert => Self::TreeConvert(serde_json::from_value(v)?),
        };
        Ok(parsed)
    }

    /// Serialize parameters to a JSON object.
    pub fn to_json_value(&self) -> Result<Value> {
        let v = match self {
            Self::AlignSummary(p) => serde_json::to_value(p),
            Self::AlignDistance(p) => serde_json::to_value(p),
            Self::AlignConcat(p) => serde_json::to_value(p),
            Self::AlignConvert(p) => serde_json::to_value(p),
            Self::AlignFilter(p) => serde_json::to_value(p),
            Self::AlignSplit(p) => serde_json::to_value(p),
            Self::AlignTrim(p) => serde_json::to_value(p),
            Self::AlignUnalign(p) => serde_json::to_value(p),
            Self::PartitionConvert(p) => serde_json::to_value(p),
            Self::ContigSummary(p) => serde_json::to_value(p),
            Self::ReadSummary(p) => serde_json::to_value(p),
            Self::MafSummary(p) => serde_json::to_value(p),
            Self::MafConvert(p) => serde_json::to_value(p),
            Self::SequenceAdd(p) => serde_json::to_value(p),
            Self::SequenceExtract(p) => serde_json::to_value(p),
            Self::SequenceFilter(p) => serde_json::to_value(p),
            Self::SequenceId(p) => serde_json::to_value(p),
            Self::SequenceRemove(p) => serde_json::to_value(p),
            Self::SequenceRename(p) => serde_json::to_value(p),
            Self::SequenceTranslate(p) => serde_json::to_value(p),
            Self::SequenceTranspose(p) => serde_json::to_value(p),
            Self::TreeConvert(p) => serde_json::to_value(p),
        };
        v.context("Failed serializing task parameters")
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlignSummaryParams {
    pub input_format: String,
    pub datatype: String,
    pub interval: usize,
    pub per_locus: bool,
    pub prefix: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "method", rename_all = "snake_case")]
pub enum DistanceGroupingSpec {
    Delimiter { delimiter: String, fields: usize },
    Regex { pattern: String },
    MapFile { path: String },
    Auto,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlignDistanceParams {
    pub input_format: String,
    pub datatype: String,
    pub model: String,
    pub triangle: String,
    pub diagonal: Option<String>,
    pub long_format: bool,
    pub grouping: Option<DistanceGroupingSpec>,
    pub resolved_grouping: Option<DistanceGroupingSpec>,
    pub prefix: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlignConcatParams {
    pub input_format: String,
    pub datatype: String,
    pub output_format: String,
    pub partition_format: String,
    #[serde(default)]
    pub sort: bool,
    pub prefix: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlignConvertParams {
    pub input_format: String,
    pub datatype: String,
    pub output_format: String,
    pub sort: bool,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum AlignFilterCriterion {
    MinTaxPercent {
        percent: f64,
        resolved_min_taxa: Option<usize>,
    },
    MinTaxCount {
        min_taxa: usize,
    },
    MinLen {
        min_len: usize,
    },
    MaxLen {
        max_len: usize,
    },
    MinParsInf {
        min_pars_inf: usize,
    },
    MaxParsInf {
        max_pars_inf: usize,
    },
    PercentInf {
        percent_inf: f64,
        resolved_min_pinf: Option<usize>,
    },
    MissingData {
        max_missing: f64,
    },
    TaxonAll {
        taxa: Vec<String>,
    },
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AlignFilterParams {
    pub input_format: String,
    pub datatype: String,
    pub criterion: AlignFilterCriterion,
    pub concat: bool,
    pub output_format: Option<String>,
    pub partition_format: Option<String>,
    pub prefix: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlignSplitParams {
    pub input_format: String,
    pub datatype: String,
    pub output_format: String,
    pub partition_format: String,
    pub skip_checking: bool,
    pub prefix: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum AlignTrimCriterion {
    MissingData { max_missing: f64 },
    ParsInf { min_pars_inf: usize },
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AlignTrimParams {
    pub input_format: String,
    pub datatype: String,
    pub output_format: String,
    pub criterion: AlignTrimCriterion,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlignUnalignParams {
    pub input_format: String,
    pub datatype: String,
    pub output_format: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PartitionConvertParams {
    pub datatype: String,
    pub input_partition_format: String,
    pub output_partition_format: String,
    pub skip_checking: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ContigSummaryParams {
    pub input_format: String,
    pub prefix: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReadSummaryParams {
    pub input_format: String,
    pub mode: String,
    pub phred: String,
    pub dup_sketch_size: usize,
    pub qc_profile: String,
    pub prefix: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct MafSummaryParams {
    pub naming: String,
    pub ref_species: Option<String>,
    pub resolved_ref_species: Option<String>,
    pub species: Option<Vec<String>>,
    pub prefix: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MafConvertParams {
    pub output_format: String,
    pub group_by: String,
    pub clip: bool,
    pub strand_aware: bool,
    pub naming: String,
    pub ref_species: Option<String>,
    pub resolved_ref_species: Option<String>,
    pub species: Option<Vec<String>>,
    pub min_block_len: Option<usize>,
    pub min_species: Option<usize>,
    pub min_score: Option<f64>,
    pub write_unmatched: bool,
    pub prefix: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SequenceAddParams {
    pub input_format: String,
    pub dest_format: String,
    pub datatype: String,
    pub output_format: String,
    pub added_only: bool,
    pub single_file: bool,
    pub include_filename: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SequenceExtractParams {
    pub input_format: String,
    pub datatype: String,
    pub output_format: String,
    pub mode: String,
    pub ids: Option<Vec<String>>,
    pub regex: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum SequenceFilterCriterion {
    MinLength { min_length: usize },
    MaxLength { max_length: usize },
    MaxGap { max_gap: f64 },
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SequenceFilterParams {
    pub input_format: String,
    pub datatype: String,
    pub output_format: String,
    pub criterion: SequenceFilterCriterion,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SequenceIdParams {
    pub input_format: String,
    pub datatype: String,
    pub map_id: bool,
    pub prefix: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SequenceRemoveParams {
    pub input_format: String,
    pub datatype: String,
    pub output_format: String,
    pub mode: String,
    pub ids: Option<Vec<String>>,
    pub regex: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SequenceRenameParams {
    pub input_format: String,
    pub datatype: String,
    pub output_format: String,
    pub dry_run: bool,
    pub mode: String,
    pub from: Option<String>,
    pub to: Option<String>,
    pub all_matches: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SequenceTranslateParams {
    pub input_format: String,
    pub datatype: String,
    pub output_format: String,
    pub table: String,
    pub reading_frame: usize,
    pub auto_frame: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SequenceTransposeParams {
    pub input_format: String,
    pub datatype: String,
    pub output_format: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TreeConvertParams {
    pub input_format: String,
    pub output_format: String,
    pub translate: bool,
    pub strip_annotations: bool,
    pub gzip: bool,
}

/// Canonical string helpers for `segul-core` enums so `report.json`
/// parameters always match `FromStr` across CLI, SEGUI, and PySEGUL.
pub mod canonical {
    use crate::helper::types::{
        ContigFmt, DataType, GeneticCodes, InputFmt, OutputFmt, PartitionFmt, PhredOffset,
        SeqReadFmt, SummaryMode,
    };

    pub fn input_fmt(fmt: &InputFmt) -> &'static str {
        match fmt {
            InputFmt::Auto => "auto",
            InputFmt::Fasta => "fasta",
            InputFmt::Nexus => "nexus",
            InputFmt::Phylip => "phylip",
        }
    }

    pub fn output_fmt(fmt: &OutputFmt) -> &'static str {
        match fmt {
            OutputFmt::Fasta => "fasta",
            OutputFmt::Nexus => "nexus",
            OutputFmt::Phylip => "phylip",
            OutputFmt::FastaInt => "fasta-int",
            OutputFmt::NexusInt => "nexus-int",
            OutputFmt::PhylipInt => "phylip-int",
        }
    }

    pub fn partition_fmt(fmt: &PartitionFmt) -> &'static str {
        match fmt {
            PartitionFmt::Charset => "charset",
            PartitionFmt::CharsetCodon => "charset-codon",
            PartitionFmt::Nexus => "nexus",
            PartitionFmt::NexusCodon => "nexus-codon",
            PartitionFmt::Raxml => "raxml",
            PartitionFmt::RaxmlCodon => "raxml-codon",
        }
    }

    pub fn datatype(dt: &DataType) -> &'static str {
        match dt {
            DataType::Dna => "dna",
            DataType::Aa => "aa",
            DataType::Ignore => "ignore",
        }
    }

    pub fn contig_fmt(fmt: &ContigFmt) -> &'static str {
        match fmt {
            ContigFmt::Auto => "auto",
            ContigFmt::Fasta => "fasta",
            ContigFmt::Gzip => "gzip",
        }
    }

    pub fn read_fmt(fmt: &SeqReadFmt) -> &'static str {
        match fmt {
            SeqReadFmt::Auto => "auto",
            SeqReadFmt::Fastq => "fastq",
            SeqReadFmt::Gzip => "gzip",
        }
    }

    pub fn phred(phred: PhredOffset) -> &'static str {
        match phred {
            PhredOffset::Auto => "auto",
            PhredOffset::Phred33 => "33",
            PhredOffset::Phred64 => "64",
        }
    }

    pub fn summary_mode(mode: &SummaryMode) -> &'static str {
        match mode {
            SummaryMode::Minimal => "minimal",
            SummaryMode::Basic => "basic",
            SummaryMode::Complete => "complete",
        }
    }

    pub fn genetic_code(code: &GeneticCodes) -> &'static str {
        match code {
            GeneticCodes::StandardCode => "1",
            GeneticCodes::VertMtDna => "2",
            GeneticCodes::YeastMtDna => "3",
            GeneticCodes::MoldProtCoelMtDna => "4",
            GeneticCodes::InvertMtDna => "5",
            GeneticCodes::CilDasHexNu => "6",
            GeneticCodes::EchiFlatwormMtDna => "9",
            GeneticCodes::EuplotidNu => "10",
            GeneticCodes::BacArchPlantPlast => "11",
            GeneticCodes::AltYeastNu => "12",
            GeneticCodes::AsciMtDna => "13",
            GeneticCodes::AltFlatwormMtDna => "14",
            GeneticCodes::ChlorMtDna => "16",
            GeneticCodes::TrematodeMtDna => "21",
            GeneticCodes::ScenedesmusMtDna => "22",
            GeneticCodes::ThrausMtDna => "23",
            GeneticCodes::RhabdopMtDna => "24",
            GeneticCodes::CaDivSR1GraciBac => "25",
            GeneticCodes::PachyNu => "26",
            GeneticCodes::MesodiniumNu => "29",
            GeneticCodes::PeritrichNu => "30",
            GeneticCodes::CephalodiscidaeMtDna => "33",
        }
    }

    #[cfg(feature = "align")]
    pub fn distance_model(model: &crate::stats::distance::DistanceModel) -> &'static str {
        use crate::stats::distance::DistanceModel;
        match model {
            DistanceModel::Uncorrected => "p",
            DistanceModel::JukesCantor => "jc69",
            DistanceModel::Kimura2P => "k2p",
            DistanceModel::TamuraNei => "tn93",
            DistanceModel::Poisson => "poisson",
        }
    }

    #[cfg(feature = "align")]
    pub fn matrix_triangle(t: crate::core::align::distance::MatrixTriangle) -> &'static str {
        use crate::core::align::distance::MatrixTriangle;
        match t {
            MatrixTriangle::Full => "full",
            MatrixTriangle::Lower => "lower",
            MatrixTriangle::Upper => "upper",
        }
    }

    #[cfg(feature = "align")]
    pub fn matrix_diagonal(d: crate::core::align::distance::MatrixDiagonal) -> &'static str {
        use crate::core::align::distance::MatrixDiagonal;
        match d {
            MatrixDiagonal::Zero => "zero",
            MatrixDiagonal::Blank => "blank",
            MatrixDiagonal::Na => "na",
            MatrixDiagonal::GroupMean => "mean",
            MatrixDiagonal::GroupMax => "max",
        }
    }

    #[cfg(feature = "maf")]
    pub fn sample_naming(n: crate::core::maf::SampleNaming) -> &'static str {
        use crate::core::maf::SampleNaming;
        match n {
            SampleNaming::AssemblyPrefix => "assembly_prefix",
            SampleNaming::FullSource => "full_source",
        }
    }

    #[cfg(feature = "tree")]
    pub fn tree_fmt(fmt: &crate::parser::tree::TreeFmt) -> &'static str {
        use crate::parser::tree::TreeFmt;
        match fmt {
            TreeFmt::Auto => "auto",
            TreeFmt::Newick => "newick",
            TreeFmt::ExtendedNewick => "enewick",
            TreeFmt::Nexus => "nexus",
        }
    }

    #[cfg(feature = "tree")]
    pub fn tree_output_fmt(fmt: &crate::parser::tree::TreeOutputFmt) -> &'static str {
        use crate::parser::tree::TreeOutputFmt;
        match fmt {
            TreeOutputFmt::Newick => "newick",
            TreeOutputFmt::ExtendedNewick => "enewick",
            TreeOutputFmt::Nexus => "nexus",
        }
    }
}
