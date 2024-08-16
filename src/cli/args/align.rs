use std::path::PathBuf;

use clap::builder::TypedValueParser as _;
use clap::Subcommand;
use clap::{builder, Args};

use super::{CommonConcatArgs, CommonSeqInput, CommonSeqOutput, IOArgs};

#[derive(Subcommand)]
pub(crate) enum AlignmentSubcommand {
    #[command(about = "Concatenate alignments", name = "concat")]
    Concat(AlignConcatArgs),
    #[command(about = "Convert sequence formats", name = "convert")]
    Convert(AlignConvertArgs),
    #[command(about = "Filter alignments", name = "filter")]
    Filter(AlignFilterArgs),
    #[command(about = "Split alignment by partitions", name = "split")]
    Split(AlignSplitArgs),
    #[command(about = "Compute Alignment Statistics", name = "summary")]
    Summary(AlignSummaryArgs),
    #[command(about = "Convert alignment to unaligned sequences", name = "unalign")]
    Unalign(UnalignArgs),
}

#[derive(Subcommand)]
pub(crate) enum PartitionSubcommand {
    #[command(about = "Convert partition formats", name = "convert")]
    Convert(PartitionArgs),
}

#[derive(Args)]
pub(crate) struct AlignConcatArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[command(flatten)]
    pub(crate) out_fmt: CommonSeqOutput,
    #[command(flatten)]
    pub(crate) concat: CommonConcatArgs,
    #[arg(short, long, help = "Output path", default_value = "Align-Concat")]
    pub(crate) output: PathBuf,
    #[arg(long = "sort", help = "Sort sequences by IDs alphabetically")]
    pub(crate) sort: bool,
}

#[derive(Args)]
pub(crate) struct AlignConvertArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[command(flatten)]
    pub(crate) out_fmt: CommonSeqOutput,
    #[arg(short, long, help = "Output path", default_value = "Align-Convert")]
    pub(crate) output: PathBuf,
    #[arg(long = "sort", help = "Sort sequences by IDs alphabetically")]
    pub(crate) sort: bool,
}

#[derive(Args)]
pub(crate) struct AlignFilterArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[command(flatten)]
    pub(crate) out_fmt: CommonSeqOutput,
    #[command(flatten)]
    pub(crate) partition: CommonConcatArgs,
    #[arg(short, long, help = "Output path", default_value = "Align-Filter")]
    pub(crate) output: String,
    #[arg(long = "concat", help = "Concat filtered alignments")]
    pub(crate) concat: bool,
    #[arg(long = "len", help = "Filter by sequence length")]
    pub(crate) len: Option<usize>,
    #[arg(long = "missing-data", help = "Filter by proportion of missing data")]
    pub(crate) missing: Option<f64>,
    #[arg(
        long = "npercent",
        help = "Filter by multiple minimal taxon percentage"
    )]
    pub(crate) npercent: Option<Vec<f64>>,
    #[arg(long = "nseq", help = "Filter by multiple minimal taxon number")]
    pub(crate) ntax: Option<usize>,
    #[arg(long = "percent", help = "Filter by minimal taxon percentage")]
    pub(crate) percent: Option<f64>,
    #[arg(
        long = "percent-inf",
        help = "Filter by minimal parsimony informative percentage"
    )]
    pub(crate) percent_inf: Option<f64>,
    #[arg(long = "pinf", help = "Filter by minimal parsimony informative sites")]
    pub(crate) pinf: Option<usize>,
    #[arg(long = "taxon-id", help = "Filter by taxon ID")]
    pub(crate) ids: Option<PathBuf>,
}

#[derive(Args)]
pub(crate) struct AlignSplitArgs {
    #[arg(short = 'i', long = "input", help = "Input partition path")]
    pub(crate) input: PathBuf,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[command(flatten)]
    pub(crate) out_fmt: CommonSeqOutput,
    #[arg(short, long, help = "Output path", default_value = "Align-Split")]
    pub(crate) output: PathBuf,
    #[arg(short = 'I', long = "input-partition", help = "Input partition format")]
    pub(crate) input_partition: Option<PathBuf>,
    #[arg(long = "skip-checking", help = "Skip checking partition format")]
    pub(crate) skip_checking: bool,
    #[arg(long, help = "Force overwriting existing output files/directory")]
    pub(crate) force: bool,
    #[arg(long = "prefix", help = "Specify prefix for output files")]
    pub(crate) prefix: Option<String>,
    #[arg(
        short = 'p',
        long = "partition-format",
        help = "Specify partition format",
        value_parser = builder::PossibleValuesParser::new(["nexus", "raxml"]),
    )]
    pub(crate) part_fmt: Option<String>,
}

#[derive(Args)]
pub(crate) struct AlignSummaryArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) fmt: CommonSeqInput,
    #[arg(short, long, help = "Output path", default_value = "Align-Summary")]
    pub(crate) output: PathBuf,
    #[arg(long = "prefix", help = "Specify prefix for output files")]
    pub(crate) prefix: Option<String>,
    #[arg(
        long = "interval",
        help = "Specify interval value for counting data matrix completeness",
        default_value_t = 5,
        value_parser = builder::PossibleValuesParser::new(["1", "2", "5", "10"])
            .map(|x| x.parse::<usize>().unwrap_or(5)),
    )]
    pub(crate) interval: usize,
    #[arg(long = "per-locus", help = "Generate summary statistic for each locus")]
    pub(crate) per_locus: bool,
}

#[derive(Args)]
pub(crate) struct PartitionArgs {
    #[arg(short, long, help = "Input a path (allow wildcard)")]
    #[cfg(target_os = "windows")]
    pub(crate) input: Option<String>,
    #[arg(short, long, help = "Input a path (allow wildcard)")]
    #[cfg(not(target_os = "windows"))]
    pub(crate) input: Option<Vec<PathBuf>>,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[command(flatten)]
    pub(crate) out_fmt: CommonSeqOutput,
    #[arg(short, long, help = "Output path", default_value = "Partition-Convert")]
    pub(crate) output: String,
    #[arg(
        short = 'p',
        long = "input-partition",
        help = "Specify partition format",
        value_parser = builder::PossibleValuesParser::new(["charset", "nexus", "raxml"]),
    )]
    pub(crate) part_fmt: Option<String>,
    #[arg(long = "codon", help = "Set codon model partition format")]
    pub(crate) codon: bool,
    #[arg(
        short = 'P',
        long = "output-partition",
        help = "Specify partition format",
        default_value = "nexus",
        value_parser = builder::PossibleValuesParser::new(["charset", "nexus", "raxml"]),
    )]
    pub(crate) out_part: String,
    #[arg(long, help = "Force overwriting existing output files/directory")]
    pub(crate) force: bool,
    #[arg(long = "skip-checking", help = "Skip checking partition formats")]
    pub(crate) skip_checking: bool,
}

#[derive(Args)]
pub(crate) struct UnalignArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[arg(
        short = 'F',
        long = "output-format",
        help = "Specify output format",
        default_value = "fasta-int",
        value_parser = builder::PossibleValuesParser::new(
            ["fasta", "fasta-int"]),
    )]
    pub(crate) output_fmt: String,
    #[arg(
        short,
        long,
        help = "Output directory path",
        default_value = "Align-Unalign"
    )]
    pub(crate) output: PathBuf,
}
