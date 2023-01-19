use std::path::{Path, PathBuf};

use clap::{builder, Args, Parser, Subcommand};

#[derive(Parser)]
#[command(author, about, version, long_about = None)]
#[command(propagate_version = true)]
pub(crate) struct Cli {
    #[command(subcommand)]
    pub(crate) subcommand: MainSubcommand,
}

#[derive(Subcommand)]
pub(crate) enum MainSubcommand {
    #[command(subcommand, about = "Raw read sequence analyses", name = "raw")]
    RawRead(RawReadSubcommand),
    #[command(subcommand, about = "Contigous sequence analyses", name = "contig")]
    Contig(ContigSubcommand),
    #[command(subcommand, about = "Alignment analyses", name = "align")]
    Alignment(AlignmentSubcommand),
    #[command(
        subcommand,
        about = "Alignment partition conversion",
        name = "partition"
    )]
    Partition(PartitionSubcommand),
    #[command(subcommand, about = "Sequence analyses", name = "sequence")]
    Sequence(SequenceSubcommand),
}

#[derive(Subcommand)]
pub(crate) enum RawReadSubcommand {
    #[command(about = "Compute raw read statistics", name = "stats")]
    RawStats(RawStatArgs),
}

#[derive(Subcommand)]
pub(crate) enum ContigSubcommand {
    #[command(about = "Compute contig statistics", name = "stats")]
    ContigStats(ContigStatArgs),
}

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
    #[command(about = "Compute Alignment Statistics", name = "stats")]
    AlignStats(AlignSummaryArgs),
}

#[derive(Subcommand)]
pub(crate) enum PartitionSubcommand {
    #[command(about = "Convert partition formats", name = "convert")]
    Convert(PartitionArgs),
}

#[derive(Subcommand)]
pub(crate) enum SequenceSubcommand {
    #[command(about = "Extract sequence from alignments", name = "extract")]
    Extract(SequenceExtractArgs),
    #[command(about = "Parse sample ID across multiple alignments", name = "id")]
    Id(SequenceIdArgs),
    #[command(about = "Remove sequence based on IDs", name = "remove")]
    Remove(SequenceRemoveArgs),
    #[command(
        about = "Batch renaming sequence IDs across multiple alignments",
        name = "stats"
    )]
    Rename(SequenceRenameArgs),
    #[command(about = "Translate DNA to amino acid sequences", name = "translate")]
    Translate(SequenceTranslateArgs),
}

#[derive(Args)]
pub(crate) struct RawStatArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
}

#[derive(Args)]
pub(crate) struct ContigStatArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
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
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Concat")]
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
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Convert")]
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
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Filter")]
    pub(crate) output: String,
    #[arg(long = "concat", help = "Concat filtered alignments")]
    pub(crate) concat: bool,
    #[arg(long = "len", help = "Filter by sequence length")]
    pub(crate) len: Option<usize>,
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
    #[arg(short, long, help = "Input path", default_value = "SEGUL-Split")]
    pub(crate) input: PathBuf,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[command(flatten)]
    pub(crate) out_fmt: CommonSeqOutput,
    #[command(flatten)]
    pub(crate) partition: CommonConcatArgs,
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Split")]
    pub(crate) output: PathBuf,
    #[arg(short = 'I', long = "input-partition", help = "Input partition file")]
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
        default_value = "nexus",
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
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Summary")]
    pub(crate) output: PathBuf,
    #[arg(long = "prefix", help = "Specify prefix for output files")]
    pub(crate) prefix: Option<String>,
    #[arg(
        long = "interval",
        help = "Specify interval value for counting data matrix completeness",
        default_value = "5",
        value_parser = builder::PossibleValuesParser::new(["1", "2", "5", "10"]),
    )]
    pub(crate) interval: usize,
    #[arg(long = "per-locus", help = "Generate summary statistic for each locus")]
    pub(crate) per_locus: bool,
}

#[derive(Args)]
pub(crate) struct PartitionArgs {
    pub(crate) input: Option<String>,
    #[arg(short, long, help = "Input a path (allow wildcard)")]
    #[cfg(not(target_os = "windows"))]
    pub(crate) input: Option<Vec<PathBuf>>,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[command(flatten)]
    pub(crate) out_fmt: CommonSeqOutput,
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Partition")]
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
pub(crate) struct SequenceExtractArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[command(flatten)]
    pub(crate) out_fmt: CommonSeqOutput,
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Extract")]
    pub(crate) output: PathBuf,
    #[arg(
        long = "re",
        help = "Specify regular expression for extracting sequences",
        require_equals = true
    )]
    pub(crate) re: Option<String>,
    #[arg(
        long = "id", 
        help = "Specify sequence ID for extracting sequences",
        required_unless_present_any(["re", "file"]),
    )]
    pub(crate) id: Option<Vec<String>>,
    #[arg(
        long = "file", 
        help = "Specify file for extracting sequences", 
        conflicts_with_all(["re", "id"]),
    )]
    pub(crate) file: Option<PathBuf>,
}

#[derive(Args)]
pub(crate) struct SequenceIdArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[command(flatten)]
    pub(crate) out_fmt: CommonSeqOutput,
    #[arg(short, long, help = "Output path", default_value = "id")]
    pub(crate) output: PathBuf,
    #[arg(long = "map", help = "Map ID across all alignments")]
    pub(crate) map: bool,
}

#[derive(Args)]
pub(crate) struct SequenceRemoveArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[command(flatten)]
    pub(crate) out_fmt: CommonSeqOutput,
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Remove")]
    pub(crate) output: PathBuf,
    #[arg(
        long = "re",
        help = "Specify regular expression for removing sequences",
        require_equals = true
    )]
    pub(crate) re: Option<String>,
    #[arg(
        long = "id",
        help = "Specify sequence ID for removing sequences",
        required_unless_present("re")
    )]
    pub(crate) id: Option<Vec<String>>,
}

#[derive(Args)]
pub(crate) struct SequenceRenameArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) out_fmt: CommonSeqOutput,
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Rename")]
    pub(crate) output: PathBuf,
}

#[derive(Args)]
pub(crate) struct SequenceTranslateArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[command(flatten)]
    pub(crate) out_fmt: CommonSeqOutput,
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Translate")]
    pub(crate) output: PathBuf,
    #[arg(long = "show-tables", help = "Show supported NCBI translation tables")]
    pub(crate) show_tables: bool,
    #[arg(
        long = "reading-frame", 
        help = "Specify reading frame", 
        default_value = "1",
        value_parser = builder::PossibleValuesParser::new(["1", "2", "3", "6", "9", "12"]),
    )]
    pub(crate) reading_frame: usize,
    #[arg(long = "auto-read", help = "Automatically detect reading frame")]
    pub(crate) auto_read: bool,
    #[arg(
        long = "table",
        help = "Specify NCBI translation table",
        default_value = "1",
        value_name = "INT",
        value_parser = builder::PossibleValuesParser::new(
            [
                "1",
                "2",
                "3",
                "4",
                "5",
                "6",
                "9",
                "10",
                "11",
                "12",
                "13",
                "14",
                "16",
                "21",
                "22",
                "23",
                "24",
                "25",
                "26",
                "29",
                "30",
                "33",
            ])
    )]
    pub(crate) table: usize,
}

#[derive(Args)]
pub(crate) struct IOArgs {
    #[arg(
        short,
        long,
        value_name = "PATH",
        help = "Input a directory",
        required_unless_present("input")
    )]
    pub(crate) dir: Option<String>,
    #[arg(short, long, help = "Input a path (allow wildcard)")]
    #[cfg(target_os = "windows")]
    pub(crate) input: Option<String>,
    #[arg(short, long, help = "Input a path (allow wildcard)")]
    #[cfg(not(target_os = "windows"))]
    pub(crate) input: Option<Vec<PathBuf>>,
    #[arg(long, help = "Force overwriting existing output files/directory")]
    pub(crate) force: bool,
}

#[derive(Args)]
pub(crate) struct CommonSeqOutput {
    #[arg(
        short = 'F',
        long = "output-format",
        help = "Specify output format",
        default_value = "nexus",
        value_parser = builder::PossibleValuesParser::new(
            ["fasta","nexus","phylip","fasta-int", "nexus-int", "phylip-int"]),
    )]
    pub(crate) output_fmt: String,
}

#[derive(Args)]
pub(crate) struct CommonSeqInput {
    #[arg(
        short = 'f',
        long = "input-format",
        value_name = "SEQUENCE INPUT FORMAT",
        help = "Specify input format",
        default_value = "auto",
        value_parser = builder::PossibleValuesParser::new(["auto","fasta","nexus","phylip"]),
    )]
    pub(crate) input_fmt: String,
    #[arg(
        long = "datatype",
        help = "Specify sequence datatype",
        default_value = "dna",
        value_parser = builder::PossibleValuesParser::new(["dna", "aa", "ignore"]),
    )]
    pub(crate) datatype: String,
}

#[derive(Args)]
pub(crate) struct CommonConcatArgs {
    #[arg(
        short = 'p',
        long = "partition-format",
        help = "Specify partition output format",
        default_value = "nexus",
        value_parser = builder::PossibleValuesParser::new(["charset", "nexus", "raxml"]),
    )]
    pub(crate) part_fmt: String,
    #[arg(long = "codon", help = "Set as a codon model partition format")]
    pub(crate) codon: bool,
    #[arg(long = "prefix", help = "Specify prefix for output files")]
    pub(crate) prefix: Option<PathBuf>,
}
