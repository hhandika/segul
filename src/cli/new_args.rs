use std::path::PathBuf;

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
    AlignStats(AlignStatArgs),
}

#[derive(Subcommand)]
pub(crate) enum PartitionSubcommand {
    #[command(about = "Convert partition formats", name = "convert")]
    Convert(PartitionConvertArgs),
}

#[derive(Subcommand)]
pub(crate) enum SequenceSubcommand {
    #[command(about = "Parse sample ID across multiple alignments", name = "stats")]
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
    pub(crate) io: CommonArgs,
}

#[derive(Args)]
pub(crate) struct ContigStatArgs {
    #[command(flatten)]
    pub(crate) io: CommonArgs,
}

#[derive(Args)]
pub(crate) struct AlignConcatArgs {
    #[command(flatten)]
    pub(crate) io: CommonSeqArgs,
    #[command(flatten)]
    pub(crate) partition: CommonPartArgs,
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Concat")]
    pub(crate) output: String,
}

#[derive(Args)]
pub(crate) struct AlignConvertArgs {
    #[command(flatten)]
    pub(crate) io: CommonArgs,
    #[command(flatten)]
    pub(crate) format: CommonSeqArgs,
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Convert")]
    pub(crate) output: String,
    #[arg(long = "sort", help = "Sort sequences by IDs alphabetically")]
    pub(crate) sort: bool,
}

#[derive(Args)]
pub(crate) struct AlignFilterArgs {
    #[command(flatten)]
    pub(crate) io: CommonArgs,
    #[command(flatten)]
    pub(crate) format: CommonSeqArgs,
}

#[derive(Args)]
pub(crate) struct AlignSplitArgs {
    #[command(flatten)]
    pub(crate) io: CommonArgs,
    #[command(flatten)]
    pub(crate) format: CommonSeqArgs,
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Split")]
    pub(crate) output: String,
}

#[derive(Args)]
pub(crate) struct AlignStatArgs {
    #[command(flatten)]
    pub(crate) io: CommonArgs,
    #[command(flatten)]
    pub(crate) format: CommonSeqArgs,
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Stats")]
    pub(crate) output: String,
}

#[derive(Args)]
pub(crate) struct PartitionConvertArgs {
    #[command(flatten)]
    pub(crate) io: CommonArgs,
    #[command(flatten)]
    pub(crate) format: CommonSeqArgs,
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Partition")]
    pub(crate) output: String,
}

#[derive(Args)]
pub(crate) struct SequenceIdArgs {
    #[command(flatten)]
    pub(crate) io: CommonArgs,
    #[command(flatten)]
    pub(crate) format: CommonSeqArgs,
    #[arg(short, long, help = "Output path", default_value = "id")]
    pub(crate) output: String,
}

#[derive(Args)]
pub(crate) struct SequenceRemoveArgs {
    #[command(flatten)]
    pub(crate) io: CommonArgs,
    #[command(flatten)]
    pub(crate) format: CommonSeqArgs,
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Remove")]
    pub(crate) output: String,
}

#[derive(Args)]
pub(crate) struct SequenceRenameArgs {
    #[command(flatten)]
    pub(crate) io: CommonArgs,
    #[command(flatten)]
    pub(crate) format: CommonSeqArgs,
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Rename")]
    pub(crate) output: String,
}

#[derive(Args)]
pub(crate) struct SequenceTranslateArgs {
    #[command(flatten)]
    pub(crate) io: CommonArgs,
    #[command(flatten)]
    pub(crate) format: CommonSeqArgs,
    #[arg(short, long, help = "Output path", default_value = "SEGUL-Translate")]
    pub(crate) output: String,
}

#[derive(Args)]
pub(crate) struct CommonArgs {
    #[arg(
        short,
        long,
        value_name = "PATH",
        help = "Input a directory",
        required_unless_present("input")
    )]
    pub(crate) dir: Option<String>,
    #[arg(short, long, help = "Input a path (allow wildcard)")]
    pub(crate) input: Option<Vec<PathBuf>>,
    #[arg(long, help = "Force overwriting existing output files/directory")]
    pub(crate) force: bool,
}

#[derive(Args)]
pub(crate) struct CommonSeqArgs {
    #[arg(
        short = 'f',
        long = "input-format",
        value_name = "SEQUENCE FORMAT",
        help = "Specify input format",
        default_value = "auto",
        value_parser = builder::PossibleValuesParser::new(["auto","fasta","nexus","phylip"]),
    )]
    pub(crate) input_format: String,
    #[arg(
        short = 'o',
        long = "output-format",
        help = "Specify output format",
        default_value = "auto",
        value_parser = builder::PossibleValuesParser::new(["fasta","nexus","phylip","fasta-int", "nexus-int", "phylip-int"]),
    )]
    pub(crate) output_format: String,
    #[arg(
        long = "datatype",
        help = "Specify sequence datatype",
        default_value = "auto",
        value_parser = builder::PossibleValuesParser::new(["dna", "aa", "ignore"]),
    )]
    pub(crate) data_type: String,
}

#[derive(Args)]
pub(crate) struct CommonPartArgs {
    #[arg(short = 'p', long = "partition", help = "Specify partition format")]
    partition: String,
    #[arg(long = "codon", help = "Set as a codon model partition format")]
    codon: bool,
    #[arg(long = "prefix", help = "Specify prefix for output files")]
    prefix: Option<String>,
}
