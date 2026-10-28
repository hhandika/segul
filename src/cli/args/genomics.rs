use std::path::PathBuf;

use clap::builder::TypedValueParser as _;
use clap::{Args, Subcommand, builder};

use crate::helper::types::{ContigFmt, PhredOffset, SeqReadFmt, SummaryMode};
use crate::stats::fastq::DEFAULT_DUP_SKETCH_SIZE;

use super::{IOArgs, PlotFlags};

#[derive(Subcommand)]
pub(crate) enum SeqReadSubcommand {
    #[command(about = "Compute sequence read statistics", name = "summary")]
    RawSummary(SeqReadSummaryArgs),
}

#[derive(Subcommand)]
pub(crate) enum ContigSubcommand {
    #[command(about = "Compute contig statistics", name = "summary")]
    ContigSummary(ContigSummaryArgs),
}

#[derive(Subcommand)]
pub(crate) enum MafSubcommand {
    #[command(about = "Convert MAF blocks to per-locus alignments", name = "convert")]
    Convert(MafConvertArgs),
    #[command(
        about = "Summarize MAF blocks, species, and reference coverage",
        name = "summary"
    )]
    Summary(MafSummaryArgs),
}

#[derive(Args)]
pub(crate) struct SeqReadSummaryArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[arg(
        short = 'f',
        long ="input-format",
        help = "Specify input format",
        default_value_t = SeqReadFmt::Auto,
        value_parser =
            builder::PossibleValuesParser::new(["auto","fastq","gzip"])
            .map(|x| x.parse::<SeqReadFmt>().expect("Invalid input format")),
    )]
    pub(crate) input_format: SeqReadFmt,
    #[arg(
        long = "mode",
        help = "Summary mode",
        default_value_t = SummaryMode::Default,
        value_parser =
            builder::PossibleValuesParser::new(["minimal", "default", "complete"])
            .map(|x| x.parse::<SummaryMode>().unwrap()))]
    pub(crate) mode: SummaryMode,
    #[arg(
        long = "phred",
        help = "Quality score offset. Auto detects it from the first reads",
        default_value_t = PhredOffset::Auto,
        value_parser =
            builder::PossibleValuesParser::new(["auto", "33", "64"])
            .map(|x| x.parse::<PhredOffset>().unwrap()))]
    pub(crate) phred: PhredOffset,
    #[arg(
        long = "dup-sketch-size",
        help = "Distinct sequences kept to estimate duplication in complete mode. \
            Larger is more accurate and uses about 30 bytes per sequence per thread. \
            0 counts every sequence exactly",
        default_value_t = DEFAULT_DUP_SKETCH_SIZE
    )]
    pub(crate) dup_sketch_size: usize,
    #[arg(
        short = 'o',
        long = "output",
        help = "Output path",
        default_value = "Read-Summary"
    )]
    pub(crate) output: PathBuf,
    #[arg(long = "prefix", help = "Specify prefix for output files")]
    pub(crate) prefix: Option<String>,
    #[command(flatten)]
    pub(crate) plot: PlotFlags,
}

#[derive(Args)]
pub(crate) struct ContigSummaryArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[arg(
        short = 'f',
        long ="input-format",
        help = "Specify input format",
        default_value_t = ContigFmt::Auto,
        value_parser =
            builder::PossibleValuesParser::new(["auto","fasta","gzip"])
            .map(|x| x.parse::<ContigFmt>().expect("Invalid input format")),
    )]
    pub(crate) input_format: ContigFmt,
    #[arg(
        short = 'o',
        long = "output",
        help = "Output path",
        default_value = "Contig-Summary"
    )]
    pub(crate) output: PathBuf,
    #[arg(long = "prefix", help = "Specify prefix for output files")]
    pub(crate) prefix: Option<String>,
    #[command(flatten)]
    pub(crate) plot: PlotFlags,
}

#[derive(Args)]
pub(crate) struct MafCommonArgs {
    #[arg(
        long = "ref-species",
        value_name = "NAME",
        help = "Reference species. Defaults to the first species in the first block"
    )]
    pub(crate) ref_species: Option<String>,
    #[arg(
        long = "full-source-name",
        help = "Use full MAF source names (e.g. hg18.chr7) as sample names"
    )]
    pub(crate) full_source_name: bool,
    #[arg(
        long = "species",
        value_name = "LIST",
        conflicts_with = "species_file",
        help = "Keep only these species (comma-separated)"
    )]
    pub(crate) species: Option<String>,
    #[arg(
        long = "species-file",
        value_name = "PATH",
        help = "Keep only the species listed in a file, one per line"
    )]
    pub(crate) species_file: Option<PathBuf>,
}

#[derive(Args)]
pub(crate) struct MafConvertArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) common: MafCommonArgs,
    #[arg(
        long,
        value_name = "PATH",
        help = "BED file of loci (.bed or .bed.gz). Blocks are assigned to every locus they overlap"
    )]
    pub(crate) bed: Option<PathBuf>,
    #[arg(
        long = "group-by",
        help = "Group blocks without a BED file: one file per block or per chromosome",
        default_value = "block",
        value_parser = builder::PossibleValuesParser::new(["block", "chrom"]),
        conflicts_with = "bed"
    )]
    pub(crate) group_by: String,
    #[arg(
        long = "no-clip",
        requires = "bed",
        help = "Keep whole blocks instead of clipping to loci"
    )]
    pub(crate) no_clip: bool,
    #[arg(
        long = "strand-aware",
        requires = "bed",
        help = "Reverse complement loci on the minus strand of the BED file"
    )]
    pub(crate) strand_aware: bool,
    #[arg(
        long = "write-unmatched",
        requires = "bed",
        help = "Write blocks that match no locus to the unmatched directory"
    )]
    pub(crate) write_unmatched: bool,
    #[arg(
        long = "min-block-len",
        value_name = "BASES",
        help = "Skip blocks with fewer reference bases"
    )]
    pub(crate) min_block_len: Option<usize>,
    #[arg(
        long = "min-species",
        value_name = "COUNT",
        help = "Skip blocks with fewer species, counting the reference"
    )]
    pub(crate) min_species: Option<usize>,
    #[arg(
        long = "min-score",
        value_name = "SCORE",
        allow_negative_numbers = true,
        help = "Skip blocks with a lower score"
    )]
    pub(crate) min_score: Option<f64>,
    #[arg(
        short = 'F',
        long = "output-format",
        help = "Specify output format",
        default_value = "fasta",
        value_parser = builder::PossibleValuesParser::new(
            ["fasta", "nexus", "phylip", "fasta-int", "nexus-int", "phylip-int"]),
    )]
    pub(crate) output_fmt: String,
    #[arg(
        short = 'o',
        long = "output",
        help = "Output path",
        default_value = "Maf-Convert"
    )]
    pub(crate) output: PathBuf,
    #[arg(
        long = "prefix",
        help = "Name of the locus report CSV written with --bed"
    )]
    pub(crate) prefix: Option<String>,
}

#[derive(Args)]
pub(crate) struct MafSummaryArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) common: MafCommonArgs,
    #[arg(
        short = 'o',
        long = "output",
        help = "Output path",
        default_value = "Maf-Summary"
    )]
    pub(crate) output: PathBuf,
    #[arg(long = "prefix", help = "Specify prefix for output files")]
    pub(crate) prefix: Option<String>,
    #[command(flatten)]
    pub(crate) plot: PlotFlags,
}
