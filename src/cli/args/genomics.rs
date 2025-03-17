

use std::path::PathBuf;

use clap::{builder, Args, Subcommand};
use clap::builder::TypedValueParser as _;

use crate::helper::types::{ContigFmt, SeqReadFmt, SummaryMode};

use super::IOArgs;


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
    #[command(about = "Convert genomic files to other formats", name = "maf")]
    Maf(MafConvertArgs),
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
    #[arg(short = 'o', long = "output", help = "Output path", default_value = "Read-Summary")]
    pub(crate) output: PathBuf,
    #[arg(long = "prefix", help = "Specify prefix for output files")]
    pub(crate) prefix: Option<String>,
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
    #[arg(short = 'o', long = "output", help = "Output path", default_value = "Contig-Summary")]
    pub(crate) output: PathBuf,
    #[arg(long = "prefix", help = "Specify prefix for output files")]
    pub(crate) prefix: Option<String>,
}

#[derive(Args)]
pub(crate) struct MafConvertArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[arg(long, help = "Path to the source of reference names")]
    pub(crate) reference: PathBuf,
    #[arg(long, help = "Source of names is a bed file")]
    pub(crate) from_bed: bool,
    #[arg(
        short = 't', 
        long = "output-format", 
        help = "Specify output format", 
        default_value = "fasta-int",
        value_parser = 
            builder::PossibleValuesParser::new(["fasta", "phylip","fasta-int", "phylip-int"]),
    )]
    pub(crate) output_fmt: String,
    #[arg(short = 'o', long = "output", help = "Output path", default_value = "Genomic-Convert")]
    pub(crate) output: PathBuf,
    #[arg(long = "prefix", help = "Specify prefix for output files")]
    pub(crate) prefix: Option<String>,
}
