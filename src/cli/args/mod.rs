use std::path::PathBuf;

use align::AlignmentSubcommand;
use align::PartitionSubcommand;
use clap::Args;
use clap::Parser;
use clap::Subcommand;
use clap::builder;
use clap::builder::TypedValueParser as _;
use genomics::ContigSubcommand;
use genomics::MafSubcommand;
use genomics::SeqReadSubcommand;
use segul_core::plot::color::{Colormap, CvdKind, Palette};
use segul_core::plot::spec::PlotFormat;

use super::args::sequence::SequenceSubcommand;
use crate::helper::logger;
use clap::{crate_authors, crate_description, crate_name, crate_version};

pub(crate) mod align;
pub(crate) mod genomics;
pub(crate) mod sequence;

#[derive(Parser)]
#[command(name = crate_name!())]
#[command(version = crate_version!())]
#[command(author = crate_authors!())]
#[command(about = crate_description!(), long_about = None)]
pub(crate) struct Cli {
    #[command(subcommand)]
    pub(crate) subcommand: MainSubcommand,
    #[arg(
        long = "log",
        help = "Log file path",
        default_value = logger::LOG_FILE,
        global = true
    )]
    pub(crate) log: PathBuf,
}

#[derive(Subcommand)]
pub(crate) enum MainSubcommand {
    #[command(subcommand, about = "Sequence read analyses", name = "read")]
    RawRead(SeqReadSubcommand),
    #[command(subcommand, about = "Contiguous sequence analyses", name = "contig")]
    Contig(ContigSubcommand),
    #[command(
        subcommand,
        about = "Multi alignment format (MAF) analyses",
        name = "maf"
    )]
    Maf(MafSubcommand),
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

#[derive(Args)]
pub(crate) struct IOArgs {
    #[arg(
        short,
        long,
        value_name = "PATH",
        help = "Input a directory",
        required_unless_present("input"),
        conflicts_with("input")
    )]
    pub(crate) dir: Option<String>,
    #[arg(short, long, help = "Input a path (allow wildcard)")]
    #[cfg(target_os = "windows")]
    pub(crate) input: Option<String>,
    #[arg(short, long, help = "Input a path (allow wildcard)", num_args(0..))]
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

#[derive(Args)]
pub(crate) struct PlotFlags {
    #[arg(
        long = "plot",
        help = "Also write charts of the summary to <OUTPUT>/plots"
    )]
    pub(crate) plot: bool,
    #[arg(
        long = "plot-format",
        help = "Chart file format",
        requires = "plot",
        default_value_t = PlotFormat::Png,
        value_parser = builder::PossibleValuesParser::new(["png", "svg"])
            .map(|x| x.parse::<PlotFormat>().expect("Invalid plot format")),
    )]
    pub(crate) plot_format: PlotFormat,
    #[arg(
        long,
        value_name = "FILE",
        requires = "plot",
        help = "TrueType or OpenType font for PNG text [default: system fonts]"
    )]
    pub(crate) font: Option<PathBuf>,
    #[arg(
        long,
        help = "Colorblind-safe palette for categories",
        requires = "plot",
        default_value_t = Palette::OkabeIto,
        value_parser = builder::PossibleValuesParser::new(["okabe-ito", "tol-bright", "tol-muted"])
            .map(|x| x.parse::<Palette>().expect("Invalid palette")),
    )]
    pub(crate) palette: Palette,
    #[arg(
        long,
        help = "Colormap for heatmaps",
        requires = "plot",
        default_value_t = Colormap::Viridis,
        value_parser = builder::PossibleValuesParser::new(["viridis", "cividis", "magma", "inferno"])
            .map(|x| x.parse::<Colormap>().expect("Invalid colormap")),
    )]
    pub(crate) colormap: Colormap,
    #[arg(
        long,
        help = "Simulate a color vision deficiency in the charts",
        requires = "plot",
        value_parser = builder::PossibleValuesParser::new(
            ["protanopia", "deuteranopia", "tritanopia", "achromatopsia"])
            .map(|x| x.parse::<CvdKind>().expect("Invalid color vision deficiency")),
    )]
    pub(crate) cvd: Option<CvdKind>,
}
