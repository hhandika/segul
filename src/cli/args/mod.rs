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
use segul_core::color::{Colormap, CvdKind, Palette};
use segul_core::plot::spec::PlotFormat;
use segul_core::report::ChecksumMode;
use segul_core::report::methods::MethodsFormat;
use segul_core::report::verify::VerifyScope;

use super::args::sequence::SequenceSubcommand;
use super::args::view::ViewSubcommand;
use crate::helper::logger;
use clap::{crate_authors, crate_description, crate_name, crate_version};

pub(crate) mod align;
pub(crate) mod genomics;
pub(crate) mod sequence;
pub(crate) mod view;

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
    #[arg(
        long = "checksum",
        value_name = "MODE",
        help = "Compute SHA-256 checksums for recorded files [env: SEGUL_CHECKSUM]",
        value_parser = builder::PossibleValuesParser::new(["none", "outputs", "inputs", "all"])
            .map(|x| x.parse::<ChecksumMode>().expect("Invalid checksum mode")),
        global = true
    )]
    pub(crate) checksum: Option<ChecksumMode>,
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
    #[command(about = "Re-run a report.json or segul-workflow.json", name = "run")]
    Run(RunArgs),
    #[command(
        subcommand,
        about = "Verify checksums, export methods, and build workflows from reports",
        name = "report"
    )]
    Report(ReportSubcommand),
    #[command(
        subcommand,
        about = "Interactive viewers in the terminal",
        name = "view"
    )]
    View(ViewSubcommand),
}

#[derive(Args)]
pub(crate) struct RunArgs {
    #[arg(
        value_name = "FILE",
        help = "Path to report.json or segul-workflow.json"
    )]
    pub(crate) file: PathBuf,
    #[arg(short, long, value_name = "DIR", help = "Override output directory")]
    pub(crate) output: Option<PathBuf>,
    #[arg(
        long = "input-dir",
        value_name = "DIR",
        help = "Override input base directory"
    )]
    pub(crate) input_dir: Option<PathBuf>,
    #[arg(long, help = "Overwrite existing output directory")]
    pub(crate) force: bool,
    #[arg(
        long = "verify-inputs",
        help = "Verify recorded input hashes and sizes before re-running"
    )]
    pub(crate) verify_inputs: bool,
}

#[derive(Subcommand)]
pub(crate) enum ReportSubcommand {
    #[command(about = "Verify SHA-256 checksums and file sizes recorded in report.json")]
    Verify(ReportVerifyArgs),
    #[command(about = "Generate a publication-ready Methods section from reports or workflows")]
    Methods(ReportMethodsArgs),
    #[command(about = "Export a segul-workflow.json recipe from chained or ordered reports")]
    Workflow(ReportWorkflowArgs),
}

#[derive(Args)]
pub(crate) struct ReportVerifyArgs {
    #[arg(value_name = "REPORT", help = "Path to report.json")]
    pub(crate) report: PathBuf,
    #[arg(
        long,
        default_value_t = VerifyScope::All,
        value_parser = builder::PossibleValuesParser::new(["all", "inputs", "outputs"])
            .map(|x| x.parse::<VerifyScope>().expect("Invalid verify scope")),
        help = "Which manifest files to verify"
    )]
    pub(crate) scope: VerifyScope,
    #[arg(
        long = "input-dir",
        value_name = "DIR",
        help = "Override input base directory"
    )]
    pub(crate) input_dir: Option<PathBuf>,
    #[arg(
        long = "output-dir",
        value_name = "DIR",
        help = "Override output directory"
    )]
    pub(crate) output_dir: Option<PathBuf>,
    #[arg(long, help = "Print machine-readable JSON verification results")]
    pub(crate) json: bool,
}

#[derive(Args)]
pub(crate) struct ReportMethodsArgs {
    #[arg(
        value_name = "FILE",
        help = "Path to report.json or segul-workflow.json"
    )]
    pub(crate) report: PathBuf,
    #[arg(
        long,
        default_value_t = MethodsFormat::Markdown,
        value_parser = builder::PossibleValuesParser::new(["markdown", "text"])
            .map(|x| x.parse::<MethodsFormat>().expect("Invalid methods format")),
        help = "Output format for generated Methods text"
    )]
    pub(crate) format: MethodsFormat,
    #[arg(
        long = "no-chain",
        help = "Only describe the specified report instead of walking parent_reports"
    )]
    pub(crate) no_chain: bool,
    #[arg(
        short,
        long,
        value_name = "FILE",
        help = "Write output to a file instead of stdout"
    )]
    pub(crate) output: Option<PathBuf>,
}

#[derive(Args)]
pub(crate) struct ReportWorkflowArgs {
    #[arg(value_name = "REPORTS", num_args(1..), required = true, help = "Leaf report.json (auto-walks lineage) or ordered report.json files")]
    pub(crate) reports: Vec<PathBuf>,
    #[arg(long, help = "Workflow name")]
    pub(crate) name: Option<String>,
    #[arg(
        short,
        long,
        value_name = "FILE",
        default_value = "segul-workflow.json",
        help = "Output workflow recipe path"
    )]
    pub(crate) output: PathBuf,
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
        help = "Palette for categories",
        requires = "plot",
        default_value_t = Palette::Dark2,
        value_parser = builder::PossibleValuesParser::new(
            ["dark2", "okabe-ito", "tol-bright", "tol-muted", "set2", "paired", "set1"])
            .map(|x| x.parse::<Palette>().expect("Invalid palette")),
    )]
    pub(crate) palette: Palette,
    #[arg(
        long,
        help = "Colormap for heatmaps",
        requires = "plot",
        default_value_t = Colormap::Viridis,
        value_parser = builder::PossibleValuesParser::new([
            "viridis", "cividis", "magma", "inferno", "plasma", "blues", "greens",
            "ylgnbu", "ylorrd", "rdbu", "puor", "brbg",
        ])
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
