use std::path::PathBuf;

use clap::builder::TypedValueParser as _;
use clap::{Args, Subcommand, builder};

use segul_core::parser::tree::{TreeFmt, TreeOutputFmt};

use super::IOArgs;

#[derive(Subcommand)]
pub(crate) enum TreeSubcommand {
    #[command(about = "Convert phylogenetic tree formats", name = "convert")]
    Convert(TreeConvertArgs),
}

#[derive(Args)]
pub(crate) struct TreeConvertArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[arg(
        short = 'f',
        long = "input-format",
        value_name = "TREE INPUT FORMAT",
        help = "Specify input tree format",
        default_value = "auto",
        value_parser = builder::PossibleValuesParser::new([
            "auto", "newick", "enewick", "nexus", "phyloxml", "nexml",
        ])
        .map(|x| x.parse::<TreeFmt>().expect("Invalid tree input format")),
    )]
    pub(crate) input_fmt: TreeFmt,
    #[arg(
        short = 'F',
        long = "output-format",
        value_name = "TREE OUTPUT FORMAT",
        help = "Specify output tree format",
        default_value = "newick",
        value_parser = builder::PossibleValuesParser::new([
            "newick", "enewick", "nexus", "phyloxml", "nexml",
        ])
        .map(|x| x.parse::<TreeOutputFmt>().expect("Invalid tree output format")),
    )]
    pub(crate) output_fmt: TreeOutputFmt,
    #[arg(short, long, help = "Output path", default_value = "Tree-Convert")]
    pub(crate) output: PathBuf,
    #[arg(
        long = "no-translate",
        help = "Omit NEXUS TRANSLATE block when writing NEXUS trees"
    )]
    pub(crate) no_translate: bool,
    #[arg(
        long = "strip-annotations",
        help = "Strip bracket annotations and comments from output trees"
    )]
    pub(crate) strip_annotations: bool,
    #[arg(long = "gzip", help = "Gzip-compress output tree files")]
    pub(crate) gzip: bool,
}
