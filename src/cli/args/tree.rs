use std::path::PathBuf;

use clap::builder::TypedValueParser as _;
use clap::{Args, Subcommand, builder};

use segul_core::core::tree::estimate::{
    BranchSwap, GapMode, ParsimonyBranchLengths, ParsimonySearch, TreeMethod,
};
use segul_core::parser::tree::{TreeFmt, TreeOutputFmt};
use segul_core::stats::distance::DistanceModel;

use super::{CommonSeqInput, IOArgs};

#[derive(Subcommand)]
pub(crate) enum TreeSubcommand {
    #[command(about = "Convert phylogenetic tree formats", name = "convert")]
    Convert(TreeConvertArgs),
    #[command(
        about = "Infer phylogenetic trees from sequence alignments",
        name = "infer"
    )]
    Infer(TreeInferArgs),
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

#[derive(Args)]
pub(crate) struct TreeInferArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[arg(
        short = 'm',
        long = "method",
        value_name = "METHOD",
        help = "Tree estimation method",
        default_value = "nj",
        value_parser = builder::PossibleValuesParser::new([
            "nj", "bionj", "upgma", "mp", "parsimony",
        ])
        .map(|x| x.parse::<TreeMethod>().expect("Invalid tree method")),
    )]
    pub(crate) method: TreeMethod,
    #[arg(
        long = "model",
        value_name = "MODEL",
        help = "Distance model for NJ, BIONJ, and UPGMA",
        default_value = "jc69",
        value_parser = builder::PossibleValuesParser::new(["p", "jc69", "k2p", "tn93", "poisson"])
            .map(|x| x.parse::<DistanceModel>().expect("Invalid distance model")),
    )]
    pub(crate) model: DistanceModel,
    #[arg(
        long = "search",
        value_name = "SEARCH",
        help = "Maximum parsimony search mode",
        default_value = "auto",
        value_parser = builder::PossibleValuesParser::new(["auto", "exact", "heuristic"])
            .map(|x| x.parse::<ParsimonySearch>().expect("Invalid parsimony search mode")),
    )]
    pub(crate) search: ParsimonySearch,
    #[arg(
        long = "swap",
        value_name = "SWAP",
        help = "Branch-swapping algorithm for heuristic parsimony",
        default_value = "tbr",
        value_parser = builder::PossibleValuesParser::new(["none", "nni", "spr", "tbr"])
            .map(|x| x.parse::<BranchSwap>().expect("Invalid branch swapping method")),
    )]
    pub(crate) swap: BranchSwap,
    #[arg(
        long = "replicates",
        value_name = "N",
        help = "Random addition replicates for heuristic parsimony",
        default_value_t = 10,
        value_parser = builder::RangedU64ValueParser::<usize>::new().range(1..),
    )]
    pub(crate) replicates: usize,
    #[arg(
        long = "gap-mode",
        value_name = "MODE",
        help = "How gaps are treated in maximum parsimony",
        default_value = "missing",
        value_parser = builder::PossibleValuesParser::new(["missing", "fifth"])
            .map(|x| x.parse::<GapMode>().expect("Invalid gap mode")),
    )]
    pub(crate) gap_mode: GapMode,
    #[arg(
        long = "mp-branch-lengths",
        value_name = "MODE",
        help = "Branch length reconstruction for maximum parsimony",
        default_value = "acctran",
        value_parser = builder::PossibleValuesParser::new(["acctran", "deltran", "none"])
            .map(|x| {
                x.parse::<ParsimonyBranchLengths>()
                    .expect("Invalid MP branch length mode")
            }),
    )]
    pub(crate) mp_branch_lengths: ParsimonyBranchLengths,
    #[arg(
        short = 'b',
        long = "bootstrap",
        value_name = "N",
        help = "Number of non-parametric bootstrap replicates",
        default_value_t = 0
    )]
    pub(crate) bootstrap: usize,
    #[arg(
        long = "seed",
        value_name = "SEED",
        help = "Random seed for parsimony search and bootstrapping",
        default_value_t = 42
    )]
    pub(crate) seed: u64,
    #[arg(
        long = "outgroup",
        value_name = "TAXA",
        value_delimiter = ',',
        help = "Comma-separated outgroup taxon or clade for rooting",
        conflicts_with = "midpoint"
    )]
    pub(crate) outgroup: Option<Vec<String>>,
    #[arg(
        long = "midpoint",
        help = "Root the estimated tree at the longest path midpoint",
        conflicts_with = "outgroup"
    )]
    pub(crate) midpoint: bool,
    #[arg(
        short = 'F',
        long = "output-format",
        value_name = "TREE OUTPUT FORMAT",
        help = "Specify output tree format",
        default_value = "newick",
        value_parser = builder::PossibleValuesParser::new(["newick", "nexus", "phyloxml", "nexml"])
            .map(|x| x.parse::<TreeOutputFmt>().expect("Invalid tree output format")),
    )]
    pub(crate) output_fmt: TreeOutputFmt,
    #[arg(short, long, help = "Output path", default_value = "Tree-Infer")]
    pub(crate) output: PathBuf,
    #[arg(long = "prefix", help = "Specify prefix for output files")]
    pub(crate) prefix: Option<String>,
    #[arg(
        long = "single-file",
        help = "Write all estimated trees into a single multi-tree file"
    )]
    pub(crate) single_file: bool,
    #[arg(
        long = "no-translate",
        help = "Omit NEXUS TRANSLATE block when writing NEXUS trees"
    )]
    pub(crate) no_translate: bool,
    #[arg(long = "gzip", help = "Gzip-compress output tree files")]
    pub(crate) gzip: bool,
}
