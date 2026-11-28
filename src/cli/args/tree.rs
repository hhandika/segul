use std::path::PathBuf;

use clap::builder::TypedValueParser as _;
use clap::{Args, Subcommand, builder};

use segul_core::core::tree::estimate::{
    BranchSwap, GapMode, ParsimonyBranchLengths, ParsimonySearch, TreeMethod,
};
use segul_core::core::tree::summary::{
    BranchLengthStrategy, BurninSpec, MccScoring, SummaryMethod, SupportFormat, SupportPlacement,
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
    #[command(about = "Build consensus and summary trees", name = "summary")]
    Summary(TreeSummaryArgs),
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

#[derive(Args)]
pub(crate) struct TreeSummaryArgs {
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
        short = 'm',
        long = "method",
        value_name = "METHOD",
        help = "Summary tree method",
        default_value = "mcc",
        value_parser = builder::PossibleValuesParser::new([
            "mcc", "mrc", "mre", "strict", "map", "target",
        ])
        .map(|x| x.parse::<SummaryMethod>().expect("Invalid tree summary method")),
    )]
    pub(crate) method: SummaryMethod,
    #[arg(
        long = "target",
        value_name = "FILE",
        help = "Target tree file to annotate when using --method target"
    )]
    pub(crate) target: Option<PathBuf>,
    #[arg(
        long = "burnin",
        value_name = "SPEC",
        help = "Per-file burn-in trees to discard (count, percentage like '10%', or fraction)",
        default_value = "0",
        value_parser = |s: &str| s.parse::<BurninSpec>(),
    )]
    pub(crate) burnin: BurninSpec,
    #[arg(
        long = "min-clade-prob",
        visible_alias = "threshold",
        value_name = "PROB",
        help = "Minimum clade frequency for MRC (>= 0.5, default 0.5) or MRE (default 0.0)"
    )]
    pub(crate) min_clade_prob: Option<f64>,
    #[arg(
        long = "rooted",
        help = "Treat input trees as rooted clades",
        conflicts_with = "unrooted"
    )]
    pub(crate) rooted: bool,
    #[arg(
        long = "unrooted",
        help = "Treat input trees as unrooted bipartitions",
        conflicts_with = "rooted"
    )]
    pub(crate) unrooted: bool,
    #[arg(
        long = "mcc-scoring",
        value_name = "SCORING",
        help = "Clade credibility scoring objective for MCC",
        default_value = "log-product",
        value_parser = builder::PossibleValuesParser::new(["log-product", "sum"])
            .map(|x| x.parse::<MccScoring>().expect("Invalid MCC scoring mode")),
    )]
    pub(crate) mcc_scoring: MccScoring,
    #[arg(
        long = "branch-lengths",
        value_name = "STRATEGY",
        help = "How to summarize branch lengths and node heights",
        default_value = "mean",
        value_parser = builder::PossibleValuesParser::new([
            "mean", "median", "ca", "common-ancestor", "keep", "none",
        ])
        .map(|x| x.parse::<BranchLengthStrategy>().expect("Invalid branch length strategy")),
    )]
    pub(crate) branch_lengths: BranchLengthStrategy,
    #[arg(
        long = "support-format",
        value_name = "FORMAT",
        help = "Format for clade support values",
        default_value = "proportion",
        value_parser = builder::PossibleValuesParser::new(["proportion", "percent"])
            .map(|x| x.parse::<SupportFormat>().expect("Invalid support format")),
    )]
    pub(crate) support_format: SupportFormat,
    #[arg(
        long = "support-placement",
        value_name = "PLACE",
        help = "Where support values are written on the summary tree",
        default_value = "auto",
        value_parser = builder::PossibleValuesParser::new(["auto", "annotation", "label", "both"])
            .map(|x| x.parse::<SupportPlacement>().expect("Invalid support placement")),
    )]
    pub(crate) support_placement: SupportPlacement,
    #[arg(
        long = "hpd",
        value_name = "LEVEL",
        help = "HPD credible interval level in (0, 1), or 0 to disable",
        default_value_t = 0.95
    )]
    pub(crate) hpd: f64,
    #[arg(
        long = "no-annotations",
        help = "Skip summarizing BEAST/FigTree/NHX node and branch annotations"
    )]
    pub(crate) no_annotations: bool,
    #[arg(
        long = "allow-partial-taxa",
        help = "Allow trees with differing taxon sets (decisive gCF for target, shared taxa otherwise)"
    )]
    pub(crate) allow_partial_taxa: bool,
    #[arg(
        long = "outgroup",
        value_name = "TAXA",
        value_delimiter = ',',
        help = "Comma-separated outgroup taxon or clade for post-summary rooting",
        conflicts_with = "midpoint"
    )]
    pub(crate) outgroup: Option<Vec<String>>,
    #[arg(
        long = "midpoint",
        help = "Root the summary tree at its longest path midpoint",
        conflicts_with = "outgroup"
    )]
    pub(crate) midpoint: bool,
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
    #[arg(short, long, help = "Output path", default_value = "Tree-Summary")]
    pub(crate) output: PathBuf,
    #[arg(long = "prefix", help = "Specify prefix for output files")]
    pub(crate) prefix: Option<String>,
    #[arg(
        long = "no-translate",
        help = "Omit NEXUS TRANSLATE block when writing NEXUS output"
    )]
    pub(crate) no_translate: bool,
    #[arg(
        long = "strip-annotations",
        help = "Strip bracket annotations and comments from the output tree"
    )]
    pub(crate) strip_annotations: bool,
    #[arg(long = "gzip", help = "Gzip-compress the output summary tree file")]
    pub(crate) gzip: bool,
}
