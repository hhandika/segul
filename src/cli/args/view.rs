use clap::builder::TypedValueParser as _;
use clap::{Args, Subcommand, builder};

use super::IOArgs;

#[derive(Subcommand)]
pub(crate) enum ViewSubcommand {
    #[command(
        about = "View alignments in an interactive terminal viewer",
        name = "align"
    )]
    Align(ViewAlignArgs),
}

#[derive(Args)]
pub(crate) struct ViewAlignArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
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
        help = "Specify sequence datatype. Ignore detects DNA or amino acids from the residues",
        default_value = "ignore",
        value_parser = builder::PossibleValuesParser::new(["dna", "aa", "ignore"]),
    )]
    pub(crate) datatype: String,
    #[arg(
        long = "theme",
        help = "Residue colors for a dark or light terminal",
        default_value = "dark",
        value_parser = builder::PossibleValuesParser::new(["dark", "light"]),
    )]
    pub(crate) theme: String,
    #[arg(
        long = "chunk",
        help = "Sites per chunk",
        default_value_t = 1000,
        value_parser = builder::PossibleValuesParser::new(["250", "500", "1000", "2000"])
            .map(|s| s.parse::<usize>().expect("Failed parsing chunk size")),
    )]
    pub(crate) chunk: usize,
}
