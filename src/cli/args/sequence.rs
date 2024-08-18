use std::path::PathBuf;

use clap::Subcommand;

use clap::builder::TypedValueParser as _;
use clap::{builder, Args};

use super::{CommonSeqInput, CommonSeqOutput, IOArgs};


#[derive(Subcommand)]
pub(crate) enum SequenceSubcommand {
    #[command(about = "Add sequences to alignments", name = "add")]
    Add(SequenceAddArgs),
    #[command(about = "Extract sequence from alignments", name = "extract")]
    Extract(SequenceExtractArgs),
    #[command(about = "Filter sequence based on selected criteria", name = "filter")]
    Filter(SequenceFilterArgs),
    #[command(about = "Parse sample ID across multiple alignments", name = "id")]
    Id(SequenceIdArgs),
    #[command(about = "Remove sequence based on IDs", name = "remove")]
    Remove(SequenceRemoveArgs),
    #[command(
        about = "Batch renaming sequence IDs across multiple alignments",
        name = "rename"
    )]
    Rename(SequenceRenameArgs),
    #[command(about = "Translate DNA to amino acid sequences", name = "translate")]
    Translate(SequenceTranslateArgs),
}

#[derive(Args)]
pub(crate) struct SequenceAddArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[arg(
        short = 'F',
        long = "output-format",
        help = "Specify output format",
        default_value = "fasta-int",
        value_parser = builder::PossibleValuesParser::new(
            ["fasta", "fasta-int"]),
    )]
    pub(crate) output_fmt: String,
    #[arg(short, long, help = "Output path", default_value = "Sequence-Addition")]
    pub(crate) output: PathBuf,
    #[arg(
        long,
        help = "Input directory for destination sequences",
        required_unless_present("destination_input"),
    )]
    pub(crate) destination_dir: Option<String>,
    #[arg(
        long,
        help = "Input a path (allow wildcard) for destination sequences",
    )]
    pub(crate) destination_input: Option<Vec<PathBuf>>,
    #[arg(
        long = "destination-format",
        help = "Specify destination sequence format",
        default_value = "auto",
        value_parser = builder::PossibleValuesParser::new(
            ["auto", "fasta", "nexus", "phylip"]
        ),
    )]
    pub(crate) destination_fmt: String,
}


#[derive(Args)]
pub(crate) struct SequenceExtractArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[command(flatten)]
    pub(crate) out_fmt: CommonSeqOutput,
    #[arg(short, long, help = "Output path", default_value = "Sequence-Extract")]
    pub(crate) output: PathBuf,
    #[arg(
        long = "re",
        help = "Specify regular expression for extracting sequences",
        require_equals = true,
        conflicts_with_all(["id", "file"]),
    )]
    pub(crate) re: Option<String>,
    #[arg(
        long = "id", 
        help = "Input sequence IDa separated by semicolon",
        required_unless_present_any(["re", "file"]),
        require_equals = true,
    )]
    pub(crate) id: Option<String>,
    #[arg(
        long = "file", 
        help = "Specify file for extracting sequences", 
        conflicts_with_all(["re", "id"]),
    )]
    pub(crate) file: Option<PathBuf>,
}

#[derive(Args)]
pub(crate) struct SequenceFilterArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[command(flatten)]
    pub(crate) out_fmt: CommonSeqOutput,
    #[arg(short, long, help = "Output directory", default_value = "Sequence-Filter")]
    pub(crate) output: PathBuf,
    #[arg(long = "min-length", help = "Filter by minimal sequence length")]
    pub(crate) min_len: Option<usize>,
    #[arg(
        long = "max-gap",
        help = "Filter by maximum gap percentage",
    )]
    pub(crate) max_gap: Option<f64>,
}

#[derive(Args)]
pub(crate) struct SequenceIdArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[command(flatten)]
    pub(crate) out_fmt: CommonSeqOutput,
    #[arg(short, long, help = "Output path", default_value = "Sequence-ID")]
    pub(crate) output: PathBuf,
    #[arg(short, long, help = "Prefix for filename")]
    pub(crate) prefix: Option<String>,
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
    #[arg(short, long, help = "Output path", default_value = "Sequence-Remove")]
    pub(crate) output: PathBuf,
    #[arg(
        long = "re",
        help = "Input regular expression for removing sequences",
        require_equals = true
    )]
    pub(crate) re: Option<String>,
    #[arg(
        long = "id",
        help = "Input sequence ID separated by semicolon",
        required_unless_present("re"),
        require_equals = true
    )]
    pub(crate) id: Option<String>,
}

#[derive(Args)]
pub(crate) struct SequenceRenameArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[command(flatten)]
    pub(crate) out_fmt: CommonSeqOutput,
    #[arg(short, long, help = "Output path", default_value = "Sequence-Rename")]
    pub(crate) output: PathBuf,
    #[arg(long = "dry-run", help = "Dry run")]
    pub(crate) dry_run: bool,
    #[arg(
        long = "replace-id", 
        help = "Rename using input IDs in a file",         
    )]
    pub(crate) replace_id: Option<PathBuf>,
    #[arg(long = "remove", help = "Remove matching input string")]
    pub(crate) remove: Option<String>,
    #[arg(
        long = "remove-re",
        help = "Remove first found matching input regular expression",
        require_equals = true
    )]
    pub(crate) remove_re: Option<String>,
    #[arg(
        long = "replace-from",
        help = "Replace matching input string with the output string",
        require_equals = true
    )]
    pub(crate) replace_from: Option<String>,
    #[arg(
        long = "remove-re-all",
        help = "Remove all found matching input regular expression",
        require_equals = true
    )]
    pub(crate) remove_re_all: Option<String>,
    #[arg(
        long = "replace-from-re",
        help = "Replace first found matching input regular expression with the output string",
        require_equals = true
    )]
    pub(crate) replace_from_re: Option<String>,
    #[arg(
        long = "replace-from-re-all",
        help = "Replace all found matching input regular expression with the output string",
        require_equals = true
    )]
    pub(crate) replace_from_re_all: Option<String>,
    #[arg(
        long = "replace-to",
        help = "Replace matching input string with the output string",
        require_equals = true
    )]
    pub(crate) replace_to: Option<String>,
}

#[derive(Args)]
pub(crate) struct SequenceTranslateArgs {
    #[command(flatten)]
    pub(crate) io: IOArgs,
    #[command(flatten)]
    pub(crate) in_fmt: CommonSeqInput,
    #[command(flatten)]
    pub(crate) out_fmt: CommonSeqOutput,
    #[arg(short, long, help = "Output path", default_value = "Sequence-Translate")]
    pub(crate) output: PathBuf,
    #[arg(long = "show-tables", help = "Show supported NCBI translation tables")]
    pub(crate) show_tables: bool,
    #[arg(
        long = "reading-frame", 
        help = "Specify reading frame", 
        default_value_t = 1,
        value_parser = builder::PossibleValuesParser::new(["1", "2", "3"])
            .map(|x| x.parse::<usize>().unwrap_or(1)),
    )]
    pub(crate) reading_frame: usize,
    #[arg(
        long = "auto-read",
        help = "Automatically detect reading frame",
    )]
    pub(crate) auto_read: bool,
    #[arg(
        long = "table",
        help = "Specify NCBI translation table",
        default_value = "1",
        value_name = "NUMBER",
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
    pub(crate) table: String,
}
