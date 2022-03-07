mod args;
mod concat;
mod convert;
mod extract;
mod filter;
mod id;
mod rename;
mod split;
mod summarize;
mod translate;

use std::ffi::OsStr;
use std::fs;
use std::io::Result;

use std::path::{Path, PathBuf};

use clap::ArgMatches;
use dialoguer::{theme::ColorfulTheme, Confirm};
use glob::glob;
use log::LevelFilter;
use log4rs::append::console::ConsoleAppender;
use log4rs::append::file::FileAppender;
use log4rs::config::{Appender, Config, Root};
use log4rs::encode::pattern::PatternEncoder;

use crate::cli::concat::ConcatParser;
use crate::cli::convert::ConvertParser;
use crate::cli::extract::ExtractParser;
use crate::cli::filter::FilterParser;
use crate::cli::id::IdParser;
use crate::cli::rename::RenameParser;
use crate::cli::split::SplitParser;
use crate::cli::summarize::SummaryParser;
use crate::cli::translate::TranslateParser;

use crate::helper::finder::Files;
use crate::helper::types::{DataType, InputFmt, OutputFmt, PartitionFmt};
use crate::helper::utils;

pub const LOG_FILE: &str = "segul.log";

pub fn parse_cli(version: &str) {
    let args = args::get_args(version);
    setup_logger().expect("Failed setting up a log file.");
    utils::print_welcome_text(version);
    match args.subcommand() {
        ("convert", Some(convert_matches)) => ConvertParser::new(convert_matches).convert(),
        ("concat", Some(concat_matches)) => ConcatParser::new(concat_matches).concat(),
        ("filter", Some(pick_matches)) => FilterParser::new(pick_matches).filter(),
        ("id", Some(id_matches)) => IdParser::new(id_matches).get_id(),
        ("extract", Some(extract_matches)) => ExtractParser::new(extract_matches).extract(),
        ("rename", Some(rename_matches)) => RenameParser::new(rename_matches).rename(),
        ("summary", Some(stats_matches)) => SummaryParser::new(stats_matches).stats(),
        ("translate", Some(trans_matches)) => TranslateParser::new(trans_matches).translate(),
        ("split", Some(split_matches)) => SplitParser::new(split_matches).split(),
        _ => unreachable!("Invalid subcommands!"),
    }
}

fn setup_logger() -> Result<()> {
    let log_dir = std::env::current_dir()?;
    let target = log_dir.join(LOG_FILE);
    let tofile = FileAppender::builder()
        .encoder(Box::new(PatternEncoder::new(
            "{d(%Y-%m-%d %H:%M:%S %Z)} - {l} - {m}\n",
        )))
        .build(target)?;

    let stdout = ConsoleAppender::builder()
        .encoder(Box::new(PatternEncoder::new("{m}\n")))
        .build();

    let config = Config::builder()
        .appender(Appender::builder().build("stdout", Box::new(stdout)))
        .appender(Appender::builder().build("logfile", Box::new(tofile)))
        .build(
            Root::builder()
                .appender("stdout")
                .appender("logfile")
                .build(LevelFilter::Info),
        )
        .expect("Failed building log configuration");

    log4rs::init_config(config).expect("Cannot initiate log configuration");

    Ok(())
}

macro_rules! check_output_path {
    ($type: ident, $execution: ident, $path: ident, $prompt: expr, $err_msg: expr) => {
        if $path.$type() {
            let selection = Confirm::with_theme(&ColorfulTheme::default())
                .with_prompt($prompt)
                .interact();
            match selection {
                Ok(yes) => {
                    if yes {
                        fs::$execution($path).expect($err_msg);
                        println!();
                    } else {
                        std::process::abort();
                    }
                }
                Err(err) => panic!("Failed parsing user input: {}", err),
            }
        }
    };
}

trait InputCli {
    fn parse_dir_input<'a>(&self, matches: &'a ArgMatches) -> &'a Path {
        Path::new(matches.value_of("dir").expect("Failed parsing a dir value"))
    }

    fn parse_input(&self, matches: &ArgMatches) -> Vec<PathBuf> {
        let inputs = matches
            .values_of("input")
            .expect("Failed parsing input values")
            .map(PathBuf::from)
            .collect::<Vec<PathBuf>>();
        if cfg!(windows) {
            let inputs = inputs
                .iter()
                .map(|t| OsStr::new(t).to_string_lossy())
                .collect::<Vec<_>>();
            let files: Vec<PathBuf> = inputs
                .iter()
                .flat_map(|i| {
                    glob(i)
                        .expect("Failed globbing files")
                        .filter_map(|ok| ok.ok())
                        .collect::<Vec<PathBuf>>()
                })
                .collect();
            assert!(!files.is_empty(), "Empty folders!");
            files
        } else {
            inputs
        }
    }

    fn get_files(&self, dir: &Path, input_fmt: &InputFmt) -> Vec<PathBuf> {
        Files::new(dir, input_fmt).get_files()
    }

    fn parse_input_fmt(&self, matches: &ArgMatches) -> InputFmt {
        let input_fmt = matches
            .value_of("input-format")
            .expect("Failed parsing an input format value");
        match input_fmt {
            "auto" => InputFmt::Auto,
            "fasta" => InputFmt::Fasta,
            "nexus" => InputFmt::Nexus,
            "phylip" => InputFmt::Phylip,
            _ => unreachable!("Unknown input format. Supported format: auto, fasta, nexus, phylip"),
        }
    }

    fn parse_datatype(&self, matches: &ArgMatches) -> DataType {
        let datatype = matches
            .value_of("datatype")
            .expect("Failed parsing dataype value");
        match datatype {
            "aa" => DataType::Aa,
            "dna" => DataType::Dna,
            "ignore" => DataType::Ignore,
            _ => unreachable!(),
        }
    }

    fn parse_overwrite_opts(&self, matches: &ArgMatches) -> bool {
        matches.is_present("overwrite")
    }
}

trait InputPrint {
    fn print_input<P: AsRef<Path>>(
        &self,
        input: &Option<P>,
        task_desc: &str,
        fcounts: usize,
        input_fmt: &InputFmt,
        datatype: &DataType,
    ) {
        if let Some(input) = input {
            log::info!("{:18}: {}", "Input dir", &input.as_ref().display());
        } else {
            log::info!("{:18}: {}", "Input path", "STDIN");
        }
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&fcounts));
        self.print_input_fmt(input_fmt);
        self.print_datatype(datatype);
        log::info!("{:18}: {}\n", "Task", task_desc);
    }

    fn print_datatype(&self, datatype: &DataType) {
        match datatype {
            DataType::Aa => log::info!("{:18}: {}", "Data type", "Amino Acid"),
            DataType::Dna => log::info!("{:18}: {}", "Data type", "DNA"),
            DataType::Ignore => log::info!("{:18}: {}", "Data type", "Ignore"),
        }
    }

    fn print_input_fmt(&self, input_fmt: &InputFmt) {
        match input_fmt {
            InputFmt::Auto => log::info!("{:18}: {}", "Input format", "Auto"),
            InputFmt::Fasta => log::info!("{:18}: {}", "Input format", "Fasta"),
            InputFmt::Nexus => log::info!("{:18}: {}", "Input format", "Nexus"),
            InputFmt::Phylip => log::info!("{:18}: {}", "Input format", "Phylip"),
        }
    }
}

trait OutputCli {
    fn parse_output<'a>(&self, matches: &'a ArgMatches) -> PathBuf {
        let output = matches
            .value_of("output")
            .expect("Failed parsing an output value");
        PathBuf::from(output)
    }

    fn parse_output_fmt(&self, matches: &ArgMatches) -> OutputFmt {
        let output_fmt = matches
            .value_of("output-format")
            .expect("Failed parsing ouput format");
        match output_fmt {
            "nexus" => OutputFmt::Nexus,
            "phylip" => OutputFmt::Phylip,
            "fasta" => OutputFmt::Fasta,
            "nexus-int" => OutputFmt::NexusInt,
            "fasta-int" => OutputFmt::FastaInt,
            "phylip-int" => OutputFmt::PhylipInt,
            _ => unreachable!("Please, specify the correct output format!"),
        }
    }

    fn check_output_file_exist(&self, path: &Path, overwrite: bool) {
        let rm_err_msg = "Failed removing existing output files";
        if overwrite {
            log::warn!("Overwriting output file: {}", path.display());
            fs::remove_file(path).expect(rm_err_msg);
        } else {
            let error_msg = format!("Output file already exists: {}. Remove?", path.display());
            check_output_path!(is_file, remove_file, path, error_msg, rm_err_msg);
        }
    }

    fn check_output_dir_exist(&self, path: &Path, overwrite: bool) {
        let rm_err_msg = "Failed removing a directory";
        if overwrite {
            log::warn!("Removing directory: {}", path.display());
            fs::remove_dir_all(path).expect(rm_err_msg);
        } else {
            let error_msg = format!("Output dir already exists: {}. Remove?", path.display());
            check_output_path!(is_dir, remove_dir_all, path, error_msg, rm_err_msg)
        }
    }
}

trait ConcatCli {
    fn parse_prefix(&self, matches: &ArgMatches, dir: &Path) -> PathBuf {
        if matches.is_present("prefix") {
            let prefix = matches
                .value_of("prefix")
                .expect("Failed parsing a prefix value");
            PathBuf::from(prefix)
        } else {
            PathBuf::from(dir)
        }
    }

    fn parse_partition_fmt(&self, matches: &ArgMatches) -> PartitionFmt {
        let part_fmt = matches
            .value_of("partition")
            .expect("Failed parsing partition format");
        if matches.is_present("codon") {
            self.parse_partition_fmt_codon(part_fmt)
        } else {
            self.parse_partition_fmt_std(part_fmt)
        }
    }

    fn parse_partition_fmt_std(&self, part_fmt: &str) -> PartitionFmt {
        match part_fmt {
            "nexus" => PartitionFmt::Nexus,
            "raxml" => PartitionFmt::Raxml,
            "charset" => PartitionFmt::Charset,
            _ => PartitionFmt::Nexus,
        }
    }

    fn parse_partition_fmt_codon(&self, part_fmt: &str) -> PartitionFmt {
        match part_fmt {
            "charset" => PartitionFmt::CharsetCodon,
            "nexus" => PartitionFmt::NexusCodon,
            "raxml" => PartitionFmt::RaxmlCodon,
            _ => PartitionFmt::NexusCodon,
        }
    }

    fn check_partition_format(&self, output_fmt: &OutputFmt, part_fmt: &PartitionFmt) {
        match output_fmt {
            OutputFmt::Nexus | OutputFmt::NexusInt => (),
            _ => {
                if let PartitionFmt::Charset | PartitionFmt::CharsetCodon = part_fmt {
                    panic!(
                        "Cannot write embedded-nexus partition 'charset' to non-nexus output. \
                Maybe you mean to write the partition to 'nexus' instead."
                    )
                }
            }
        }
    }
}
