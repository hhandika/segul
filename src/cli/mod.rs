//! Command line interface for parsing and executing commands.
mod args;
mod cli;
mod concat;
mod convert;
mod extract;
mod filter;
mod id;
mod partition;
mod raw;
mod remove;
mod rename;
mod split;
mod summarize;
mod translate;

#[cfg(target_os = "windows")]
use glob::glob;
#[cfg(target_os = "windows")]
use std::ffi::OsStr;

use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use clap::Parser;
use colored::Colorize;
use dialoguer::{theme::ColorfulTheme, Confirm};

use crate::cli::args::Cli;
use crate::helper::finder::Files;
use crate::helper::types::{DataType, InputFmt, OutputFmt, PartitionFmt, RawReadFmt};
use crate::helper::{logger, utils};

/// Parse command line arguments and execute commands.
pub fn parse_cli() {
    let time = Instant::now();
    let args = Cli::parse();
    logger::setup_logger(&args.log).expect("Failed setting up a log file.");
    utils::print_welcome_text(clap::crate_version!());
    cli::match_cli_subcommand(&args.subcommand);
    log::info!("{:18}: {}", "Log file", &args.log.display());
    let duration = time.elapsed();
    println!();
    if duration.as_secs() < 60 {
        log::info!("{:18}: {:?}", "Execution time", duration);
    } else {
        let time = utils::parse_duration(duration.as_secs());
        log::info!("{:18}: {}", "Execution time (HH:MM:SS)", time);
    }
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
                        println!("{}", "Aborted!".red());
                        std::process::exit(0);
                    }
                }
                Err(err) => panic!("Failed parsing user input: {}", err),
            }
        }
    };
}

macro_rules! collect_paths {
    ($self: ident, $dir: ident, $input_fmt: ident) => {
        match &$dir {
            Some(dir) => {
                $self.input_dir = Some(PathBuf::from(&dir));
                $self.glob_paths(&dir, &$input_fmt)
            }
            None => $self.collect_paths(&$self.args.io.input),
        }
    };
}

pub(crate) use collect_paths;

trait InputCli {
    #[cfg(target_os = "windows")]
    fn collect_paths(&self, input: &Option<String>) -> Vec<PathBuf> {
        let inputs = input
            .iter()
            .map(|t| OsStr::new(t).to_string_lossy())
            .collect::<Vec<_>>();
        let files: Vec<PathBuf> = inputs
            .iter()
            .flat_map(|i| {
                glob(i)
                    .expect("Failed finding files")
                    .filter_map(|ok| ok.ok())
                    .collect::<Vec<PathBuf>>()
            })
            .collect();
        assert!(!files.is_empty(), "Empty folders!");
        files
    }

    #[cfg(not(target_os = "windows"))]
    fn collect_paths(&self, input: &Option<Vec<PathBuf>>) -> Vec<PathBuf> {
        match input {
            Some(paths) => paths.to_vec(),
            None => panic!("No input files!"),
        }
    }
}

trait AlignSeqInput {
    fn parse_input_fmt(&self, input_fmt: &str) -> InputFmt {
        match input_fmt {
            "auto" => InputFmt::Auto,
            "fasta" => InputFmt::Fasta,
            "nexus" => InputFmt::Nexus,
            "phylip" => InputFmt::Phylip,
            _ => unreachable!("Unknown input format. Supported format: auto, fasta, nexus, phylip"),
        }
    }

    fn parse_datatype(&self, datatype: &str) -> DataType {
        match datatype {
            "aa" => DataType::Aa,
            "dna" => DataType::Dna,
            "ignore" => DataType::Ignore,
            _ => unreachable!(),
        }
    }

    fn glob_paths(&self, dir: &str, input_fmt: &InputFmt) -> Vec<PathBuf> {
        Files::new(Path::new(dir)).find(input_fmt)
    }
}

trait InputPrint {
    fn print_input_info(&self) {}
}

impl InputPrint for RawReadPrint<'_> {
    fn print_input_info(&self) {
        if let Some(input) = self.input {
            log::info!("{:18}: {}", "Input dir", &input.display());
        } else {
            log::info!("{:18}: {}", "Input path", "STDIN");
        }
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&self.fcounts));
    }
}

trait RawInputCli {
    fn glob_paths(&self, dir: &str, input_fmt: &RawReadFmt) -> Vec<PathBuf> {
        Files::new(Path::new(dir)).find_raw_read(input_fmt)
    }
}

struct RawReadPrint<'a> {
    input: &'a Option<PathBuf>,
    input_fmt: &'a RawReadFmt,
    task_desc: &'a str,
    fcounts: usize,
}

impl<'a> RawReadPrint<'a> {
    fn new(
        input: &'a Option<PathBuf>,
        input_fmt: &'a RawReadFmt,
        task_desc: &'a str,
        fcounts: usize,
    ) -> Self {
        Self {
            input,
            input_fmt,
            task_desc,
            fcounts,
        }
    }

    fn print(&self) {
        self.print_input_info();
        log::info!("{:18}: {}\n", "Input format:", self.input_fmt);
        log::info!("{:18}: {}\n", "Task", self.task_desc);
    }
}

impl InputPrint for AlignSeqPrint<'_> {
    fn print_input_info(&self) {
        if let Some(input) = self.input {
            log::info!("{:18}: {}", "Input dir", &input.display());
        } else {
            log::info!("{:18}: {}", "Input path", "STDIN");
        }
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&self.fcounts));
    }
}

struct AlignSeqPrint<'a> {
    input: &'a Option<PathBuf>,
    input_fmt: &'a InputFmt,
    datatype: &'a DataType,
    task_desc: &'a str,
    fcounts: usize,
}

impl<'a> AlignSeqPrint<'a> {
    fn new(
        input: &'a Option<PathBuf>,
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        task_desc: &'a str,
        fcounts: usize,
    ) -> Self {
        Self {
            input,
            input_fmt,
            datatype,
            task_desc,
            fcounts,
        }
    }

    fn print(&self) {
        self.print_input_info();
        self.print_seq_input_fmt();
        self.print_seq_datatype();
        log::info!("{:18}: {}\n", "Task", self.task_desc);
    }

    fn print_seq_datatype(&self) {
        match self.datatype {
            DataType::Aa => log::info!("{:18}: {}", "Data type", "Amino Acid"),
            DataType::Dna => log::info!("{:18}: {}", "Data type", "DNA"),
            DataType::Ignore => log::info!("{:18}: {}", "Data type", "Ignore"),
        }
    }

    fn print_seq_input_fmt(&self) {
        match self.input_fmt {
            InputFmt::Auto => log::info!("{:18}: {}", "Input format", "Auto"),
            InputFmt::Fasta => log::info!("{:18}: {}", "Input format", "FASTA"),
            InputFmt::Nexus => log::info!("{:18}: {}", "Input format", "NEXUS"),
            InputFmt::Phylip => log::info!("{:18}: {}", "Input format", "PHYLIP"),
        }
    }
}

trait OutputCli {
    fn parse_output_fmt(&self, output_fmt: &str) -> OutputFmt {
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

    fn check_output_dir_exist(&self, path: &Path, overwrite: bool) {
        let rm_err_msg = "Failed removing a directory";
        if overwrite {
            if path.is_dir() {
                log::warn!(
                    "{} Removing existing directory: {}\n",
                    "WARNING!".red(),
                    path.display()
                );
                fs::remove_dir_all(path).expect(rm_err_msg);
            }
        } else {
            let error_msg = format!("Output dir already exists: {}. Remove?", path.display());
            check_output_path!(is_dir, remove_dir_all, path, error_msg, rm_err_msg)
        }
    }

    fn check_output_file_exist(&self, path: &Path, overwrite: bool) {
        let rm_err_msg = "Failed removing existing output files";
        if overwrite {
            if path.is_file() {
                log::warn!(
                    "{} Overwriting existing files: {}\n",
                    "WARNING!".red(),
                    path.display()
                );
                fs::remove_file(path).expect(rm_err_msg);
            }
        } else {
            let error_msg = format!("Output file already exists: {}. Remove?", path.display());
            check_output_path!(is_file, remove_file, path, error_msg, rm_err_msg);
        }
    }
}

trait ConcatCli {
    fn parse_prefix(&self, prefix: &Option<PathBuf>, output_dir: &Path) -> PathBuf {
        match prefix {
            Some(prefix) => prefix.to_path_buf(),
            None => output_dir.to_path_buf(),
        }
    }

    fn parse_partition_fmt(&self, part_fmt: &str, codon: bool) -> PartitionFmt {
        if codon {
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
