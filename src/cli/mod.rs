mod args;
mod concat;
mod convert;
mod filter;
mod id;
mod summary;

use std::io::{self, BufWriter, Result, Write};
use std::path::{Path, PathBuf};

use clap::ArgMatches;
use indexmap::IndexSet;
use rayon::prelude::*;

use crate::cli::concat::ConcatParser;
use crate::cli::convert::ConvertParser;
use crate::cli::filter::FilterParser;
use crate::cli::id::IdParser;
use crate::cli::summary::SummaryParser;

use crate::core::msa;
use crate::core::summary::SeqStats;
use crate::helper::common::{DataType, InputFmt, OutputFmt, PartitionFmt};
use crate::helper::finder::{Files, IDs};
use crate::helper::utils;

pub fn parse_cli(version: &str) {
    let args = args::get_args(version);
    let text = format!("SEGUL v{}", version);
    utils::print_title(&text);
    match args.subcommand() {
        ("convert", Some(convert_matches)) => ConvertParser::new(convert_matches).convert(),
        ("concat", Some(concat_matches)) => ConcatParser::new(concat_matches).concat(),
        ("filter", Some(pick_matches)) => FilterParser::new(pick_matches).filter(),
        ("id", Some(id_matches)) => IdParser::new(id_matches).get_id(),
        ("summary", Some(stats_matches)) => SummaryParser::new(stats_matches).stats(),
        _ => unreachable!(),
    }
}

enum InputType {
    File,
    Dir,
    Wildcard,
}

trait Cli {
    fn get_file_input<'a>(&self, matches: &'a ArgMatches) -> &'a str {
        matches
            .value_of("input")
            .expect("CANNOT FIND AN INPUT FILE")
    }

    fn get_dir_input<'a>(&self, matches: &'a ArgMatches) -> &'a str {
        matches.value_of("dir").expect("CANNOT READ DIR PATH")
    }

    fn parse_input_wcard(&self, matches: &ArgMatches) -> Vec<PathBuf> {
        matches
            .values_of("wildcard")
            .expect("FAILED PARSING npercent")
            .map(PathBuf::from)
            .collect()
    }

    fn get_files(&self, dir: &str, input_fmt: &InputFmt) -> Vec<PathBuf> {
        Files::new(dir, input_fmt).get_files()
    }

    fn get_input_type(&self, matches: &ArgMatches) -> InputType {
        if matches.is_present("input") {
            InputType::File
        } else if matches.is_present("dir") {
            InputType::Dir
        } else {
            InputType::Wildcard
        }
    }

    fn get_output<'a>(&self, matches: &'a ArgMatches) -> &'a str {
        matches.value_of("output").expect("CANNOT READ OUTPUT PATH")
    }

    fn get_output_path(&self, matches: &ArgMatches) -> PathBuf {
        if matches.is_present("output") {
            let output = self.get_output(matches);
            PathBuf::from(output)
        } else {
            PathBuf::from(".")
        }
    }

    fn get_input_fmt(&self, matches: &ArgMatches) -> InputFmt {
        let input_fmt = matches
            .value_of("format")
            .expect("CANNOT READ FORMAT INPUT");
        match input_fmt {
            "auto" => InputFmt::Auto,
            "fasta" => InputFmt::Fasta,
            "nexus" => InputFmt::Nexus,
            "phylip" => InputFmt::Phylip,
            _ => panic!(
                "UNSUPPORTED FORMAT. \
        THE PROGRAM ONLY ACCEPT fasta, fasta-int, nexus, nexus-int, and phylip. \
        ALL IN lowercase. \
        YOUR INPUT: {} ",
                input_fmt
            ),
        }
    }

    fn get_output_fmt(&self, matches: &ArgMatches) -> OutputFmt {
        let output_fmt = matches
            .value_of("output-format")
            .expect("CANNOT READ FORMAT INPUT");
        match output_fmt {
            "nexus" => OutputFmt::Nexus,
            "phylip" => OutputFmt::Phylip,
            "fasta" => OutputFmt::Fasta,
            "nexus-int" => OutputFmt::NexusInt,
            "fasta-int" => OutputFmt::FastaInt,
            "phylip-int" => OutputFmt::PhylipInt,
            _ => panic!(
                "UNSUPPORTED FORMAT. \
        THE PROGRAM ONLY ACCEPT fasta, fasta-int, nexus, nexus-int, phylip, and phylip-int. ALL IN lowercase. \
        YOUR INPUT: {} ",
                output_fmt
            ),
        }
    }

    fn get_datatype(&self, matches: &ArgMatches) -> DataType {
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
}

trait PartCLi {
    fn parse_partition_fmt(&self, matches: &ArgMatches) -> PartitionFmt {
        let part_fmt = matches
            .value_of("partition")
            .expect("CANNOT READ PARTITION FORMAT");
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
                if let PartitionFmt::Nexus | PartitionFmt::NexusCodon = part_fmt {
                    panic!(
                        "CANNOT WRITE EMBEDDED-NEXUS PARTITION TO NON-NEXUS OUTPUT. \
                MAYBE YOU MEAN TO WRITE THE PARTITION TO 'charset' INSTEAD."
                    )
                }
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use clap::{App, Arg};

    #[test]
    fn get_id_output_path_test() {
        let arg = App::new("segul-test")
            .arg(Arg::with_name("dir").default_value("./test_dir/"))
            .get_matches();
        let id = IdParser::new(&arg);
        let res = PathBuf::from("./test_dir.txt");
        assert_eq!(res, id.get_output_path(&arg));
    }
}
