use std::path::{Path, PathBuf};

use clap::ArgMatches;

use crate::cli::{ConcatCli, InputCli, InputPrint, OutputCli};
use crate::handler::split::Splitter;
use crate::helper::types::PartitionFmt;

impl InputPrint for SplitParser<'_> {}
impl OutputCli for SplitParser<'_> {}
impl InputCli for SplitParser<'_> {}
impl ConcatCli for SplitParser<'_> {
    fn parse_partition_fmt(&self, matches: &ArgMatches) -> PartitionFmt {
        let part_fmt = matches
            .value_of("partition")
            .expect("Failed parsing partition format");
        match part_fmt {
            "nexus" => PartitionFmt::Nexus,
            "raxml" => PartitionFmt::Raxml,
            _ => unreachable!("Failed parsing partition format"),
        }
    }
}

pub(in crate::cli) struct SplitParser<'a> {
    matches: &'a ArgMatches,
}

impl<'a> SplitParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches) -> Self {
        Self { matches }
    }

    pub(in crate::cli) fn split(&mut self) {
        let input = self.parse_input_path();
        let input_fmt = self.parse_input_fmt(self.matches);
        let datatype = self.parse_datatype(self.matches);
        let output_fmt = self.parse_output_fmt(self.matches);
        let output = self.parse_output(self.matches);

        // If users do not specify input partition.
        // Assume partition is in the sequence file.
        let partitions = if self.matches.is_present("input-partition") {
            self.parse_partition_path()
        } else {
            PathBuf::from(&input)
        };
        let part_fmt = self.parse_part_fmt(&partitions);
        let prefix = self.parse_prefix_input();
        let task_desc = "Alignment splitting";
        self.print_input(&None::<PathBuf>, task_desc, 1, &input_fmt, &datatype);
        let is_overwrite = self.parse_overwrite_opts(self.matches);
        let is_uncheck = self.parse_uncheck_part_flag(self.matches);
        self.check_output_dir_exist(&output, is_overwrite);
        let split = Splitter::new(&input, &datatype, &input_fmt, &output, &output_fmt);
        split.split_alignment(&partitions, &part_fmt, &prefix, is_uncheck);
    }

    fn parse_input_path(&self) -> PathBuf {
        let input_file = self
            .matches
            .value_of("input")
            .expect("Input file is required");
        PathBuf::from(input_file)
    }

    // Because this prefix is used for all output files
    // We have to parse it separately instead of using
    // the concat prefix that return PathBuf
    fn parse_prefix_input(&self) -> Option<String> {
        if self.matches.is_present("prefix") {
            let prefix = self
                .matches
                .value_of("prefix")
                .expect("Failed parsing prefix");
            Some(prefix.to_string())
        } else {
            None
        }
    }

    fn parse_partition_path(&self) -> PathBuf {
        PathBuf::from(
            self.matches
                .value_of("input-partition")
                .expect("No partition file provided"),
        )
    }

    fn parse_part_fmt(&self, part_path: &Path) -> PartitionFmt {
        if !self.matches.is_present("partition") {
            let ext = part_path
                .extension()
                .expect("Failed getting file extension")
                .to_str()
                .expect("Failed getting file extension as string");
            match ext {
                "txt" | "raxml" => PartitionFmt::Raxml,
                "nex" | "nexus" | "charset" => PartitionFmt::Nexus,
                _ => panic!(
                    "Cannot infer partition format from the file extension.\
                Please, specify using the --partition (or -p in short version) option"
                ),
            }
        } else {
            self.parse_partition_fmt(self.matches)
        }
    }
}
