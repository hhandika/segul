use std::path::{Path, PathBuf};

use clap::ArgMatches;

use crate::cli::{ConcatCli, InputCli, InputPrint, OutputCli};
use crate::core::split::Splitter;
use crate::helper::types::PartitionFmt;

impl InputPrint for SplitParser<'_> {}
impl OutputCli for SplitParser<'_> {
    fn parse_output<'a>(&self, matches: &'a ArgMatches) -> PathBuf {
        if !matches.is_present("output") {
            let output = matches.value_of("input").expect("Failed parsing input");
            let output_path = Path::new(output).file_stem().expect("Failed parsing input");
            PathBuf::from(output_path)
        } else {
            let output = matches
                .value_of("output")
                .expect("Failed parsing an output value");
            PathBuf::from(output)
        }
    }
}
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
    matches: &'a ArgMatches<'a>,
}

impl<'a> SplitParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self { matches }
    }

    pub(in crate::cli) fn split(&mut self) {
        let input = self.parse_input_path();
        let input_fmt = self.parse_input_fmt(self.matches);
        let datatype = self.parse_datatype(self.matches);
        let output_fmt = self.parse_output_fmt(self.matches);
        let output = self.parse_output(self.matches);
        let partitions = self.parse_partition_path();
        let part_fmt = self.parse_part_fmt(&partitions);
        let task_desc = "Alignment splitting";
        self.print_input(&None::<PathBuf>, task_desc, 1, &input_fmt, &datatype);
        self.check_output_dir_exist(&output);
        let split = Splitter::new(&input, &datatype, &input_fmt, &output, &output_fmt);
        split.split_alignment(&partitions, &part_fmt);
    }

    fn parse_input_path(&self) -> PathBuf {
        let input_file = self
            .matches
            .value_of("input")
            .expect("Input file is required");
        PathBuf::from(input_file)
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
                _ => panic!("Unsupported partition file format"),
            }
        } else {
            self.parse_partition_fmt(self.matches)
        }
    }
}
