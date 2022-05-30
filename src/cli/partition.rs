use std::ffi::OsStr;
use std::path::{Path, PathBuf};

use clap::ArgMatches;

use crate::cli::{ConcatCli, InputCli, InputPrint, OutputCli};
use crate::handler::partition::PartConverter;
use crate::helper::types::{DataType, PartitionFmt};
use crate::helper::utils;

impl InputPrint for PartParser<'_> {}
impl InputCli for PartParser<'_> {}
impl ConcatCli for PartParser<'_> {}
impl OutputCli for PartParser<'_> {}

pub(in crate::cli) struct PartParser<'a> {
    matches: &'a ArgMatches,
}

impl<'a> PartParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches) -> Self {
        Self { matches }
    }

    pub(in crate::cli) fn convert(&self) {
        let inputs = self.parse_input(self.matches);
        let input_counts = inputs.len();
        let in_part_fmt = if self.matches.is_present("partition") {
            self.parse_partition_fmt(self.matches)
        } else {
            PartitionFmt::Charset
        };
        let datatype = self.parse_datatype(self.matches);
        let out_part_fmt = self.parse_out_part_fmt();
        let is_overwrite = self.parse_overwrite_opts(self.matches);
        let is_uncheck = self.parse_uncheck_part_flag(self.matches);
        let task_desc = "Converting partitions";
        inputs.iter().for_each(|input| {
            self.print_input_info(input, task_desc, input_counts, &datatype);
            let output = self.construct_output_path(input, &out_part_fmt);
            self.check_output_file_exist(&output, is_overwrite);
            let converter = PartConverter::new(input, &in_part_fmt, &output, &out_part_fmt);
            converter.convert(&datatype, is_uncheck);
            if input_counts > 1 {
                utils::print_divider();
            }
        });
    }

    fn construct_output_path(&self, input: &Path, out_part_fmt: &PartitionFmt) -> PathBuf {
        let fstem = input
            .file_stem()
            .and_then(OsStr::to_str)
            .expect("Failed to parse input file stem");
        let mut fname = PathBuf::from(format!("{}_partition", fstem));
        match *out_part_fmt {
            PartitionFmt::Nexus | PartitionFmt::NexusCodon => {
                fname.set_extension("nex");
            }
            PartitionFmt::Raxml | PartitionFmt::RaxmlCodon => {
                fname.set_extension("txt");
            }
            _ => unreachable!("Failed to parse partition format"),
        }
        let parent_path = input.parent().expect("Failed to parse input parent path");
        parent_path.join(fname)
    }

    fn parse_out_part_fmt(&self) -> PartitionFmt {
        let part_fmt = self
            .matches
            .value_of("output-partition")
            .expect("Failed parsing partition format");
        if self.matches.is_present("codon") {
            self.parse_partition_fmt_codon(part_fmt)
        } else {
            self.parse_partition_fmt_std(part_fmt)
        }
    }

    fn print_input_info(&self, input: &Path, task_desc: &str, fcounts: usize, datatype: &DataType) {
        log::info!("{:18}: {}", "Input path", input.display());
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&fcounts));
        self.print_datatype(datatype);
        log::info!("{:18}: {}\n", "Task", task_desc);
    }
}
