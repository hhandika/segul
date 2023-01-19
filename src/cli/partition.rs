use std::ffi::OsStr;
use std::path::{Path, PathBuf};

use crate::cli::{ConcatCli, InputCli, InputPrint, OutputCli};
use crate::handler::partition::PartConverter;
use crate::helper::types::{DataType, PartitionFmt};
use crate::helper::utils;

use super::args::PartitionArgs;

impl InputPrint for PartParser<'_> {}
impl InputCli for PartParser<'_> {}
impl ConcatCli for PartParser<'_> {}
impl OutputCli for PartParser<'_> {}

pub(in crate::cli) struct PartParser<'a> {
    args: &'a PartitionArgs,
}

impl<'a> PartParser<'a> {
    pub(in crate::cli) fn new(args: &'a PartitionArgs) -> Self {
        Self { args }
    }

    pub(in crate::cli) fn convert(&self) {
        let inputs = self.collect_paths(&self.args.input);
        let input_counts = inputs.len();
        let in_part_fmt = if let Some(part_fmt) = &self.args.part_fmt {
            self.parse_partition_fmt_std(part_fmt)
        } else {
            PartitionFmt::Charset
        };
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let out_part_fmt = self.parse_partition_fmt(&self.args.out_part, self.args.codon);

        let task_desc = "Converting partitions";
        inputs.iter().for_each(|input| {
            self.print_input_info(input, task_desc, input_counts, &datatype);
            let output = self.construct_output_path(input, &out_part_fmt);
            self.check_output_file_exist(&output, self.args.force);
            let converter = PartConverter::new(input, &in_part_fmt, &output, &out_part_fmt);
            converter.convert(&datatype, self.args.skip_checking);
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

    fn print_input_info(&self, input: &Path, task_desc: &str, fcounts: usize, datatype: &DataType) {
        log::info!("{:18}: {}", "Input path", input.display());
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&fcounts));
        self.print_datatype(datatype);
        log::info!("{:18}: {}\n", "Task", task_desc);
    }
}
