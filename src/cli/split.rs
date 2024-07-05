use std::path::Path;

use crate::core::align::split::AlignmentSplitting;
use crate::helper::logger::AlignSeqLogger;
use crate::helper::types::PartitionFmt;

use super::args::AlignSplitArgs;
use super::{AlignSeqInput, ConcatCli, InputCli, OutputCli};

impl OutputCli for SplitParser<'_> {}
impl InputCli for SplitParser<'_> {}
impl ConcatCli for SplitParser<'_> {}
impl AlignSeqInput for SplitParser<'_> {}

pub(in crate::cli) struct SplitParser<'a> {
    args: &'a AlignSplitArgs,
}

impl<'a> SplitParser<'a> {
    pub(in crate::cli) fn new(args: &'a AlignSplitArgs) -> Self {
        Self { args }
    }

    pub(in crate::cli) fn split(&mut self) {
        let input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let output_fmt = self.parse_output_fmt(&self.args.out_fmt.output_fmt);

        // If users do not specify input partition.
        // Assume partition is in the sequence file.
        let partitions = match &&self.args.input_partition {
            Some(path) => path,
            None => &self.args.input,
        };

        let part_fmt = self.parse_part_fmt(partitions);
        let task = "Alignment splitting";
        AlignSeqLogger::new(None, &input_fmt, &datatype, 1).log(task);
        self.check_output_dir_exist(&self.args.output, self.args.force);
        let split = AlignmentSplitting::new(
            &self.args.input,
            &datatype,
            &input_fmt,
            &self.args.output,
            &output_fmt,
        );
        split.split(
            partitions,
            &part_fmt,
            &self.args.prefix,
            self.args.skip_checking,
        );
    }

    fn parse_part_fmt(&self, part_path: &Path) -> PartitionFmt {
        match &self.args.part_fmt {
            Some(fmt) => self.parse_partition_fmt(fmt),
            None => self.infer_part_fmt(part_path),
        }
    }

    fn infer_part_fmt(&self, part_path: &Path) -> PartitionFmt {
        let ext = part_path
            .extension()
            .expect("Failed getting file extension")
            .to_str()
            .expect("Failed getting file extension as string");
        self.parse_partition_ext(ext)
    }

    fn parse_partition_fmt(&self, fmt: &str) -> PartitionFmt {
        match fmt {
            "raxml" => PartitionFmt::Raxml,
            "nexus" => PartitionFmt::Nexus,
            _ => unreachable!(
                "Cannot infer partition format from the file extension.\
                Please, specify using the --partition (or -p in short version) option"
            ),
        }
    }

    fn parse_partition_ext(&self, ext: &str) -> PartitionFmt {
        match ext {
            "txt" | "raxml" => PartitionFmt::Raxml,
            "nex" | "nexus" | "charset" => PartitionFmt::Nexus,
            _ => panic!(
                "Cannot infer partition format from the file extension.\
                Please, specify using the --partition (or -p in short version) option"
            ),
        }
    }
}
