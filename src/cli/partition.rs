use crate::handler::sequence::partition::PartConverter;
use crate::helper::partition::construct_output_path;
use crate::helper::types::PartitionFmt;
use crate::helper::{logger, utils};

use super::args::PartitionArgs;
use super::{AlignSeqInput, ConcatCli, InputCli, OutputCli};

impl InputCli for PartParser<'_> {}
impl ConcatCli for PartParser<'_> {}
impl OutputCli for PartParser<'_> {}
impl AlignSeqInput for PartParser<'_> {}

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

        inputs.iter().for_each(|input| {
            logger::log_input_partition(input, input_counts);
            let output = construct_output_path(input, &out_part_fmt);
            self.check_output_file_exist(&output, self.args.force);
            let converter = PartConverter::new(input, &in_part_fmt, &output, &out_part_fmt);
            converter.convert(&datatype, self.args.skip_checking);
            if input_counts > 1 {
                utils::print_divider();
            }
        });
    }
}
