use std::path::PathBuf;

use super::{collect_paths, AlignSeqInput, ConcatCli, InputCli, OutputCli};
use crate::cli::args::AlignConcatArgs;
use crate::handler::align::concat::ConcatHandler;
use crate::helper::logger::AlignSeqLogger;

impl ConcatCli for ConcatParser<'_> {}
impl OutputCli for ConcatParser<'_> {}
impl InputCli for ConcatParser<'_> {}
impl AlignSeqInput for ConcatParser<'_> {}

pub(in crate::cli) struct ConcatParser<'a> {
    args: &'a AlignConcatArgs,
    input_dir: Option<PathBuf>,
}

impl<'a> ConcatParser<'a> {
    pub(in crate::cli) fn new(args: &'a AlignConcatArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn concat(&mut self) {
        let input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let output_fmt = self.parse_output_fmt(&self.args.out_fmt.output_fmt);
        let prefix = self.parse_prefix(&self.args.concat.prefix, &self.args.output);
        let part_fmt = self.parse_partition_fmt(&self.args.concat.part_fmt, self.args.concat.codon);
        self.check_partition_format(&output_fmt, &part_fmt);
        let task = "Alignment concatenation";
        let dir = &self.args.io.dir;
        let mut files = collect_paths!(self, dir, input_fmt);
        AlignSeqLogger::new(
            self.input_dir.as_deref(),
            &input_fmt,
            &datatype,
            files.len(),
        )
        .log(task);
        let is_overwrite = self.args.io.force;
        self.check_output_dir_exist(&self.args.output, is_overwrite);

        let mut concat = ConcatHandler::new(
            &input_fmt,
            &self.args.output,
            &output_fmt,
            &part_fmt,
            &prefix,
        );
        concat.concat_alignment(&mut files, &datatype);
    }
}
