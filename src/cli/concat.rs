use std::path::PathBuf;

use super::{collect_paths, AlignSeqInput, AlignSeqPrint, ConcatCli, InputCli, OutputCli};
use crate::cli::args::AlignConcatArgs;
use crate::handler::align::concat::ConcatHandler;
use crate::helper::filenames;

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
        let output = filenames::create_output_fname(&self.args.output, &prefix, &output_fmt);
        let part_fmt = self.parse_partition_fmt(&self.args.concat.part_fmt, self.args.concat.codon);
        self.check_partition_format(&output_fmt, &part_fmt);
        let task_desc = "Alignment concatenation";
        let dir = &self.args.io.dir;
        let mut files = collect_paths!(self, dir, input_fmt);
        AlignSeqPrint::new(
            &self.input_dir,
            &input_fmt,
            &datatype,
            task_desc,
            files.len(),
        )
        .print();

        let is_overwrite = self.args.io.force;
        self.check_output_dir_exist(&self.args.output, is_overwrite);

        let mut concat = ConcatHandler::new(&input_fmt, &output, &output_fmt, &part_fmt);
        concat.concat_alignment(&mut files, &datatype);
    }
}
