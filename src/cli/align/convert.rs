use std::path::PathBuf;

use crate::{core::align::convert::AlignmentConversion, helper::logger::AlignSeqLogger};

use crate::cli::args::align::AlignConvertArgs;

use crate::cli::{collect_paths, AlignSeqInput, InputCli, OutputCli};

impl InputCli for ConvertParser<'_> {}
impl OutputCli for ConvertParser<'_> {}
impl AlignSeqInput for ConvertParser<'_> {}

pub(in crate::cli) struct ConvertParser<'a> {
    args: &'a AlignConvertArgs,
    input_dir: Option<PathBuf>,
}

impl<'a> ConvertParser<'a> {
    pub(in crate::cli) fn new(args: &'a AlignConvertArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn convert(&mut self) {
        let input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let output_fmt = self.parse_output_fmt(&self.args.out_fmt.output_fmt);
        let task = "Sequence format conversion";
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);
        AlignSeqLogger::new(
            self.input_dir.as_deref(),
            &input_fmt,
            &datatype,
            files.len(),
        )
        .log(task);
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        let convert = AlignmentConversion::new(&input_fmt, &output_fmt, &datatype, self.args.sort);
        convert.convert(&files, &self.args.output);
    }
}
