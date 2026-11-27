use std::path::PathBuf;

use segul_core::report::ReportContext;

use crate::{core::align::convert::AlignmentConversion, helper::logger::AlignSeqLogger};

use crate::cli::args::align::AlignConvertArgs;

use crate::cli::{AlignSeqInput, InputCli, OutputCli, collect_paths};

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

    pub(in crate::cli) fn convert(&mut self, ctx: &ReportContext) {
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
        let convert = AlignmentConversion::new(&input_fmt, &output_fmt, &datatype, self.args.sort)
            .with_report_context(ctx.clone().with_input_dir(self.input_dir.clone()));
        convert
            .convert_with_report(&files, &self.args.output)
            .unwrap_or_else(|err| panic!("{err}"));
    }
}
