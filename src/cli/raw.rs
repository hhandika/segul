use crate::handler::raw::summarize::RawSummaryHandler;

use super::{args::RawSummaryArgs, InputCli, OutputCli, RawReadPrint};

impl InputCli for RawSummaryParser<'_> {}
impl OutputCli for RawSummaryParser<'_> {}

pub(in crate::cli) struct RawSummaryParser<'a> {
    args: &'a RawSummaryArgs,
}

impl<'a> RawSummaryParser<'a> {
    pub(in crate::cli) fn new(args: &'a RawSummaryArgs) -> Self {
        Self { args }
    }

    pub(in crate::cli) fn summarize(&self) {
        let inputs = self.collect_paths(&self.args.io.input);
        let fcounts = inputs.len();
        let input_fmt = &self.args.input_format;
        let task = "Summarize raw read sequences";
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        RawReadPrint::new(&None, input_fmt, task, fcounts).print();
        RawSummaryHandler::new(
            &inputs,
            &self.args.input_format,
            &self.args.mode,
            &self.args.output,
        )
        .summarize();
    }
}
