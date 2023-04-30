use std::path::PathBuf;

use crate::handler::raw::summarize::RawSummaryHandler;

use super::{args::RawSummaryArgs, collect_paths, InputCli, OutputCli, RawInputCli, RawReadPrint};

impl InputCli for RawSummaryParser<'_> {}
impl OutputCli for RawSummaryParser<'_> {}
impl RawInputCli for RawSummaryParser<'_> {}

pub(in crate::cli) struct RawSummaryParser<'a> {
    args: &'a RawSummaryArgs,
    input_dir: Option<PathBuf>,
}

impl<'a> RawSummaryParser<'a> {
    pub(in crate::cli) fn new(args: &'a RawSummaryArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn summarize(&mut self) {
        let input_fmt = &self.args.input_format;
        let dir = &self.args.io.dir;
        let mut files = collect_paths!(self, dir, input_fmt);
        let fcounts = files.len();
        let task = "Summarize raw read sequences";
        RawReadPrint::new(&self.input_dir, input_fmt, task, fcounts).print();
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        RawSummaryHandler::new(
            &mut files,
            &self.args.input_format,
            &self.args.mode,
            &self.args.output,
        )
        .summarize(self.args.low_mem);
    }
}
