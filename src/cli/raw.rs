use crate::{
    handler::raw::summarize::RawSummaryHandler,
    helper::types::{RawReadFmt, SummaryMode},
};

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
        RawReadPrint::new(&None, &input_fmt, task, fcounts).print();
        RawSummaryHandler::new(&inputs, &RawReadFmt::Auto, &SummaryMode::Default).summarize();
    }
}
