use super::{args::RawSummaryArgs, InputCli, OutputCli};

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
        inputs.iter().for_each(|input| println!("{:?}", input));
    }
}
