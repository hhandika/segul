use std::path::PathBuf;

use crate::{core::read::summarize::GenomicReadSummary, helper::logger::ReadLogger};

use super::{args::SeqReadSummaryArgs, collect_paths, InputCli, OutputCli, RawInputCli};

impl InputCli for ReadSummaryCliParser<'_> {}
impl OutputCli for ReadSummaryCliParser<'_> {}
impl RawInputCli for ReadSummaryCliParser<'_> {}

pub(in crate::cli) struct ReadSummaryCliParser<'a> {
    args: &'a SeqReadSummaryArgs,
    input_dir: Option<PathBuf>,
}

impl<'a> ReadSummaryCliParser<'a> {
    pub(in crate::cli) fn new(args: &'a SeqReadSummaryArgs) -> Self {
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
        ReadLogger::new(self.input_dir.as_deref(), input_fmt, fcounts).log(task);
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        GenomicReadSummary::new(
            &mut files,
            &self.args.input_format,
            &self.args.mode,
            &self.args.output,
            self.args.prefix.as_deref(),
        )
        .summarize();
    }
}
