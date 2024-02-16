use std::path::PathBuf;

use super::{args::ContigSummaryArgs, collect_paths, ContigInputCli, InputCli, OutputCli};
use crate::handler::contig::summarize::ContigSummaryHandler;
use crate::helper::logger::ContigLogger;

pub(in crate::cli) struct ContigCliParser<'a> {
    args: &'a ContigSummaryArgs,
    input_dir: Option<PathBuf>,
}

impl InputCli for ContigCliParser<'_> {}
impl OutputCli for ContigCliParser<'_> {}
impl ContigInputCli for ContigCliParser<'_> {}

impl<'a> ContigCliParser<'a> {
    pub(in crate::cli) fn new(args: &'a ContigSummaryArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn summarize(&mut self) {
        let input_fmt = &self.args.input_format;
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);
        let fcounts = files.len();
        let task = "Summarize contig sequences";
        ContigLogger::new(self.input_dir.as_deref(), input_fmt, fcounts).log(task);
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        ContigSummaryHandler::new(
            &files,
            input_fmt,
            &self.args.output,
            self.args.prefix.as_deref(),
        )
        .summarize();
    }
}
