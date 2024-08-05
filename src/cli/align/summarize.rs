use std::path::PathBuf;

use crate::core::align::summarize::AlignmentSummary;
use crate::helper::logger::AlignSeqLogger;

use crate::cli::args::AlignSummaryArgs;
use crate::cli::{collect_paths, AlignSeqInput, InputCli, OutputCli};

impl InputCli for SummaryParser<'_> {}
impl OutputCli for SummaryParser<'_> {}
impl AlignSeqInput for SummaryParser<'_> {}

pub(in crate::cli) struct SummaryParser<'a> {
    args: &'a AlignSummaryArgs,
    input_dir: Option<PathBuf>,
}

impl<'a> SummaryParser<'a> {
    pub(in crate::cli) fn new(args: &'a AlignSummaryArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn summarize(&mut self) {
        let input_fmt = self.parse_input_fmt(&self.args.fmt.input_fmt);
        let datatype = self.parse_datatype(&self.args.fmt.datatype);
        let task = "Sequence summary statistics";
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
        let mut summary =
            AlignmentSummary::new(&input_fmt, &self.args.output, self.args.interval, &datatype);
        if self.args.per_locus {
            summary.summarize_locus(&files, self.args.prefix.as_deref());
        } else {
            summary.summarize_all(&files, self.args.prefix.as_deref());
        }
    }
}
