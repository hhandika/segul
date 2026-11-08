use std::path::PathBuf;

use segul_core::plot::contig::plot_contig_summary;

use crate::cli::args::genomics::ContigSummaryArgs;
use crate::cli::plot::PlotCli;
use crate::cli::{ContigInputCli, InputCli, OutputCli, collect_paths};
use crate::core::contig::summarize::ContigSummaryHandler;
use crate::helper::logger::ContigLogger;

pub(in crate::cli) struct ContigCliParser<'a> {
    args: &'a ContigSummaryArgs,
    input_dir: Option<PathBuf>,
}

impl InputCli for ContigCliParser<'_> {}
impl OutputCli for ContigCliParser<'_> {}
impl ContigInputCli for ContigCliParser<'_> {}
impl PlotCli for ContigCliParser<'_> {}

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
        let prefix = self.args.prefix.as_deref();
        let summary =
            ContigSummaryHandler::new(&files, input_fmt, &self.args.output, prefix).summarize();
        self.write_plots(&self.args.plot, &self.args.output, prefix, |writer| {
            plot_contig_summary(&summary, writer)
        });
    }
}
