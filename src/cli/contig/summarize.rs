use std::path::PathBuf;

use segul_core::plot::contig::plot_contig_summary;
use segul_core::report::ReportContext;

use crate::cli::args::genomics::ContigSummaryArgs;
use crate::cli::plot::{PlotCli, append_plots_to_report};
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

    pub(in crate::cli) fn summarize(&mut self, ctx: &ReportContext) {
        let input_fmt = &self.args.input_format;
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);
        let fcounts = files.len();
        let task = "Summarize contig sequences";
        ContigLogger::new(self.input_dir.as_deref(), input_fmt, fcounts).log(task);
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        let prefix = self.args.prefix.as_deref();
        let written = ContigSummaryHandler::new(&files, input_fmt, &self.args.output, prefix)
            .with_report_context(ctx.clone().with_input_dir(self.input_dir.clone()))
            .summarize();
        let plots = self.write_plots(&self.args.plot, &self.args.output, prefix, |writer| {
            plot_contig_summary(&written.summary, writer)
        });
        append_plots_to_report(&written.report, &self.args.output, &plots);
    }
}
