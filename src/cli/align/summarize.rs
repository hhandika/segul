use std::path::PathBuf;

use crate::core::align::summarize::AlignmentSummary;
use crate::helper::logger::AlignSeqLogger;

use segul_core::plot::align::plot_alignment_summary;
use segul_core::report::ReportContext;

use crate::cli::args::align::AlignSummaryArgs;
use crate::cli::plot::{PlotCli, append_plots_to_report};
use crate::cli::{AlignSeqInput, InputCli, OutputCli, collect_paths};

impl InputCli for SummaryParser<'_> {}
impl OutputCli for SummaryParser<'_> {}
impl AlignSeqInput for SummaryParser<'_> {}
impl PlotCli for SummaryParser<'_> {}

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

    pub(in crate::cli) fn summarize(&mut self, ctx: &ReportContext) {
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
            AlignmentSummary::new(&input_fmt, &self.args.output, self.args.interval, &datatype)
                .with_report_context(ctx.clone().with_input_dir(self.input_dir.clone()));
        let prefix = self.args.prefix.as_deref();
        if self.args.per_locus {
            summary
                .summarize_locus_with_report(&files, prefix)
                .unwrap_or_else(|err| panic!("{err}"));
        } else {
            let outputs = summary.summarize_all(&files, prefix);
            let plots = self.write_plots(&self.args.plot, &self.args.output, prefix, |writer| {
                plot_alignment_summary(&outputs, writer)
            });
            append_plots_to_report(&outputs.report, &self.args.output, &plots);
        }
    }
}
