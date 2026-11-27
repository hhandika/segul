use segul_core::plot::maf::plot_maf_summary;
use segul_core::report::ReportContext;

use crate::{core::maf::summary::MafSummary, helper::utils};

use super::MafInputCli;
use crate::cli::plot::{PlotCli, append_plots_to_report};
use crate::cli::{InputCli, OutputCli, args::genomics::MafSummaryArgs};

pub(in crate::cli) struct MafSummaryParser<'a> {
    args: &'a MafSummaryArgs,
}

impl InputCli for MafSummaryParser<'_> {}
impl OutputCli for MafSummaryParser<'_> {}
impl MafInputCli for MafSummaryParser<'_> {}
impl PlotCli for MafSummaryParser<'_> {}

impl<'a> MafSummaryParser<'a> {
    pub(in crate::cli) fn new(args: &'a MafSummaryArgs) -> Self {
        Self { args }
    }

    pub(in crate::cli) fn summarize(&mut self, ctx: &ReportContext) {
        let files = self.collect_maf_files(&self.args.io);
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&files.len()));
        log::info!("{:18}: {}", "Input format", "MAF");
        log::info!("{:18}: {}\n", "Task", "MAF summary");
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        let input_dir = self.args.io.dir.as_ref().map(std::path::PathBuf::from);
        let summary = MafSummary::new(&files, &self.args.output, self.args.prefix.as_deref())
            .with_naming(self.parse_naming(&self.args.common))
            .with_ref_species(self.args.common.ref_species.clone())
            .with_species(self.parse_species(&self.args.common))
            .with_report_context(ctx.clone().with_input_dir(input_dir));
        match summary.try_summarize_with_report() {
            Ok(written) => {
                let plots = self.write_plots(
                    &self.args.plot,
                    &self.args.output,
                    self.args.prefix.as_deref(),
                    |writer| plot_maf_summary(&written.per_file, &written.per_species, writer),
                );
                append_plots_to_report(&written.report, &self.args.output, &plots);
            }
            Err(e) => {
                log::error!("{:#}", e);
                std::process::exit(1);
            }
        }
    }
}
