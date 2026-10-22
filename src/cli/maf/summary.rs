use segul_core::plot::maf::plot_maf_summary;

use crate::{core::maf::summary::MafSummary, helper::utils};

use super::MafInputCli;
use crate::cli::plot::PlotCli;
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

    pub(in crate::cli) fn summarize(&mut self) {
        let files = self.collect_maf_files(&self.args.io);
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&files.len()));
        log::info!("{:18}: {}", "Input format", "MAF");
        log::info!("{:18}: {}\n", "Task", "MAF summary");
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        let summary = MafSummary::new(&files, &self.args.output, self.args.prefix.as_deref())
            .with_naming(self.parse_naming(&self.args.common))
            .with_ref_species(self.args.common.ref_species.clone())
            .with_species(self.parse_species(&self.args.common));
        match summary.summarize() {
            Ok((per_file, per_species)) => {
                self.write_plots(
                    &self.args.plot,
                    &self.args.output,
                    self.args.prefix.as_deref(),
                    |writer| plot_maf_summary(&per_file, &per_species, writer),
                );
            }
            Err(e) => {
                log::error!("{:#}", e);
                std::process::exit(1);
            }
        }
    }
}
