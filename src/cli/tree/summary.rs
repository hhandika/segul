use std::path::{Path, PathBuf};

use segul_core::core::tree::estimate::RootingMethod;
use segul_core::core::tree::summary::{RootedMode, TreeSummarization, TreeSummaryOptions};
use segul_core::parser::tree::{TreeFileFinder, TreeFmt};
use segul_core::report::ReportContext;
use segul_core::writer::tree::TreeWriteOptions;

use crate::cli::args::tree::TreeSummaryArgs;
use crate::cli::{InputCli, OutputCli, collect_paths};
use crate::helper::utils;

pub(in crate::cli) struct TreeSummaryParser<'a> {
    args: &'a TreeSummaryArgs,
    input_dir: Option<PathBuf>,
}

impl InputCli for TreeSummaryParser<'_> {}
impl OutputCli for TreeSummaryParser<'_> {}

impl<'a> TreeSummaryParser<'a> {
    pub(in crate::cli) fn new(args: &'a TreeSummaryArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn summarize(&mut self, ctx: &ReportContext) {
        let dir = &self.args.io.dir;
        let input_fmt = self.args.input_fmt;
        let files = collect_paths!(self, dir, input_fmt);

        match &self.input_dir {
            Some(input_dir) => log::info!("{:18}: {}", "Input dir", input_dir.display()),
            None => log::info!("{:18}: {}", "Input path", "STDIN"),
        }
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&files.len()));
        log::info!("{:18}: {}", "Input format", input_fmt);
        log::info!("{:18}: {}", "Summary method", self.args.method);
        log::info!("{:18}: {}\n", "Task", "Phylogenetic tree summary");

        self.check_output_dir_exist(&self.args.output, self.args.io.force);

        let rooted_mode = if self.args.rooted {
            RootedMode::Rooted
        } else if self.args.unrooted {
            RootedMode::Unrooted
        } else {
            RootedMode::Auto
        };

        let rooting = if self.args.midpoint {
            RootingMethod::Midpoint
        } else if let Some(ref outgroup) = self.args.outgroup {
            RootingMethod::Outgroup(outgroup.clone())
        } else {
            RootingMethod::None
        };

        let options = TreeSummaryOptions {
            method: self.args.method,
            target_tree: self.args.target.clone(),
            burnin: self.args.burnin,
            min_clade_prob: self.args.min_clade_prob,
            rooted_mode,
            mcc_scoring: self.args.mcc_scoring,
            branch_lengths: self.args.branch_lengths,
            support_format: self.args.support_format,
            support_placement: self.args.support_placement,
            hpd_level: self.args.hpd,
            summarize_annotations: !self.args.no_annotations,
            allow_partial_taxa: self.args.allow_partial_taxa,
            rooting,
            output_fmt: self.args.output_fmt,
            write_opts: TreeWriteOptions {
                use_translate: !self.args.no_translate,
                write_annotations: !self.args.strip_annotations,
                gzip: self.args.gzip,
                ..Default::default()
            },
        };

        let summarization = TreeSummarization::new(&self.args.input_fmt, options)
            .with_report_context(ctx.clone().with_input_dir(self.input_dir.clone()));

        if let Err(err) = summarization.summarize_with_report(
            &files,
            &self.args.output,
            self.args.prefix.as_deref(),
        ) {
            log::error!("{:#}", err);
            std::process::exit(1);
        }
    }

    fn glob_paths(&self, dir: &str, input_fmt: &TreeFmt) -> Vec<PathBuf> {
        TreeFileFinder::new(Path::new(dir))
            .find(input_fmt)
            .unwrap_or_else(|err| {
                log::error!("{:#}", err);
                std::process::exit(1);
            })
    }
}
