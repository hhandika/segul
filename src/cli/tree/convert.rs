use std::path::{Path, PathBuf};

use segul_core::core::tree::convert::TreeConversion;
use segul_core::parser::tree::{ParseOptions, TreeFileFinder, TreeFmt};
use segul_core::report::ReportContext;
use segul_core::writer::tree::TreeWriteOptions;

use crate::cli::args::tree::TreeConvertArgs;
use crate::cli::{InputCli, OutputCli, collect_paths};
use crate::helper::utils;

pub(in crate::cli) struct TreeConvertParser<'a> {
    args: &'a TreeConvertArgs,
    input_dir: Option<PathBuf>,
}

impl InputCli for TreeConvertParser<'_> {}
impl OutputCli for TreeConvertParser<'_> {}

impl<'a> TreeConvertParser<'a> {
    pub(in crate::cli) fn new(args: &'a TreeConvertArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn convert(&mut self, ctx: &ReportContext) {
        let dir = &self.args.io.dir;
        let input_fmt = self.args.input_fmt;
        let files = collect_paths!(self, dir, input_fmt);

        match &self.input_dir {
            Some(input_dir) => log::info!("{:18}: {}", "Input dir", input_dir.display()),
            None => log::info!("{:18}: {}", "Input path", "STDIN"),
        }
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&files.len()));
        log::info!("{:18}: {}", "Input format", input_fmt);
        log::info!("{:18}: {}\n", "Task", "Tree format conversion");

        self.check_output_dir_exist(&self.args.output, self.args.io.force);

        let parse_opts = ParseOptions {
            parse_annotations: !self.args.strip_annotations,
            ..Default::default()
        };
        let write_opts = TreeWriteOptions {
            use_translate: !self.args.no_translate,
            write_annotations: !self.args.strip_annotations,
            gzip: self.args.gzip,
            ..Default::default()
        };

        let conversion = TreeConversion::new(&self.args.input_fmt, &self.args.output_fmt)
            .with_parse_options(parse_opts)
            .with_write_options(write_opts)
            .with_report_context(ctx.clone().with_input_dir(self.input_dir.clone()));

        if let Err(err) = conversion.convert_with_report(&files, &self.args.output) {
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
