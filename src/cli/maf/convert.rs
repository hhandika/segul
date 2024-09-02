use std::path::Path;

use crate::{
    cli::{args::genomics::GenomicConvertArgs, InputCli, OutputCli},
    core::maf::convert::MafConverter,
    helper::{finder::MafFileFinder, utils},
};

pub(in crate::cli) struct MafConvertParser<'a> {
    args: &'a GenomicConvertArgs,
}

impl InputCli for MafConvertParser<'_> {}

impl OutputCli for MafConvertParser<'_> {}

impl<'a> MafConvertParser<'a> {
    pub(in crate::cli) fn new(args: &'a GenomicConvertArgs) -> Self {
        Self { args }
    }

    pub(in crate::cli) fn convert(&mut self) {
        let output_fmt = self.parse_output_fmt(&self.args.output_fmt);
        let files = match &self.args.io.dir {
            Some(dir) => {
                log::info!("{:18}: {}", "Input dir", &dir);
                MafFileFinder::new(Path::new(dir)).find_recursive()
            }
            None => {
                log::info!("{:18}: {}", "Input path", "STDIN");
                self.collect_paths(&self.args.io.input)
            }
        };
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&files.len()));
        let task = "MAF format conversion";
        log::info!("{:18}: {}", "Input format:", "MAF");
        log::info!("{:18}: {}\n", "Task", task);
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        let convert = MafConverter::new(
            &files,
            &self.args.reference,
            self.args.from_bed,
            &self.args.output,
            &output_fmt,
        );
        convert.convert();
    }
}
