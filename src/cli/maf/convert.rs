use std::path::PathBuf;

use crate::{
    cli::{args::genomics::MafConvertArgs, OutputCli},
    core::maf::convert::MafConverter,
};

pub(in crate::cli) struct MafConvertParser<'a> {
    args: &'a MafConvertArgs,
}

impl OutputCli for MafConvertParser<'_> {}

impl<'a> MafConvertParser<'a> {
    pub(in crate::cli) fn new(args: &'a MafConvertArgs) -> Self {
        Self { args }
    }

    pub(in crate::cli) fn convert(&mut self) {
        // let task = "Sequence format conversion";
        let input: Vec<PathBuf> = vec![self.args.input.clone()];
        let convert = MafConverter::new(
            &input,
            &self.args.reference_path,
            self.args.from_bed,
            &self.args.output,
            &self.args.output_fmt,
        );
        convert.convert();
    }
}
