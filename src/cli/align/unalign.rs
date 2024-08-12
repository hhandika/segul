use std::path::PathBuf;

use crate::{
    cli::{args::align::UnalignArgs, collect_paths, AlignSeqInput, InputCli, OutputCli},
    core::align::unalign::UnalignAlignment,
    helper::logger::AlignSeqLogger,
};

pub(in crate::cli) struct UnalignParser<'a> {
    args: &'a UnalignArgs,
    input_dir: Option<PathBuf>,
}

impl InputCli for UnalignParser<'_> {}
impl OutputCli for UnalignParser<'_> {}
impl AlignSeqInput for UnalignParser<'_> {}

impl<'a> UnalignParser<'a> {
    pub(in crate::cli) fn new(args: &'a UnalignArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn unalign(&mut self) {
        let input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let output_fmt = self.parse_output_fmt(&self.args.output_fmt);
        let task = "Unalign sequence alignment";
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
        let unalign = UnalignAlignment::new(
            &files,
            &input_fmt,
            &datatype,
            &self.args.output,
            &output_fmt,
        );
        unalign.unalign();
    }
}
