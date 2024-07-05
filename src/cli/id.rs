use std::path::PathBuf;

use crate::core::sequence::id::SequenceID;
use crate::helper::logger::AlignSeqLogger;

use super::args::SequenceIdArgs;
use super::{collect_paths, AlignSeqInput, InputCli, OutputCli};

impl InputCli for IdParser<'_> {}
impl OutputCli for IdParser<'_> {}
impl AlignSeqInput for IdParser<'_> {}

pub(in crate::cli) struct IdParser<'a> {
    args: &'a SequenceIdArgs,
    input_dir: Option<PathBuf>,
}

impl<'a> IdParser<'a> {
    pub(in crate::cli) fn new(args: &'a SequenceIdArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn extract(&mut self) {
        let input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);
        let log = AlignSeqLogger::new(
            self.input_dir.as_deref(),
            &input_fmt,
            &datatype,
            files.len(),
        );
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        let id = SequenceID::new(
            &files,
            &input_fmt,
            &datatype,
            &self.args.output,
            self.args.prefix.as_deref(),
        );
        if self.args.map {
            let task = "Sequence ID Mapping";
            log.log(task);
            id.map_id();
        } else {
            let task = "Sequence ID Generation";
            log.log(task);
            id.get_unique();
        }
    }
}
