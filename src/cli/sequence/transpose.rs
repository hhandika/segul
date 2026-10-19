use std::path::PathBuf;

use crate::core::sequence::transpose::SequenceTranspose;
use crate::helper::logger::AlignSeqLogger;

use crate::cli::args::sequence::SequenceTransposeArgs;
use crate::cli::{AlignSeqInput, InputCli, OutputCli, collect_paths};

impl InputCli for TransposeParser<'_> {}
impl OutputCli for TransposeParser<'_> {}
impl AlignSeqInput for TransposeParser<'_> {}

pub(in crate::cli) struct TransposeParser<'a> {
    args: &'a SequenceTransposeArgs,
    input_dir: Option<PathBuf>,
}

impl<'a> TransposeParser<'a> {
    pub(in crate::cli) fn new(args: &'a SequenceTransposeArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn transpose(&mut self) {
        let input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let output_fmt = self.parse_output_fmt(&self.args.output_fmt);
        let task = "Sequence Transposition";
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
        SequenceTranspose::new(&input_fmt, &datatype, &self.args.output, &output_fmt)
            .transpose(&files);
    }
}
