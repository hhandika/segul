use std::path::PathBuf;

use segul_core::report::ReportContext;

use crate::cli::args::align::AlignConcatArgs;
use crate::cli::{AlignSeqInput, ConcatCli, InputCli, OutputCli, collect_paths};
use crate::core::align::concat::AlignmentConcatenation;
use crate::helper::logger::AlignSeqLogger;

impl ConcatCli for ConcatParser<'_> {}
impl OutputCli for ConcatParser<'_> {}
impl InputCli for ConcatParser<'_> {}
impl AlignSeqInput for ConcatParser<'_> {}

pub(in crate::cli) struct ConcatParser<'a> {
    args: &'a AlignConcatArgs,
    input_dir: Option<PathBuf>,
}

impl<'a> ConcatParser<'a> {
    pub(in crate::cli) fn new(args: &'a AlignConcatArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn concat(&mut self, ctx: &ReportContext) {
        let input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let output_fmt = self.parse_output_fmt(&self.args.out_fmt.output_fmt);
        let prefix = self.parse_prefix(&self.args.concat.prefix, &self.args.output);
        let part_fmt = self.parse_partition_fmt(&self.args.concat.part_fmt, self.args.concat.codon);
        self.check_partition_format(&output_fmt, &part_fmt);
        let task = "Alignment concatenation";
        let dir = &self.args.io.dir;
        let mut files = collect_paths!(self, dir, input_fmt);
        AlignSeqLogger::new(
            self.input_dir.as_deref(),
            &input_fmt,
            &datatype,
            files.len(),
        )
        .log(task);
        let is_overwrite = self.args.io.force;
        self.check_output_dir_exist(&self.args.output, is_overwrite);

        let prefix_name = self
            .args
            .concat
            .prefix
            .as_deref()
            .and_then(|p| p.file_name())
            .and_then(|s| s.to_str());
        let mut concat = AlignmentConcatenation::new(
            &input_fmt,
            &self.args.output,
            &output_fmt,
            &part_fmt,
            &prefix,
        )
        .with_sort(self.args.sort)
        .with_report_context(ctx.clone().with_input_dir(self.input_dir.clone()));
        concat
            .concat_with_report(&mut files, &datatype, prefix_name)
            .unwrap_or_else(|err| panic!("{err}"));
    }
}
