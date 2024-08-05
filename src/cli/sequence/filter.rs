use std::path::PathBuf;

use colored::Colorize;

use crate::{
    cli::{args::sequence::SequenceFilterArgs, collect_paths, AlignSeqInput, InputCli, OutputCli},
    core::sequence::filter::{SeqFilteringParameters, SequenceFiltering},
    helper::logger::AlignSeqLogger,
};

impl InputCli for SequenceFilterParser<'_> {}
impl OutputCli for SequenceFilterParser<'_> {}
impl AlignSeqInput for SequenceFilterParser<'_> {}

pub(in crate::cli) struct SequenceFilterParser<'a> {
    args: &'a SequenceFilterArgs,
    input_dir: Option<PathBuf>,
    params: SeqFilteringParameters,
}

impl<'a> SequenceFilterParser<'a> {
    pub(in crate::cli) fn new(args: &'a SequenceFilterArgs) -> Self {
        Self {
            args,
            input_dir: None,
            params: SeqFilteringParameters::None,
        }
    }

    pub(in crate::cli) fn filter(&mut self) {
        let input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let output_fmt = self.parse_output_fmt(&self.args.out_fmt.output_fmt);
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);
        let task = "Filter sequences";
        AlignSeqLogger::new(
            self.input_dir.as_deref(),
            &input_fmt,
            &datatype,
            files.len(),
        )
        .log(task);
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        log::info!("{}", "Filtering Parameters".yellow());
        self.parse_params();
        let filter = SequenceFiltering::new(
            &files,
            &input_fmt,
            &datatype,
            &self.args.output,
            &output_fmt,
            &self.params,
        );
        filter.filter();
    }

    fn parse_params(&mut self) {
        if let Some(min_len) = self.args.min_len {
            log::info!("{:18}: {}\n", "Minimum length", min_len);
            self.params = SeqFilteringParameters::MinSequenceLength(min_len);
        }
        if let Some(max_gap) = self.args.max_gap {
            let percent_max_gap = max_gap * 100.0;
            log::info!("{:18}: {}\n", "Max gaps", format!("{}%", percent_max_gap));
            self.params = SeqFilteringParameters::PercentMaxGap(percent_max_gap);
        }
    }
}
