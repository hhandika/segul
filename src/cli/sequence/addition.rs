//! Parser for sequence Addition

use std::path::PathBuf;

use colored::Colorize;

use crate::{
    cli::{args::sequence::SequenceAddArgs, collect_paths, AlignSeqInput, InputCli, OutputCli},
    core::sequence::addition::SequenceAddition,
    helper::{
        logger::AlignSeqLogger,
        types::{InputFmt, OutputFmt},
    },
};

pub(in crate::cli) struct AdditionParser<'a> {
    args: &'a SequenceAddArgs,
    input_dir: Option<PathBuf>,
    dest_dir: Option<PathBuf>,
}

impl InputCli for AdditionParser<'_> {}
impl OutputCli for AdditionParser<'_> {
    fn parse_output_fmt(&self, output_fmt: &str) -> OutputFmt {
        match output_fmt {
            "fasta" => OutputFmt::Fasta,
            "fasta-int" => OutputFmt::FastaInt,
            _ => unreachable!("Output format is not supported"),
        }
    }
}
impl AlignSeqInput for AdditionParser<'_> {}

impl<'a> AdditionParser<'a> {
    pub(in crate::cli) fn new(args: &'a SequenceAddArgs) -> Self {
        Self {
            args,
            input_dir: None,
            dest_dir: None,
        }
    }

    pub(in crate::cli) fn add(&mut self) {
        let input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        let dest_fmt = self.parse_input_fmt(&self.args.destination_fmt);
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let output_fmt = self.parse_output_fmt(&self.args.output_fmt);
        let task = "Sequence addition";
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);
        let dest_files = self.collect_destination_paths(&dest_fmt);
        AlignSeqLogger::new(
            self.input_dir.as_deref(),
            &input_fmt,
            &datatype,
            files.len(),
        )
        .log(task);
        self.log_destination_info(&dest_files);
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        let add = SequenceAddition::new(
            &files,
            &input_fmt,
            &datatype,
            &self.args.output,
            &output_fmt,
        );
        add.add(&dest_files, &dest_fmt);
    }

    fn collect_destination_paths(&mut self, dest_fmt: &InputFmt) -> Vec<PathBuf> {
        match &self.args.destination_dir {
            Some(dir) => {
                self.dest_dir = Some(PathBuf::from(&dir));
                self.glob_paths(&dir, dest_fmt)
            }
            None => self.collect_paths(&self.args.destination_input),
        }
    }

    fn log_destination_info(&self, dest_files: &[PathBuf]) {
        log::info!("{}", "Destination Files".yellow());
        log::info!("{:18}: {}\n", "File counts", dest_files.len());
    }
}
