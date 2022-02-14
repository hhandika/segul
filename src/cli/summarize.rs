use std::path::PathBuf;

use clap::ArgMatches;

use crate::cli::{InputCli, InputPrint, OutputCli};
use crate::core::summarize::SeqStats;
use crate::helper::types::{DataType, InputFmt};

impl InputCli for SummaryParser<'_> {}
impl InputPrint for SummaryParser<'_> {}
impl OutputCli for SummaryParser<'_> {}

pub(in crate::cli) struct SummaryParser<'a> {
    matches: &'a ArgMatches<'a>,
    interval: usize,
    input_fmt: InputFmt,
    datatype: DataType,
    input_dir: Option<PathBuf>,
}

impl<'a> SummaryParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            interval: 0,
            input_fmt: InputFmt::Fasta,
            datatype: DataType::Dna,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn stats(&mut self) {
        self.input_fmt = self.parse_input_fmt(self.matches);
        self.interval = self.parse_interval();
        self.datatype = self.parse_datatype(self.matches);
        let prefix = self.parse_prefix();
        let task_desc = "Sequence summary statistics";
        let files = if self.matches.is_present("dir") {
            let dir = self.parse_dir_input(self.matches);
            self.input_dir = Some(PathBuf::from(dir));
            self.get_files(dir, &self.input_fmt)
        } else {
            self.parse_input(self.matches)
        };

        self.print_input::<PathBuf>(
            &self.input_dir,
            task_desc,
            files.len(),
            &self.input_fmt,
            &self.datatype,
        );

        let output = self.parse_output(self.matches);
        self.check_output_dir_exist(&output);
        SeqStats::new(&self.input_fmt, &output, self.interval, &self.datatype)
            .get_stats_all(&files, &prefix);
    }

    fn parse_prefix(&self) -> Option<String> {
        if self.matches.is_present("prefix") {
            Some(
                self.matches
                    .value_of("prefix")
                    .expect("Failed parsing prefix input")
                    .to_string(),
            )
        } else {
            None
        }
    }

    fn parse_interval(&self) -> usize {
        let interval = self
            .matches
            .value_of("percent-interval")
            .expect("Failed parsing the interval command");
        interval
            .parse::<usize>()
            .expect("Failed parsing interval values to integer")
    }
}
