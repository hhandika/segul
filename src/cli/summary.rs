use clap::ArgMatches;

use crate::cli::*;
use crate::helper::types::{DataType, InputFmt};

impl InputCli for SummaryParser<'_> {
    // We can't ignore datatype here because
    /// we generate different summary statistics for each data type.
    fn parse_datatype(&self, matches: &ArgMatches) -> DataType {
        let datatype = matches
            .value_of("datatype")
            .expect("Failed parsing dataype value");
        match datatype {
            "aa" => DataType::Aa,
            "dna" => DataType::Dna,
            _ => unreachable!("Please, define the data type of the file"),
        }
    }
}

impl InputPrint for SummaryParser<'_> {}

impl OutputCli for SummaryParser<'_> {
    fn parse_output<'a>(&self, matches: &'a ArgMatches) -> PathBuf {
        let output = matches.value_of("output").expect("CANNOT READ OUTPUT PATH");
        let csv = format!("{}_per_locus", output);
        PathBuf::from(csv).with_extension("csv")
    }
}

pub(in crate::cli) struct SummaryParser<'a> {
    matches: &'a ArgMatches<'a>,
    interval: usize,
    input_fmt: InputFmt,
    datatype: DataType,
}

impl<'a> SummaryParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            interval: 0,
            input_fmt: InputFmt::Fasta,
            datatype: DataType::Dna,
        }
    }

    pub(in crate::cli) fn stats(&mut self) {
        self.input_fmt = self.parse_input_fmt(self.matches);
        self.interval = self.parse_interval();
        self.datatype = self.parse_datatype(self.matches);
        let input_type = self.parse_input_type(self.matches);
        let task_desc = "Sequence summary statistics";
        match input_type {
            InputType::File => self.get_stats_file(task_desc),
            InputType::Dir => {
                let dir = self.parse_dir_input(self.matches);
                let files = self.get_files(&dir, &self.input_fmt);
                self.print_input_multi(&Some(dir), task_desc, files.len(), &self.input_fmt);
                self.get_stats_multiple(&files);
            }
            InputType::Wildcard => {
                let files = self.parse_input_wcard(self.matches);
                self.print_input_multi::<PathBuf>(&None, task_desc, files.len(), &self.input_fmt);
                self.get_stats_multiple(&files)
            }
        }
    }

    fn get_stats_multiple(&self, files: &[PathBuf]) {
        let output = self.parse_output(self.matches);
        SeqStats::new(&self.input_fmt, &output, self.interval, &self.datatype).get_stats_dir(files);
    }

    fn get_stats_file(&self, task_desc: &str) {
        self.parse_input_fmt(self.matches);
        let input = Path::new(self.parse_file_input(self.matches));
        let output = self.parse_output(self.matches);
        self.print_input_file(input, task_desc, &self.input_fmt);
        SeqStats::new(&self.input_fmt, &output, self.interval, &self.datatype)
            .get_seq_stats_file(input);
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
