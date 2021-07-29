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
            _ => unreachable!(),
        }
    }
}

impl OutputCli for SummaryParser<'_> {}

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
        self.input_fmt = self.parse_input_fmt(&self.matches);
        self.interval = self.parse_interval();
        self.datatype = self.parse_datatype(self.matches);
        let input_type = self.parse_input_type(&self.matches);
        match input_type {
            InputType::File => self.get_stats_file(),
            InputType::Dir => {
                let dir = self.parse_dir_input(self.matches);
                let files = self.get_files(dir, &self.input_fmt);
                self.print_input_multi(Some(dir), files.len());
                self.get_stats_multiple(&files);
            }
            InputType::Wildcard => {
                let files = self.parse_input_wcard(&self.matches);
                self.print_input_multi::<PathBuf>(None, files.len());
                self.get_stats_multiple(&files)
            }
        }
    }

    fn get_stats_multiple(&self, files: &[PathBuf]) {
        let output = self.parse_output(&self.matches);
        SeqStats::new(&self.input_fmt, output, self.interval, &self.datatype).get_stats_dir(&files);
    }

    fn get_stats_file(&self) {
        self.parse_input_fmt(&self.matches);
        let input = Path::new(self.parse_file_input(self.matches));
        let output = self.parse_output(&self.matches);
        self.print_input_file(&input);
        SeqStats::new(&self.input_fmt, output, self.interval, &self.datatype)
            .get_seq_stats_file(input);
    }

    fn parse_interval(&self) -> usize {
        let interval = self
            .matches
            .value_of("comp-interval")
            .expect("CAN'T GET INTERVAL VALUES");
        interval
            .parse::<usize>()
            .expect("FAIL PARSING INTERVAL VALUES")
    }

    fn print_input_file(&self, input: &Path) {
        log::info!("{:18}: {}", "Input", &input.display());
        self.print_input_fmt();
        log::info!("{:18}: Sequence summary statistics\n", "Task");
    }

    fn print_input_multi<P: AsRef<Path>>(&self, input: Option<P>, fcounts: usize) {
        if let Some(input) = input {
            log::info!("{:18}: {}", "Input dir", &input.as_ref().display());
        } else {
            log::info!("{:18}: {}", "Input dir", "WILDCARD");
        }
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&fcounts));
        self.print_input_fmt();
        log::info!("{:18}: Sequence summary statistics\n", "Task");
    }

    fn print_input_fmt(&self) {
        match self.input_fmt {
            InputFmt::Auto => log::info!("{:18}: {}", "Input format", "Auto"),
            InputFmt::Fasta => log::info!("{:18}: {}", "Input format", "Fasta"),
            InputFmt::Nexus => log::info!("{:18}: {}", "Input format", "Nexus"),
            InputFmt::Phylip => log::info!("{:18}: {}", "Input format", "Phylip"),
        }
    }
}
