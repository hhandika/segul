use clap::ArgMatches;

use crate::cli::*;
use crate::helper::common::{DataType, InputFmt};

impl InputCli for SummaryParser<'_> {}
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
        self.print_input_file().unwrap();
        let input_type = self.parse_input_type(&self.matches);
        match input_type {
            InputType::File => self.get_stats_file(),
            InputType::Dir => {
                let dir = self.parse_dir_input(self.matches);
                let files = self.get_files(dir, &self.input_fmt);
                self.get_stats_multiple(&files);
            }
            InputType::Wildcard => {
                let files = self.parse_input_wcard(&self.matches);
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
        self.print_input_file().unwrap();
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

    fn is_input_dir(&self) -> bool {
        self.matches.is_present("dir")
    }

    // FIXME: This does not cover all cases
    fn print_input_file(&self) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "Command\t\t: segul summary")?;
        if self.is_input_dir() {
            writeln!(
                writer,
                "Input dir\t: {}\n",
                self.parse_dir_input(self.matches)
            )?;
        } else {
            writeln!(writer, "Input\t\t: WILDCARD\n",)?;
        }

        Ok(())
    }
}
