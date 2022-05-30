use std::path::PathBuf;

use clap::ArgMatches;

use crate::cli::{InputCli, InputPrint, OutputCli};
use crate::handler::convert::Converter;
use crate::helper::types::{DataType, InputFmt, OutputFmt};

impl InputCli for ConvertParser<'_> {}
impl InputPrint for ConvertParser<'_> {}
impl OutputCli for ConvertParser<'_> {}

pub(in crate::cli) struct ConvertParser<'a> {
    matches: &'a ArgMatches,
    input_fmt: InputFmt,
    output: PathBuf,
    output_fmt: OutputFmt,
    datatype: DataType,
    sort: bool,
    input_dir: Option<PathBuf>,
}

impl<'a> ConvertParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches) -> Self {
        Self {
            matches,
            input_fmt: InputFmt::Auto,
            datatype: DataType::Dna,
            output: PathBuf::new(),
            output_fmt: OutputFmt::Nexus,
            sort: false,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn convert(&mut self) {
        self.input_fmt = self.parse_input_fmt(self.matches);
        self.output = self.parse_output(self.matches);
        self.output_fmt = self.parse_output_fmt(self.matches);
        self.datatype = self.parse_datatype(self.matches);
        self.is_sort();
        let task_desc = "Sequence format conversion";
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

        let is_overwrite = self.parse_overwrite_opts(self.matches);
        self.check_output_dir_exist(&self.output, is_overwrite);
        let mut convert =
            Converter::new(&self.input_fmt, &self.output_fmt, &self.datatype, self.sort);
        convert.convert(&files, &self.output);
    }

    fn is_sort(&mut self) {
        self.sort = self.matches.is_present("sort");
    }
}
