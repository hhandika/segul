use std::path::PathBuf;

use clap::ArgMatches;

use crate::cli::{InputCli, InputPrint, OutputCli};
use crate::core::id::Id;

impl InputCli for IdParser<'_> {}
impl OutputCli for IdParser<'_> {}
impl InputPrint for IdParser<'_> {}

pub(in crate::cli) struct IdParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_dir: Option<PathBuf>,
}

impl<'a> IdParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches) -> Self {
        Self {
            matches,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn get_id(&mut self) {
        let input_fmt = self.parse_input_fmt(self.matches);
        let datatype = self.parse_datatype(self.matches);
        let task_desc = "IDs finding";
        let mut files = if self.is_input_dir() {
            let dir = self.parse_dir_input(self.matches);
            self.input_dir = Some(PathBuf::from(dir));
            self.get_files(dir, &input_fmt)
        } else {
            self.parse_input_wcard(self.matches)
        };

        self.print_input_multi(
            &self.input_dir,
            task_desc,
            files.len(),
            &input_fmt,
            &datatype,
        );
        let output = self.parse_output(self.matches);
        if self.matches.is_present("map") {
            let output = output.with_extension("csv");
            self.check_output_file_exist(&output);
            let id = Id::new(&output, &input_fmt, &datatype);
            alphanumeric_sort::sort_path_slice(&mut files);
            id.map_id(&files);
        } else {
            let output = output.with_extension("txt");
            self.check_output_file_exist(&output);
            let id = Id::new(&output, &input_fmt, &datatype);
            id.generate_id(&files);
        }
    }

    fn is_input_dir(&self) -> bool {
        self.matches.is_present("dir")
    }
}
