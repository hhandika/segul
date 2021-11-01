use std::path::PathBuf;

use clap::ArgMatches;

use crate::cli::*;
use crate::core::id::Id;

impl InputCli for IdParser<'_> {}

impl OutputCli for IdParser<'_> {}

impl InputPrint for IdParser<'_> {}

pub(in crate::cli) struct IdParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_dir: Option<PathBuf>,
    output: PathBuf,
    files: Vec<PathBuf>,
}

impl<'a> IdParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches) -> Self {
        Self {
            matches,
            input_dir: None,
            output: PathBuf::new(),
            files: Vec::new(),
        }
    }

    pub(in crate::cli) fn get_id(&mut self) {
        let input_fmt = self.parse_input_fmt(self.matches);
        let datatype = self.parse_datatype(self.matches);
        let task_desc = "IDs finding";
        self.files = if self.is_input_dir() {
            let dir = self.parse_dir_input(self.matches);
            self.input_dir = Some(PathBuf::from(dir));
            self.get_files(dir, &input_fmt)
        } else {
            self.parse_input_wcard(self.matches)
        };

        self.print_input_multi(
            &self.input_dir,
            task_desc,
            self.files.len(),
            &input_fmt,
            &datatype,
        );

        self.output = self.parse_output(self.matches).with_extension("txt");
        let id = Id::new(&self.output, &input_fmt, &datatype);
        self.check_output_file_exist(&self.output);
        // if self.matches.is_present("map") {
        //     // id.map_id(&self.files);
        // } else {
        id.generate_id(&self.files);
        // }
    }

    fn is_input_dir(&self) -> bool {
        self.matches.is_present("dir")
    }
}
