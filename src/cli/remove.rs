use std::path::PathBuf;

use ansi_term::Colour::Yellow;
use clap::ArgMatches;

// use crate::helper::utils;
use crate::{
    cli::{InputCli, InputPrint, OutputCli},
    handler::remove::{Remove, RemoveOpts},
};

impl InputCli for RemoveParser<'_> {}
impl InputPrint for RemoveParser<'_> {}
impl OutputCli for RemoveParser<'_> {}

pub(in crate::cli) struct RemoveParser<'a> {
    matches: &'a ArgMatches,
    input_dir: Option<PathBuf>,
}

impl<'a> RemoveParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches) -> Self {
        Self {
            matches,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn remove(&mut self) {
        let input_fmt = self.parse_input_fmt(self.matches);
        let datatype = self.parse_datatype(self.matches);
        let output_fmt = self.parse_output_fmt(self.matches);
        let outdir = self.parse_output(self.matches);
        let task_desc = "Sequence Renaming";
        let files = if self.matches.is_present("dir") {
            let dir = self.parse_dir_input(self.matches);
            self.input_dir = Some(PathBuf::from(dir));
            self.get_files(dir, &input_fmt)
        } else {
            self.parse_input(self.matches)
        };

        self.print_input(
            &self.input_dir,
            task_desc,
            files.len(),
            &input_fmt,
            &datatype,
        );
        let opts = self.parse_remove_opts();

        let is_overwrite = self.parse_overwrite_opts(self.matches);
        self.check_output_dir_exist(&outdir, is_overwrite);
        Remove::new(&input_fmt, &datatype, &outdir, &output_fmt, &opts).remove(&files);
    }

    fn parse_remove_opts(&self) -> RemoveOpts {
        log::info!("{}", Yellow.paint("Params"));
        match self.matches {
            m if m.is_present("id") => {
                let ids = self
                    .matches
                    .values_of("id")
                    .expect("Failed parsing ids")
                    .map(String::from)
                    .collect();
                log::info!("{:18}: id", "Options");
                log::info!("{:18}, {:?}", "Values", ids);
                RemoveOpts::Id(ids)
            }

            m if m.is_present("regex") => {
                let input_re = self
                    .matches
                    .value_of("re")
                    .expect("Failed parsing regex values")
                    .to_string();
                RemoveOpts::Regex(input_re)
            }
            _ => unimplemented!("Unknown errors in parsing command line input!"),
        }
    }
}
