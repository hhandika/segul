use std::path::PathBuf;

use ansi_term::Colour::Yellow;
use clap::ArgMatches;

use crate::cli::{InputCli, InputPrint, OutputCli};
use crate::handler::extract::{Extract, ExtractOpts};
use crate::parser::txt;

impl InputCli for ExtractParser<'_> {}
impl InputPrint for ExtractParser<'_> {}
impl OutputCli for ExtractParser<'_> {}

pub(in crate::cli) struct ExtractParser<'a> {
    matches: &'a ArgMatches,
    input_dir: Option<PathBuf>,
    params: ExtractOpts,
}

impl<'a> ExtractParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches) -> Self {
        Self {
            matches,
            input_dir: None,
            params: ExtractOpts::None,
        }
    }

    pub(in crate::cli) fn extract(&mut self) {
        let input_fmt = self.parse_input_fmt(self.matches);
        let datatype = self.parse_datatype(self.matches);
        let output_fmt = self.parse_output_fmt(self.matches);
        let outdir = self.parse_output(self.matches);
        let task_desc = "Sequence extraction";
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

        let is_overwrite = self.parse_overwrite_opts(self.matches);
        self.check_output_dir_exist(&outdir, is_overwrite);
        log::info!("{}", Yellow.paint("ExtractOpts"));
        self.parse_params();
        let extract = Extract::new(&self.params, &input_fmt, &datatype);
        extract.extract_sequences(&files, &outdir, &output_fmt);
    }

    fn parse_params(&mut self) {
        match self.matches {
            m if m.is_present("regex") => {
                let re = self.parse_regex();
                log::info!("{:18}: {}\n", "Regex", re);
                self.params = ExtractOpts::Regex(re);
            }
            m if m.is_present("id") => {
                let ids = self.parse_id();
                log::info!("{:18}: {:?}\n", "IDs", ids);
                self.params = ExtractOpts::Id(ids);
            }
            m if m.is_present("file") => {
                let ids = self.parse_file();
                log::info!(
                    "{:18}: {}\n",
                    "File",
                    self.matches
                        .value_of("file")
                        .expect("Failed parsing file path")
                );
                self.params = ExtractOpts::Id(ids);
            }
            _ => unreachable!("Unknown parameters!"),
        }
    }

    fn parse_regex(&self) -> String {
        let re = self
            .matches
            .value_of("regex")
            .expect("Failed parsing regex string");
        String::from(re)
    }

    fn parse_file(&self) -> Vec<String> {
        let file = PathBuf::from(
            self.matches
                .value_of("file")
                .expect("Failed parsing file path"),
        );
        assert!(file.is_file(), "File does not exist: {}", file.display());
        txt::parse_text_file(&file)
    }

    fn parse_id(&self) -> Vec<String> {
        self.matches
            .values_of("id")
            .expect("Failed parsing IDs input")
            .map(String::from)
            .collect()
    }
}
