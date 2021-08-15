use ansi_term::Colour::Yellow;

use crate::cli::*;
use crate::core::extract::Extract;

impl InputCli for ExtractParser<'_> {}
impl InputPrint for ExtractParser<'_> {}
impl OutputCli for ExtractParser<'_> {}

pub(in crate::cli) struct ExtractParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_dir: Option<PathBuf>,
    regex: Option<String>,
}

impl<'a> ExtractParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            input_dir: None,
            regex: None,
        }
    }

    pub(in crate::cli) fn extract(&mut self) {
        let input_fmt = self.parse_input_fmt(self.matches);
        let datatype = self.parse_datatype(self.matches);
        // let output_fmt = self.parse_output_fmt(self.matches);
        // let dir = self.parse_output(self.matches);
        self.regex = self.parse_regex();
        let task_desc = "Sequence extraction";
        let files = if self.is_input_wcard() {
            self.parse_input_wcard(self.matches)
        } else {
            let dir = self.parse_dir_input(self.matches);
            self.input_dir = Some(PathBuf::from(dir));
            self.get_files(dir, &input_fmt)
        };

        self.print_input_multi(
            &self.input_dir,
            task_desc,
            files.len(),
            &input_fmt,
            &datatype,
        );
        self.print_input_info();
        Extract::new(&self.regex).extract_sequences(&files);
    }

    fn parse_regex(&self) -> Option<String> {
        if self.matches.is_present("regex") {
            let re = self
                .matches
                .value_of("regex")
                .expect("Failed parsing regex string");
            Some(String::from(re))
        } else {
            None
        }
    }

    fn is_input_wcard(&self) -> bool {
        self.matches.is_present("wildcard")
    }

    fn print_input_info(&self) {
        log::info!("{}", Yellow.paint("Params"));
        match &self.regex {
            Some(re) => log::info!("{:18}: {}\n", "Regex", re),
            None => log::info!("{:18}: {}\n", "Regex", "None"),
        };
    }
}
