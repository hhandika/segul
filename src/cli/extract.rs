use crate::cli::*;
use crate::core::extract::Extract;

impl InputCli for ExtractParser<'_> {}

pub(in crate::cli) struct ExtractParser<'a> {
    matches: &'a ArgMatches<'a>,
}

impl<'a> ExtractParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self { matches }
    }

    pub(in crate::cli) fn extract(&self) {
        println!("Extract {}", self.parse_dir_input(self.matches).display());
        let re = self.parse_regex();
        Extract::new(re).print_input_info();
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
}
