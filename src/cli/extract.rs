use ansi_term::Colour::Yellow;

use crate::cli::*;
use crate::core::extract::{Extract, Params};

impl InputCli for ExtractParser<'_> {}
impl InputPrint for ExtractParser<'_> {}
impl OutputCli for ExtractParser<'_> {}

pub(in crate::cli) struct ExtractParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_dir: Option<PathBuf>,
    params: Params,
}

impl<'a> ExtractParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            input_dir: None,
            params: Params::None,
        }
    }

    pub(in crate::cli) fn extract(&mut self) {
        let input_fmt = self.parse_input_fmt(self.matches);
        let datatype = self.parse_datatype(self.matches);
        let output_fmt = self.parse_output_fmt(self.matches);
        let outdir = self.parse_output(self.matches);
        self.parse_params();
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
        self.check_output_dir_exist(&outdir);
        self.print_input_info();
        let extract = Extract::new(&self.params, &input_fmt, &datatype);
        extract.extract_sequences(&files, &outdir, &output_fmt);
    }

    fn parse_params(&mut self) {
        match self.matches {
            m if m.is_present("regex") => self.params = Params::Regex(self.parse_regex()),
            m if m.is_present("id") => self.params = Params::Id(self.parse_id()),
            m if m.is_present("file") => self.params = Params::File(self.parse_file()),
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

    fn parse_file(&self) -> PathBuf {
        let file = PathBuf::from(
            self.matches
                .value_of("file")
                .expect("Failed parsing file path"),
        );
        assert!(file.is_file(), "File does not exist: {}", file.display());
        file
    }

    fn parse_id(&self) -> Vec<String> {
        self.matches
            .values_of("id")
            .expect("Failed parsing IDs input")
            .map(|s| String::from(s))
            .collect()
    }

    fn is_input_wcard(&self) -> bool {
        self.matches.is_present("wildcard")
    }

    fn print_input_info(&self) {
        log::info!("{}", Yellow.paint("Params"));
        match &self.params {
            Params::Regex(re) => log::info!("{:18}: {}\n", "Regex", re),
            Params::File(path) => log::info!("{:18}: {}\n", "File", path.display()),
            Params::Id(ids) => log::info!("{:18}: {:?}\n", "IDs", ids),
            Params::None => panic!("Please, specify a matching parameter!"),
        };
    }
}
