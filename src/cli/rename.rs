use crate::cli::*;
use crate::core::rename::Rename;

impl InputCli for RenameParser<'_> {}
impl InputPrint for RenameParser<'_> {}
impl OutputCli for RenameParser<'_> {}

pub(in crate::cli) struct RenameParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_dir: Option<PathBuf>,
}

impl<'a> RenameParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn rename(&mut self) {
        let input_fmt = self.parse_input_fmt(self.matches);
        let datatype = self.parse_datatype(self.matches);
        let output_fmt = self.parse_output_fmt(self.matches);
        let outdir = self.parse_output(self.matches);
        let task_desc = "Sequence Renaming";
        let files = if self.matches.is_present("wildcard") {
            self.parse_input_wcard(self.matches)
        } else {
            let dir = self.parse_dir_input(self.matches);
            self.input_dir = Some(PathBuf::from(dir));
            self.get_files(dir, &input_fmt)
        };

        let ids = self
            .matches
            .value_of("names")
            .expect("Failed parsing path to id names");

        self.print_input_multi(
            &self.input_dir,
            task_desc,
            files.len(),
            &input_fmt,
            &datatype,
        );

        self.check_output_dir_exist(&outdir);
        Rename::new(&input_fmt, &datatype).rename(&ids, &outdir, &output_fmt);
    }
}
