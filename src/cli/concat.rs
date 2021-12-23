use crate::cli::*;
use crate::helper::filenames;

impl ConcatCLi for ConcatParser<'_> {}
impl InputPrint for ConcatParser<'_> {}

impl InputCli for ConcatParser<'_> {
    fn parse_input_type(&self, matches: &ArgMatches) -> InputType {
        if matches.is_present("dir") {
            InputType::Dir
        } else {
            InputType::Wildcard
        }
    }
}

impl OutputCli for ConcatParser<'_> {}

pub(in crate::cli) struct ConcatParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_dir: Option<PathBuf>,
}

impl<'a> ConcatParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn concat(&mut self) {
        let input_fmt = self.parse_input_fmt(self.matches);
        let datatype = self.parse_datatype(self.matches);
        let output_fmt = self.parse_output_fmt(self.matches);
        let dir = self.parse_output(self.matches);
        let prefix = self.parse_prefix(self.matches, &dir);
        // let final_path = dir.join(prefix);
        let output = filenames::create_output_fname(&dir, &prefix, &output_fmt);
        let part_fmt = self.parse_partition_fmt(self.matches);
        self.check_partition_format(&output_fmt, &part_fmt);
        let task_desc = "Alignment concatenation";
        let mut files = if self.is_input_wcard() {
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

        self.check_output_dir_exist(&dir);
        let mut concat = ConcatHandler::new(&input_fmt, &output, &output_fmt, &part_fmt);

        concat.concat_alignment(&mut files, &datatype);
    }

    fn is_input_wcard(&self) -> bool {
        self.matches.is_present("wildcard")
    }
}
