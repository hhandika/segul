use std::ffi::OsStr;
use std::path::{Path, PathBuf};

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

        let output = self.parse_output(self.matches);
        let output = output.with_extension("txt");
        let overwrite = self.parse_overwrite_opts(self.matches);
        self.check_output_file_exist(&output, overwrite);
        let id = Id::new(&output, &input_fmt, &datatype);
        if self.matches.is_present("map") {
            let map_fname = self.create_map_fname(&output);
            self.check_output_file_exist(&map_fname, overwrite);
            id.map_id(&files, &map_fname);
        } else {
            id.generate_id(&files);
        }
    }

    fn create_map_fname(&self, output: &Path) -> PathBuf {
        let parent = output.parent().expect("Failed getting parent dir");
        let fstem = output
            .file_stem()
            .and_then(OsStr::to_str)
            .expect("Failed getting file stem for mapping IDs");
        parent.join(format!("{}_map", fstem)).with_extension("csv")
    }
}
