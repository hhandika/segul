use std::path::{Path, PathBuf};

use ansi_term::Colour::Yellow;
use clap::ArgMatches;

use crate::cli::{InputCli, InputPrint, OutputCli};
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
        let files = if self.matches.is_present("dir") {
            let dir = self.parse_dir_input(self.matches);
            self.input_dir = Some(PathBuf::from(dir));
            self.get_files(dir, &input_fmt)
        } else {
            self.parse_input(self.matches)
        };

        let ids = Path::new(
            self.matches
                .value_of("names")
                .expect("Failed parsing path to id names"),
        );

        self.print_input(
            &self.input_dir,
            task_desc,
            files.len(),
            &input_fmt,
            &datatype,
        );

        self.check_output_dir_exist(&outdir);
        log::info!("{}", Yellow.paint("Names"));
        log::info!(
            "{:18}: {}",
            "File",
            ids.file_name()
                .expect("Failed parsing name path")
                .to_string_lossy()
        );
        if self.matches.is_present("dry-run") {
            Rename::new(&input_fmt, &datatype, ids).dry_run();
        } else {
            Rename::new(&input_fmt, &datatype, ids).rename(&files, &outdir, &output_fmt);
        }
    }
}
