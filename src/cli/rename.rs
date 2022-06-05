use std::path::{Path, PathBuf};

use ansi_term::Colour::Yellow;
use clap::ArgMatches;

use crate::cli::{InputCli, InputPrint, OutputCli};
use crate::handler::rename::{Rename, RenameOpts};
use crate::helper::utils;
use crate::parser::delimited;

impl InputCli for RenameParser<'_> {}
impl InputPrint for RenameParser<'_> {}
impl OutputCli for RenameParser<'_> {}

pub(in crate::cli) struct RenameParser<'a> {
    matches: &'a ArgMatches,
    input_dir: Option<PathBuf>,
}

impl<'a> RenameParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches) -> Self {
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

        self.print_input(
            &self.input_dir,
            task_desc,
            files.len(),
            &input_fmt,
            &datatype,
        );
        let params = self.parse_rename_opts();
        let is_overwrite = self.parse_overwrite_opts(self.matches);
        self.check_output_dir_exist(&outdir, is_overwrite);
        if self.matches.is_present("dry-run") {
            Rename::new(&input_fmt, &datatype, &outdir, &output_fmt).dry_run(&params);
        } else {
            Rename::new(&input_fmt, &datatype, &outdir, &output_fmt).rename(&files, &params);
        }
    }

    fn parse_rename_opts(&self) -> RenameOpts {
        log::info!("{}", Yellow.paint("Params"));
        match self.matches {
            m if m.is_present("names") => {
                let id_path = Path::new(
                    self.matches
                        .value_of("names")
                        .expect("Failed parsing path to id names"),
                );
                let names = self.parse_names(id_path);
                self.print_rename_id_info(id_path, &names.len());
                RenameOpts::RnId(names)
            }
            m if m.is_present("rm-string") => {
                let input_str = self
                    .matches
                    .value_of("rm-string")
                    .expect("Failed parsing input string");
                self.print_remove_str_info(input_str);
                RenameOpts::RmStr(input_str.to_string())
            }
            _ => unreachable!("Unknown errors in parsing command line input!"),
        }
    }

    fn parse_names(&self, id_path: &Path) -> Vec<(String, String)> {
        delimited::parse_delimited_text(id_path)
    }

    fn print_rename_id_info(&self, id_path: &Path, id_count: &usize) {
        log::info!("{:18}: --names", "Options");
        log::info!(
            "{:18}: {}",
            "File",
            id_path
                .file_name()
                .expect("Failed parsing name path")
                .to_string_lossy()
        );
        log::info!("{:18}: {}\n", "New ID count", utils::fmt_num(id_count));
    }

    fn print_remove_str_info(&self, input_str: &str) {
        log::info!("{:18}: --rm-string", "Options");
        log::info!("{:18}: {}\n", "Input string", input_str);
    }
}
