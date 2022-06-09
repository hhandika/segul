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
        let opts = self.parse_rename_opts();
        if self.matches.is_present("dry-run") {
            Rename::new(&input_fmt, &datatype, &opts).dry_run(&files);
        } else {
            let is_overwrite = self.parse_overwrite_opts(self.matches);
            self.check_output_dir_exist(&outdir, is_overwrite);
            Rename::new(&input_fmt, &datatype, &opts).rename(&files, &outdir, &output_fmt);
        }
    }

    fn parse_rename_opts(&self) -> RenameOpts {
        log::info!("{}", Yellow.paint("Params"));
        match self.matches {
            m if m.is_present("replace-id") => {
                let id_path = Path::new(
                    self.matches
                        .value_of("replace-id")
                        .expect("Failed parsing path to id names"),
                );
                let names = self.parse_names(id_path);
                self.print_rename_id_info(id_path, &names.len());
                RenameOpts::RnId(names)
            }
            m if m.is_present("remove") => {
                let input_str = self
                    .matches
                    .value_of("rm-string")
                    .expect("Failed parsing input string");
                self.print_remove_str_info(input_str);
                RenameOpts::RmStr(input_str.to_string())
            }
            m if m.is_present("remove-re") => {
                let input_re = self
                    .matches
                    .value_of("remove-re")
                    .expect("Failed parsing input regex");
                let is_all = false;
                self.print_remove_re_info(input_re, "--remove-re");
                RenameOpts::RmRegex(input_re.to_string(), is_all)
            }
            m if m.is_present("remove-re-all") => {
                let input_re = self
                    .matches
                    .value_of("remove-re-all")
                    .expect("Failed parsing input regex");
                let is_all = true;
                self.print_remove_re_info(input_re, "--remove-re-all");
                RenameOpts::RmRegex(input_re.to_string(), is_all)
            }
            m if m.is_present("replace-from") => {
                let input_str = self
                    .matches
                    .value_of("replace-from")
                    .expect("Failed parsing input string");
                let output_str = self
                    .matches
                    .value_of("replace-to")
                    .expect("Failed parsing output string");
                self.print_replace_str_info(input_str, output_str);
                RenameOpts::RpStr(input_str.to_string(), output_str.to_string())
            }
            m if m.is_present("replace-from-re") => {
                let input_re = self
                    .matches
                    .value_of("replace-from-re")
                    .expect("Failed parsing input regex");
                let output_str = self
                    .matches
                    .value_of("replace-to")
                    .expect("Failed parsing output string");
                let is_all = false;
                self.print_replace_re_info(input_re, output_str, "--replace-from-re");
                RenameOpts::RpRegex(input_re.to_string(), output_str.to_string(), is_all)
            }
            _ => unreachable!("Unknown errors in parsing command line input!"),
        }
    }

    fn parse_names(&self, id_path: &Path) -> Vec<(String, String)> {
        delimited::parse_delimited_text(id_path)
    }

    fn print_rename_id_info(&self, id_path: &Path, id_count: &usize) {
        log::info!("{:18}: --replace", "Options");
        log::info!(
            "{:18}: {}",
            "File",
            id_path
                .file_name()
                .expect("Failed parsing name path")
                .to_string_lossy()
        );
        log::info!("{:18}: {}\n", "New ID counts", utils::fmt_num(id_count));
    }

    fn print_remove_str_info(&self, input_str: &str) {
        log::info!("{:18}: --remove", "Options");
        log::info!("{:18}: {}\n", "Input string", input_str);
    }

    fn print_replace_str_info(&self, input_str: &str, output_str: &str) {
        log::info!("{:18}: --replace", "Options");
        log::info!("{:18}: {}", "Replace from", input_str);
        log::info!("{:18}: {}\n", "Replace to", output_str);
    }

    fn print_remove_re_info(&self, input_re: &str, options: &str) {
        log::info!("{:18}: {}", "Options", options);
        log::info!("{:18}: {}\n", "Input regex", input_re);
    }

    fn print_replace_re_info(&self, input_re: &str, output_str: &str, options: &str) {
        log::info!("{:18}: {}", "Options", options);
        log::info!("{:18}: {}", "Replace from", input_re);
        log::info!("{:18}: {}\n", "Replace to", output_str);
    }
}
