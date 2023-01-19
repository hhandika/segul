use std::path::{Path, PathBuf};

use colored::Colorize;

use crate::cli::{InputCli, InputPrint, OutputCli};
use crate::handler::rename::{Rename, RenameDry, RenameOpts};
use crate::helper::utils;
use crate::parser::delimited;

use super::args::SequenceRenameArgs;
use super::collect_paths;

impl InputCli for RenameParser<'_> {}
impl InputPrint for RenameParser<'_> {}
impl OutputCli for RenameParser<'_> {}

pub(in crate::cli) struct RenameParser<'a> {
    args: &'a SequenceRenameArgs,
    input_dir: Option<PathBuf>,
}

impl<'a> RenameParser<'a> {
    pub(in crate::cli) fn new(args: &'a SequenceRenameArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn rename(&mut self) {
        let input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let output_fmt = self.parse_output_fmt(&self.args.out_fmt.output_fmt);
        let task_desc = "Sequence Renaming";
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);
        self.print_input(
            &self.input_dir,
            task_desc,
            files.len(),
            &input_fmt,
            &datatype,
        );
        let opts = self.parse_rename_opts();
        if self.args.dry_run {
            RenameDry::new(&input_fmt, &datatype, &opts).dry_run(&files);
        } else {
            self.check_output_dir_exist(&self.args.output, self.args.io.force);
            Rename::new(&input_fmt, &datatype, &self.args.output, &output_fmt, &opts)
                .rename(&files);
        }
    }

    fn parse_rename_opts(&self) -> RenameOpts {
        log::info!("{}", "Params".yellow());

        if let Some(path) = &self.args.replace_id {
            let id_path = Path::new(&path);
            let names = self.parse_names(id_path);
            self.print_rename_id_info(id_path, &names.len());
            RenameOpts::RnId(names)
        } else if let Some(input_str) = &self.args.remove {
            self.print_remove_str_info(&input_str);
            RenameOpts::RmStr(input_str.to_string())
        } else if let Some(input_re) = &self.args.remove_re {
            self.print_remove_re_info(&input_re, "--remove-re");
            RenameOpts::RmRegex(input_re.to_string(), false)
        } else if let Some(re) = &self.args.remove_re_all {
            let is_all = true;
            self.print_remove_re_info(&re, "--remove-re-all");
            RenameOpts::RmRegex(re.to_string(), is_all)
        } else if let Some(id) = &self.args.replace_from {
            match &self.args.replace_to {
                Some(to) => {
                    self.print_replace_str_info(&id, &to);
                    RenameOpts::RpStr(id.to_string(), to.to_string())
                }
                None => unreachable!("Missing replace-to"),
            }
        } else if let Some(re) = &self.args.replace_from_re {
            match &self.args.replace_to {
                Some(to) => {
                    self.print_replace_re_info(&re, &to, "--replace-from-re");
                    RenameOpts::RpRegex(re.to_string(), to.to_string(), false)
                }
                None => unreachable!("Missing replace-to"),
            }
        } else if let Some(re) = &self.args.replace_from_re_all {
            match &self.args.replace_to {
                Some(to) => {
                    self.print_replace_re_info(&re, &to, "--replace-from-re-all");
                    RenameOpts::RpRegex(re.to_string(), to.to_string(), true)
                }
                None => unreachable!("Missing replace-to"),
            }
        } else {
            unreachable!("Missing rename options")
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
