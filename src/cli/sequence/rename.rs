use std::path::{Path, PathBuf};

use colored::Colorize;

use crate::core::sequence::rename::{SeqRenamingParameters, SequenceRenaming, SequenceRenamingDry};
use crate::helper::logger::AlignSeqLogger;
use crate::helper::utils;
use crate::parser::delimited;

use crate::cli::args::sequence::SequenceRenameArgs;
use crate::cli::{collect_paths, AlignSeqInput, InputCli, OutputCli};

impl InputCli for RenameParser<'_> {}
impl OutputCli for RenameParser<'_> {}
impl AlignSeqInput for RenameParser<'_> {}

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
        let task = "Sequence Renaming";
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);
        AlignSeqLogger::new(
            self.input_dir.as_deref(),
            &input_fmt,
            &datatype,
            files.len(),
        )
        .log(task);
        let opts = self.parse_rename_opts();
        if self.args.dry_run {
            SequenceRenamingDry::new(&input_fmt, &datatype, &opts).dry_run(&files);
        } else {
            self.check_output_dir_exist(&self.args.output, self.args.io.force);
            SequenceRenaming::new(&input_fmt, &datatype, &self.args.output, &output_fmt, &opts)
                .rename(&files);
        }
    }

    fn parse_rename_opts(&self) -> SeqRenamingParameters {
        log::info!("{}", "Renaming Parameters".yellow());

        if let Some(path) = &self.args.replace_id {
            let id_path = Path::new(&path);
            let names = self.parse_names(id_path);
            self.print_rename_id_info(id_path, &names.len());
            SeqRenamingParameters::RnId(names)
        } else if let Some(input_str) = &self.args.remove {
            self.print_remove_str_info(input_str);
            SeqRenamingParameters::RmStr(input_str.to_string())
        } else if let Some(re) = &self.args.remove_re {
            self.print_remove_re_info(re, "--remove-re");
            SeqRenamingParameters::RmRegex(re.to_string(), false)
        } else if let Some(re) = &self.args.remove_re_all {
            let is_all = true;
            self.print_remove_re_info(re, "--remove-re-all");
            SeqRenamingParameters::RmRegex(re.to_string(), is_all)
        } else if let Some(id) = &self.args.replace_from {
            match &self.args.replace_to {
                Some(to) => {
                    self.print_replace_str_info(id, to);
                    SeqRenamingParameters::RpStr(id.to_string(), to.to_string())
                }
                None => unreachable!("Missing replace-to"),
            }
        } else if let Some(re) = &self.args.replace_from_re {
            match &self.args.replace_to {
                Some(to) => {
                    self.print_replace_re_info(re, to, "--replace-from-re");
                    SeqRenamingParameters::RpRegex(re.to_string(), to.to_string(), false)
                }
                None => unreachable!("Missing replace-to"),
            }
        } else if let Some(re) = &self.args.replace_from_re_all {
            match &self.args.replace_to {
                Some(to) => {
                    self.print_replace_re_info(re, to, "--replace-from-re-all");
                    SeqRenamingParameters::RpRegex(re.to_string(), to.to_string(), true)
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
