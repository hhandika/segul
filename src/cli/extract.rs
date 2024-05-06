use std::path::{Path, PathBuf};

use colored::Colorize;

use crate::handler::sequence::extract::{Extract, ExtractOpts};
use crate::helper::logger::AlignSeqLogger;
use crate::parser::txt;

use super::args::SequenceExtractArgs;
use super::{collect_paths, AlignSeqInput, InputCli, OutputCli};

impl InputCli for ExtractParser<'_> {}
impl OutputCli for ExtractParser<'_> {}
impl AlignSeqInput for ExtractParser<'_> {}

pub(in crate::cli) struct ExtractParser<'a> {
    args: &'a SequenceExtractArgs,
    input_dir: Option<PathBuf>,
    params: ExtractOpts,
}

impl<'a> ExtractParser<'a> {
    pub(in crate::cli) fn new(args: &'a SequenceExtractArgs) -> Self {
        Self {
            args,
            input_dir: None,
            params: ExtractOpts::None,
        }
    }

    pub(in crate::cli) fn extract(&mut self) {
        let input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let output_fmt = self.parse_output_fmt(&self.args.out_fmt.output_fmt);
        let task = "Sequence extraction";
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);
        AlignSeqLogger::new(
            self.input_dir.as_deref(),
            &input_fmt,
            &datatype,
            files.len(),
        )
        .log(task);
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        log::info!("{}", "ExtractOpts".yellow());
        self.parse_params();
        let extract = Extract::new(&self.params, &input_fmt, &datatype);
        extract.extract_sequences(&files, &self.args.output, &output_fmt);
    }

    fn parse_params(&mut self) {
        if let Some(re) = &self.args.re {
            log::info!("{:18}: {}\n", "Regex", re);
            self.params = ExtractOpts::Regex(re.clone());
        }

        if let Some(id) = &self.args.id {
            let id_list = self.parse_id_opts(id);
            log::info!("{:18}: {:?}\n", "IDs", &id_list);
            self.params = ExtractOpts::Id(id_list);
        }

        if let Some(file) = &self.args.file {
            let ids = self.parse_file(file);
            log::info!(
                "{:18}: {}\n",
                "File",
                self.args
                    .file
                    .as_ref()
                    .expect("Failed parsing file path")
                    .display()
            );
            self.params = ExtractOpts::Id(ids);
        }
    }

    fn parse_file(&self, file: &Path) -> Vec<String> {
        assert!(file.is_file(), "File does not exist: {}", file.display());
        txt::parse_text_file(file)
    }

    fn parse_id_opts(&self, id_input: &str) -> Vec<String> {
        let id_list: Vec<String> = id_input.split(';').map(|s| s.trim().to_string()).collect();
        if id_list.is_empty() {
            panic!("Failed parsing the ID input. Make sure you use semicolon to separate the IDs");
        }

        id_list
    }
}
