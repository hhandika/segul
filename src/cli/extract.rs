use std::path::{Path, PathBuf};

use colored::Colorize;

use crate::cli::AlignSeqPrint;
use crate::handler::sequence::extract::{Extract, ExtractOpts};
use crate::parser::txt;

use super::args::SequenceExtractArgs;
use super::{collect_paths, AlignSeqInput, InputCli, InputPrint, OutputCli};

impl InputCli for ExtractParser<'_> {}
impl InputPrint for ExtractParser<'_> {}
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
        let task_desc = "Sequence extraction";
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);
        AlignSeqPrint::new(
            &self.input_dir,
            &input_fmt,
            &datatype,
            task_desc,
            files.len(),
        )
        .print();
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
            log::info!("{:18}: {:?}\n", "IDs", id);
            self.params = ExtractOpts::Id(id.clone());
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
}
