use std::path::PathBuf;

use colored::Colorize;

use crate::{
    core::sequence::remove::{SeqRemovalParameters, SequenceRemoval},
    helper::logger::AlignSeqLogger,
};

use crate::cli::{
    args::sequence::SequenceRemoveArgs, collect_paths, AlignSeqInput, InputCli, OutputCli,
};

impl InputCli for RemoveParser<'_> {}
impl OutputCli for RemoveParser<'_> {}
impl AlignSeqInput for RemoveParser<'_> {}

pub(in crate::cli) struct RemoveParser<'a> {
    args: &'a SequenceRemoveArgs,
    input_dir: Option<PathBuf>,
}

impl<'a> RemoveParser<'a> {
    pub(in crate::cli) fn new(args: &'a SequenceRemoveArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn remove(&mut self) {
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
        let opts = self.parse_remove_opts();
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        SequenceRemoval::new(&input_fmt, &datatype, &self.args.output, &output_fmt, &opts)
            .remove(&files);
    }

    fn parse_remove_opts(&self) -> SeqRemovalParameters {
        log::info!("{}", "Params".yellow());
        if let Some(re) = &self.args.re {
            log::info!("{:18}: {}\n", "Regex", "Options");
            log::info!("{:18}, {}\n", "Values", re);
            SeqRemovalParameters::Regex(re.clone())
        } else if let Some(ids) = &self.args.id {
            let id_list = self.parse_id_opts(ids);
            log::info!("{:18}: id", "Options");
            log::info!("{:18}: {:?}", "Values", &id_list);
            SeqRemovalParameters::Id(id_list)
        } else {
            unimplemented!("SeqRemovalParameters::None is not implemented yet")
        }
    }

    fn parse_id_opts(&self, id_input: &str) -> Vec<String> {
        let id_list: Vec<String> = id_input.split(';').map(|s| s.trim().to_string()).collect();
        if id_list.is_empty() {
            panic!("Failed parsing the ID input. Make sure you use semicolon to separate the IDs");
        }

        id_list
    }
}
