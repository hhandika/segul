use std::path::PathBuf;

use colored::Colorize;

use crate::{
    cli::{InputCli, InputPrint, OutputCli},
    handler::remove::{Remove, RemoveOpts},
};

use super::{args::SequenceRemoveArgs, collect_paths};

impl InputCli for RemoveParser<'_> {}
impl InputPrint for RemoveParser<'_> {}
impl OutputCli for RemoveParser<'_> {}

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
        let opts = self.parse_remove_opts();

        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        Remove::new(&input_fmt, &datatype, &self.args.output, &output_fmt, &opts).remove(&files);
    }

    fn parse_remove_opts(&self) -> RemoveOpts {
        log::info!("{}", "Params".yellow());
        if let Some(re) = &self.args.re {
            log::info!("{:18}: {}\n", "Regex", "Options");
            log::info!("{:18}, {}\n", "Values", re);
            return RemoveOpts::Regex(re.clone());
        }
        if let Some(ids) = &self.args.id {
            log::info!("{:18}: id", "Options");
            log::info!("{:18}, {:?}", "Values", ids);
            return RemoveOpts::Id(ids.clone());
        } else {
            unimplemented!("RemoveOpts::None is not implemented yet")
        }
    }
}
