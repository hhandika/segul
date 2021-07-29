use clap::ArgMatches;

use crate::cli::*;
use crate::helper::types::{DataType, InputFmt, OutputFmt};

impl PartCLi for ConcatParser<'_> {}
impl InputPrint for ConcatParser<'_> {}

impl InputCli for ConcatParser<'_> {
    fn parse_input_type(&self, matches: &ArgMatches) -> InputType {
        if matches.is_present("dir") {
            InputType::Dir
        } else {
            InputType::Wildcard
        }
    }
}

impl OutputCli for ConcatParser<'_> {
    fn parse_output_path(&self, matches: &ArgMatches) -> PathBuf {
        PathBuf::from(self.parse_output(matches))
    }
}

pub(in crate::cli) struct ConcatParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_fmt: InputFmt,
    input_dir: Option<PathBuf>,
    input_type: InputType,
    output_fmt: OutputFmt,
    part_fmt: PartitionFmt,
    datatype: DataType,
}

impl<'a> ConcatParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            input_fmt: InputFmt::Fasta,
            input_dir: None,
            input_type: InputType::Dir,
            output_fmt: OutputFmt::Nexus,
            part_fmt: PartitionFmt::Charset,
            datatype: DataType::Dna,
        }
    }

    pub(in crate::cli) fn concat(&mut self) {
        self.input_fmt = self.parse_input_fmt(self.matches);
        self.input_type = self.parse_input_type(self.matches);
        self.datatype = self.parse_datatype(self.matches);
        let output = self.parse_output(self.matches);
        self.output_fmt = self.parse_output_fmt(self.matches);
        self.part_fmt = self.parse_partition_fmt(self.matches);
        self.check_partition_format(&self.output_fmt, &self.part_fmt);
        let task_desc = "Alignment concatenation";
        let mut files = if self.is_input_wcard() {
            self.parse_input_wcard(&self.matches)
        } else {
            let dir = self.parse_dir_input(self.matches);
            self.input_dir = Some(PathBuf::from(dir));
            self.get_files(dir, &self.input_fmt)
        };
        self.print_input_multi(&self.input_dir, task_desc, files.len(), &self.input_fmt);
        let concat =
            msa::MSAlignment::new(&self.input_fmt, output, &self.output_fmt, &self.part_fmt);

        concat.concat_alignment(&mut files, &self.datatype);
    }

    fn is_input_wcard(&self) -> bool {
        self.matches.is_present("wildcard")
    }
}
