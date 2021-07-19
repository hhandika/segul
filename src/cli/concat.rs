use clap::ArgMatches;

use crate::cli::*;
use crate::helper::common::{DataType, InputFmt, OutputFmt};

impl PartCLi for ConcatParser<'_> {}

impl Cli for ConcatParser<'_> {
    fn get_input_type(&self, matches: &ArgMatches) -> InputType {
        if matches.is_present("dir") {
            InputType::Dir
        } else {
            InputType::Wildcard
        }
    }

    fn get_output_path(&self, matches: &ArgMatches) -> PathBuf {
        PathBuf::from(self.get_output(matches))
    }
}

pub(in crate::cli) struct ConcatParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_fmt: InputFmt,
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
            input_type: InputType::Dir,
            output_fmt: OutputFmt::Nexus,
            part_fmt: PartitionFmt::Charset,
            datatype: DataType::Dna,
        }
    }

    pub(in crate::cli) fn concat(&mut self) {
        self.input_fmt = self.get_input_fmt(self.matches);
        self.input_type = self.get_input_type(self.matches);
        self.datatype = self.get_datatype(self.matches);
        let output = self.get_output(self.matches);
        self.output_fmt = self.get_output_fmt(self.matches);
        self.part_fmt = self.parse_partition_fmt(self.matches);
        self.check_partition_format(&self.output_fmt, &self.part_fmt);
        let mut files = if self.is_input_wcard() {
            self.parse_input_wcard(&self.matches)
        } else {
            let dir = self.get_dir_input(self.matches);
            self.get_files(dir, &self.input_fmt)
        };
        self.print_user_input().unwrap();
        let concat =
            msa::MSAlignment::new(&self.input_fmt, output, &self.output_fmt, &self.part_fmt);

        concat.concat_alignment(&mut files, &self.datatype);
    }

    fn is_input_wcard(&self) -> bool {
        self.matches.is_present("wildcard")
    }

    fn print_user_input(&self) -> Result<()> {
        let io = io::stdout();
        let mut writer = io::BufWriter::new(io);
        writeln!(writer, "Command\t\t: segul concat")?;
        if !self.is_input_wcard() {
            writeln!(
                writer,
                "Input dir\t: {}\n",
                self.get_dir_input(self.matches)
            )?;
        } else {
            writeln!(writer, "Input\t\t: WILDCARD\n",)?;
        }
        Ok(())
    }
}
