use crate::core::converter::Converter;
use clap::ArgMatches;

use crate::cli::*;
use crate::helper::types::{DataType, InputFmt, OutputFmt};

impl InputCli for ConvertParser<'_> {}
impl OutputCli for ConvertParser<'_> {}

pub(in crate::cli) struct ConvertParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_fmt: InputFmt,
    output: PathBuf,
    output_fmt: OutputFmt,
    datatype: DataType,
    is_dir: bool,
}

impl<'a> ConvertParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            input_fmt: InputFmt::Auto,
            datatype: DataType::Dna,
            output: PathBuf::new(),
            output_fmt: OutputFmt::Nexus,
            is_dir: false,
        }
    }

    pub(in crate::cli) fn convert(&mut self) {
        self.input_fmt = self.parse_input_fmt(self.matches);
        self.output = self.parse_output_path(self.matches);
        self.output_fmt = self.parse_output_fmt(self.matches);
        self.datatype = self.parse_datatype(self.matches);
        let input_type = self.parse_input_type(self.matches);
        match input_type {
            InputType::File => self.convert_file(),
            InputType::Dir => {
                let dir = self.parse_dir_input(self.matches);
                let files = self.get_files(dir, &self.input_fmt);
                self.print_input_dir(Some(dir), files.len(), &self.output);
                self.convert_multiple_files(&files);
            }
            InputType::Wildcard => {
                let files = self.parse_input_wcard(&self.matches);
                self.print_input_dir::<PathBuf>(None, files.len(), &self.output);
                self.convert_multiple_files(&files)
            }
        }
    }

    fn convert_file(&mut self) {
        let input = Path::new(self.parse_file_input(self.matches));
        self.print_input_file(input);
        self.convert_any(input, &self.output, &self.output_fmt);
    }

    fn convert_multiple_files(&mut self, files: &[PathBuf]) {
        self.is_dir = true;
        let spin = utils::set_spinner();
        spin.set_message("Converting alignments...");
        files.par_iter().for_each(|file| {
            let output = self.output.join(file.file_stem().unwrap());
            self.convert_any(file, &output, &self.output_fmt);
        });
        spin.finish_with_message("DONE!");
    }

    fn convert_any(&self, input: &Path, output: &Path, output_fmt: &OutputFmt) {
        let mut convert = Converter::new(input, output, output_fmt, &self.datatype);
        convert.set_isdir(self.is_dir);
        if self.is_sort() {
            convert.convert_sorted(&self.input_fmt);
        } else {
            convert.convert_unsorted(&self.input_fmt);
        }
    }

    fn is_sort(&self) -> bool {
        self.matches.is_present("sort")
    }

    fn print_input_file(&self, input: &Path) {
        log::info!("{:18}: {}", "Input", &input.display());
        log::info!("{:18}: Sequence format conversion", "Task");
    }

    fn print_input_dir<P: AsRef<Path>>(&self, input: Option<P>, nfile: usize, output: &Path) {
        if let Some(input) = input {
            log::info!("{:18}: {}", "Input dir", &input.as_ref().display());
        } else {
            log::info!("{:18}: {}", "Input dir", "WILDCARD");
        }
        log::info!("{:18}: {}", "Total files", utils::fmt_num(&nfile));
        log::info!("{:18}: Sequence format conversion\n", "Task");
        log::info!("{:18}: {}", "Output dir", output.display());
    }
}
