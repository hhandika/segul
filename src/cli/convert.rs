use crate::core::converter::Converter;
use clap::ArgMatches;

use crate::cli::*;
use crate::helper::types::{DataType, InputFmt, OutputFmt};

impl InputCli for ConvertParser<'_> {}
impl InputPrint for ConvertParser<'_> {}
impl OutputCli for ConvertParser<'_> {
    fn parse_output<'a>(&self, matches: &'a ArgMatches) -> PathBuf {
        if matches.is_present("output") {
            PathBuf::from(
                matches
                    .value_of("output")
                    .expect("Failed parsing an output values"),
            )
        } else {
            PathBuf::from(self.parse_file_input(matches).file_name().expect(
                "Faile parsing input file to get the ouput name. Please specify output names!",
            ))
        }
    }
}

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
        self.output = self.parse_output(self.matches);
        self.output_fmt = self.parse_output_fmt(self.matches);
        self.datatype = self.parse_datatype(self.matches);
        let input_type = self.parse_input_type(self.matches);
        let task_desc = "Sequence format conversion";
        match input_type {
            InputType::File => self.convert_file(task_desc),
            InputType::Dir => {
                let dir = self.parse_dir_input(self.matches);
                let files = self.get_files(dir, &self.input_fmt);
                self.print_input_multi(
                    &Some(dir),
                    task_desc,
                    files.len(),
                    &self.input_fmt,
                    &self.datatype,
                );
                self.convert_multiple_files(&files);
                log::info!("{:18}: {}", "Output dir", self.output.display());
            }
            InputType::Wildcard => {
                let files = self.parse_input_wcard(self.matches);
                self.print_input_multi::<PathBuf>(
                    &None,
                    task_desc,
                    files.len(),
                    &self.input_fmt,
                    &self.datatype,
                );
                self.convert_multiple_files(&files);
                log::info!("{:18}: {}", "Output dir", self.output.display());
            }
        }
    }

    fn convert_file(&self, task_desc: &str) {
        let input = Path::new(self.parse_file_input(self.matches));
        let output = self.create_output_fname(&self.output, &self.output_fmt);
        self.check_output_file_exist(&output);
        self.print_input_file(input, task_desc, &self.input_fmt, &self.datatype);
        self.convert_any(input, &output, &self.output_fmt);
    }

    fn convert_multiple_files(&mut self, files: &[PathBuf]) {
        self.is_dir = true;
        self.check_output_dir_exist(&self.output);
        let spin = utils::set_spinner();
        spin.set_message("Converting sequence format...");
        files.par_iter().for_each(|file| {
            let output = self.output.join(file.file_stem().unwrap());
            self.convert_any(file, &output, &self.output_fmt);
        });
        spin.finish_with_message("Finished converting sequence format!\n");
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
}
