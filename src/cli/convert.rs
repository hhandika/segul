use crate::core::converter::Converter;
use clap::ArgMatches;

use crate::cli::*;
use crate::helper::common::{DataType, InputFmt, OutputFmt};

impl Cli for ConvertParser<'_> {}

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
        self.input_fmt = self.get_input_fmt(&self.matches);
        self.output = self.get_output_path(self.matches);
        self.output_fmt = self.get_output_fmt(self.matches);
        self.datatype = self.get_datatype(&self.matches);
        let input_type = self.get_input_type(&self.matches);
        match input_type {
            InputType::File => self.convert_file(),
            InputType::Dir => {
                let dir = self.get_dir_input(self.matches);
                let files = self.get_files(dir, &self.input_fmt);
                self.convert_multiple_files(&files);
                self.print_input_dir(Path::new(dir), files.len(), &self.output)
                    .unwrap();
            }
            InputType::Wildcard => {
                let files = self.parse_input_wcard(&self.matches);
                self.convert_multiple_files(&files)
            }
        }
    }

    fn convert_file(&mut self) {
        let input = Path::new(self.get_file_input(self.matches));
        self.print_input_file(input).unwrap();
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

    fn print_input_file(&self, input: &Path) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "Command\t\t: segul convert")?;
        writeln!(writer, "Input\t\t: {}", &input.display())?;
        Ok(())
    }

    fn print_input_dir(&self, input: &Path, nfile: usize, output: &Path) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "Command\t\t: segul convert")?;
        writeln!(writer, "Input dir\t: {}", &input.display())?;
        writeln!(writer, "Total files\t: {}", utils::fmt_num(&nfile))?;
        writeln!(writer, "Output dir \t: {}\n", output.display())?;
        Ok(())
    }
}
