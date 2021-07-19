use std::fs::File;

use clap::ArgMatches;

use crate::cli::*;

impl InputCli for IdParser<'_> {
    fn parse_input_type(&self, matches: &ArgMatches) -> InputType {
        if matches.is_present("dir") {
            InputType::Dir
        } else {
            InputType::Wildcard
        }
    }
}

impl OutputCli for IdParser<'_> {
    fn parse_output_path(&self, matches: &ArgMatches) -> PathBuf {
        if matches.is_present("output") {
            let output = self.parse_output(matches);
            PathBuf::from(output).with_extension("txt")
        } else {
            let input = Path::new(self.parse_dir_input(matches));
            input.with_extension("txt")
        }
    }
}

pub(in crate::cli) struct IdParser<'a> {
    matches: &'a ArgMatches<'a>,
}

impl<'a> IdParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches) -> Self {
        Self { matches }
    }

    pub(in crate::cli) fn get_id(&self) {
        let input_fmt = self.parse_input_fmt(&self.matches);
        let datatype = self.parse_datatype(self.matches);
        let files = if self.is_input_dir() {
            let dir = self.parse_dir_input(self.matches);
            self.get_files(dir, &input_fmt)
        } else {
            self.parse_input_wcard(&self.matches)
        };
        self.print_input().unwrap();
        let spin = utils::set_spinner();
        spin.set_message("Indexing IDs..");
        let ids = IDs::new(&files, &input_fmt, &datatype).get_id_all();
        spin.finish_with_message("DONE!");
        self.write_results(&ids);
    }

    fn is_input_dir(&self) -> bool {
        self.matches.is_present("dir")
    }

    fn write_results(&self, ids: &IndexSet<String>) {
        let fname = self.parse_output_path(&self.matches);
        let file = File::create(&fname).expect("CANNOT CREATE AN OUTPUT FILE");
        let mut writer = BufWriter::new(file);
        ids.iter().for_each(|id| {
            writeln!(writer, "{}", id).unwrap();
        });
        writer.flush().unwrap();
        self.print_output(&fname, ids.len()).unwrap();
    }

    fn print_input(&self) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "Command\t\t\t: segul id")?;
        if self.is_input_dir() {
            writeln!(
                writer,
                "Input dir\t: {}\n",
                self.parse_dir_input(self.matches)
            )?;
        } else {
            writeln!(writer, "Input\t\t: WILDCARD",)?;
        }
        Ok(())
    }

    fn print_output(&self, output: &Path, ids: usize) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "\nTotal unique IDs\t: {}", ids)?;
        writeln!(writer, "Output\t\t\t: {}", output.display())?;

        Ok(())
    }
}
