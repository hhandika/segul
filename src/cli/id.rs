use std::fs::File;
use std::path::PathBuf;

use ansi_term::Colour::Yellow;
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
    fn parse_output(&self, matches: &ArgMatches) -> PathBuf {
        if matches.is_present("output") {
            let output = self.parse_output(matches);
            output.with_extension("txt")
        } else {
            let input = self.parse_dir_input(matches);
            input.with_extension("txt")
        }
    }
}

impl InputPrint for IdParser<'_> {}

pub(in crate::cli) struct IdParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_dir: Option<PathBuf>,
}

impl<'a> IdParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches) -> Self {
        Self {
            matches,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn get_id(&mut self) {
        let input_fmt = self.parse_input_fmt(self.matches);
        let datatype = self.parse_datatype(self.matches);
        let task_desc = "IDs finding";
        let files = if self.is_input_dir() {
            let dir = self.parse_dir_input(self.matches);
            self.input_dir = Some(PathBuf::from(dir));
            self.get_files(dir, &input_fmt)
        } else {
            self.parse_input_wcard(self.matches)
        };
        self.print_input_multi(
            &self.input_dir,
            task_desc,
            files.len(),
            &input_fmt,
            &datatype,
        );
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
        let fname = self.parse_output(self.matches);
        let file = File::create(&fname).expect("CANNOT CREATE AN OUTPUT FILE");
        let mut writer = BufWriter::new(file);
        ids.iter().for_each(|id| {
            writeln!(writer, "{}", id).unwrap();
        });
        writer.flush().unwrap();
        self.print_output(&fname, ids.len());
    }

    fn print_output(&self, output: &Path, ids: usize) {
        log::info!("\n{}", Yellow.paint("Output"));
        log::info!("{:18}: {}", "Total unique IDs", ids);
        log::info!("{:18}: {}", "File output", output.display());
    }
}
