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

impl OutputCli for IdParser<'_> {}

impl InputPrint for IdParser<'_> {}

pub(in crate::cli) struct IdParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_dir: Option<PathBuf>,
    output: PathBuf,
}

impl<'a> IdParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches) -> Self {
        Self {
            matches,
            input_dir: None,
            output: PathBuf::new(),
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
        self.output = self.parse_output(self.matches).with_extension("txt");
        self.check_output_file_exist(&self.output);
        let spin = utils::set_spinner();
        spin.set_message("Indexing IDs..");
        let ids = IDs::new(&files, &input_fmt, &datatype).get_id_all();
        spin.finish_with_message("DONE!");
        self.write_results(&ids).expect("Failed writing results");
        self.print_output(ids.len());
    }

    fn is_input_dir(&self) -> bool {
        self.matches.is_present("dir")
    }

    fn write_results(&self, ids: &IndexSet<String>) -> Result<()> {
        let file = File::create(&self.output).expect("CANNOT CREATE AN OUTPUT FILE");
        let mut writer = BufWriter::new(file);
        ids.iter().for_each(|id| {
            writeln!(writer, "{}", id).unwrap();
        });
        writer.flush()?;
        Ok(())
    }

    fn print_output(&self, ids: usize) {
        log::info!("\n{}", Yellow.paint("Output"));
        log::info!("{:18}: {}", "Total unique IDs", ids);
        log::info!("{:18}: {}", "File output", self.output.display());
    }
}

// #[cfg(test)]
// mod test {
//     use super::*;
//     use clap::{App, Arg};

//     #[test]
//     fn get_id_output_path_test() {
//         let arg = App::new("segul-test")
//             .arg(Arg::with_name("dir").default_value("./test_dir/"))
//             .get_matches();
//         let id = IdParser::new(&arg);
//         let res = PathBuf::from("./test_dir.txt");
//         assert_eq!(res, id.output);
//     }
// }
