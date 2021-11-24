use crate::cli::*;
use crate::core::translate::{NCBITable, Translate};

impl InputCli for TranslateParser<'_> {}
impl InputPrint for TranslateParser<'_> {}
impl OutputCli for TranslateParser<'_> {}

pub(in crate::cli) struct TranslateParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_dir: Option<PathBuf>,
    trans_table: NCBITable,
}

#[allow(unused_variables)]
impl<'a> TranslateParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            input_dir: None,
            trans_table: NCBITable::StandardCode,
        }
    }

    pub(in crate::cli) fn translate(&mut self) {
        let input_fmt = self.parse_input_fmt(self.matches);
        let datatype = self.parse_datatype(self.matches);
        let output_fmt = self.parse_output_fmt(self.matches);
        let outdir = self.parse_output(self.matches);
        let task_desc = "Sequence Translation";
        let files = if self.matches.is_present("wildcard") {
            self.parse_input_wcard(self.matches)
        } else {
            let dir = self.parse_dir_input(self.matches);
            self.input_dir = Some(PathBuf::from(dir));
            self.get_files(dir, &input_fmt)
        };

        self.print_input_multi(
            &self.input_dir,
            task_desc,
            files.len(),
            &input_fmt,
            &datatype,
        );

        self.check_output_dir_exist(&outdir);
        self.parse_trans_table();
        let translate = Translate::new(&self.trans_table, &input_fmt, &datatype);
        translate.translate_sequences(&files, &outdir, &output_fmt);
    }

    fn parse_trans_table(&mut self) {
        match self.matches {
            m if m.is_present("2") => self.trans_table = NCBITable::MtDna,
            _ => (),
        }
    }
}
