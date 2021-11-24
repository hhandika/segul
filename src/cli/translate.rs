use ansi_term::Colour::Yellow;

use crate::cli::*;
use crate::core::translate::Translate;
use crate::helper::types::GeneticCodes;

impl InputCli for TranslateParser<'_> {}
impl InputPrint for TranslateParser<'_> {}
impl OutputCli for TranslateParser<'_> {}

pub(in crate::cli) struct TranslateParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_dir: Option<PathBuf>,
    trans_table: GeneticCodes,
}

#[allow(unused_variables)]
impl<'a> TranslateParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            input_dir: None,
            trans_table: GeneticCodes::StandardCode,
        }
    }

    pub(in crate::cli) fn parse_translation_cli(&mut self) {
        if self.matches.is_present("show-tables") {
            self.show_ncbi_tables();
        } else {
            self.translate();
        }
    }

    fn translate(&mut self) {
        let input_fmt = self.parse_input_fmt(self.matches);
        let datatype = self.parse_datatype(self.matches);
        let output_fmt = self.parse_output_fmt(self.matches);
        let outdir = self.parse_output(self.matches);
        let frame = self.get_reading_frame();
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
        let translate = Translate::new(&self.trans_table, &input_fmt, &datatype, frame);
        translate.translate_sequences(&files, &outdir, &output_fmt);
    }

    fn parse_trans_table(&mut self) {
        match self.matches {
            m if m.is_present("2") => self.trans_table = GeneticCodes::VertMtDna,
            _ => (),
        }
    }

    fn get_reading_frame(&self) -> usize {
        self.matches
            .value_of("reading-frame")
            .expect("Failed getting reading frame values")
            .parse::<usize>()
            .expect("Failed parsing reading frame values")
    }

    fn show_ncbi_tables(&self) {
        println!("{}", Yellow.paint("NCBI Genetic Code Tables"));
        println!(
            "1. The Standard Code\n\
            2. The Vertebrate Mitochondrial Code\n\
            3. The Yeast Mitochondrial Code\n\
            4. TheMold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code\n\
            5. The Invertebrate Mitochondrial Code\n\
            6. The Ciliate, Dasycladacean and Hexamita Nuclear Code\n\
            9. The Echinoderm and Flatworm Mitochondrial Code\n\
            10. The Euplotid Nuclear Code\n\
            11. The Bacterial, Archaeal and Plant Plastid Code\n\
            12. The Alternative Yeast Nuclear Code\n\
            13. The Ascidian Mitochondrial Code\n\
            14. The Alternative Flatworm Mitochondrial Code\n\
            16. Chlorophycea nMitochondrial Code\n\
            21. Trematode Mitochondrial Code\n\
            22. Scenedesmus obliquus Mitochondrial Code\n\
            23. Thraustochytrium Mitochondrial Code\n\
            24. Rhabdopleuridae Mitochondrial Code\n\
            25. Candidate Division SR1 and Gracilibacteria Code\n\
            26. Pachysolen tannophilus Nuclear Code\n\
            27. Karyorelict Nuclear Code\n\
            28. Condylostoma Nuclear Code\n\
            29. Mesodinium Nuclear Code\n\
            30. Peritrich Nuclear Code\n\
            31. Blastocrithidia Nuclear Code\n\
            33. Cephalodiscidae Mitochondrial UAA-Tyr Code\n\
            "
        );
    }
}
