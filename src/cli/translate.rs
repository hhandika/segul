use ansi_term::Colour::Yellow;

use crate::cli::*;
use crate::core::translate::Translate;
use crate::helper::types::GeneticCodes;
use crate::parse_table_args;

impl InputCli for TranslateParser<'_> {}
impl InputPrint for TranslateParser<'_> {}
impl OutputCli for TranslateParser<'_> {}

pub(in crate::cli) struct TranslateParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_dir: Option<PathBuf>,
    trans_table: GeneticCodes,
}

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
        log::info!("{}", Yellow.paint("Params"));
        self.parse_trans_table();
        let translate = Translate::new(&self.trans_table, &input_fmt, &datatype);
        match frame {
            Some(num) => {
                log::info!("{:18}: {}\n", "Reading frame", &num);
                translate.translate_all(&files, num, &outdir, &output_fmt);
            }
            None => {
                log::info!("{:18}: Auto\n", "Reading frame");
                translate.translate_all_autoframe(&files, &outdir, &output_fmt);
            }
        }
    }

    fn parse_trans_table(&mut self) {
        let table = self
            .matches
            .value_of("table")
            .expect("Failed parsing table input");
        match table {
            "1" => log::info!("{:18}: {}", "Translation Table", table),
            "2" => parse_table_args!(self, VertMtDna, table),
            "3" => parse_table_args!(self, YeastMtDna, table),
            "4" => parse_table_args!(self, MoldProtCoelMtDna, table),
            "5" => parse_table_args!(self, InvertMtDna, table),
            "6" => parse_table_args!(self, CilDasHexNu, table),
            "9" => parse_table_args!(self, EchiFlatwormMtDna, table),
            "10" => parse_table_args!(self, EuplotidNu, table),
            "11" => parse_table_args!(self, BacArchPlantPlast, table),
            "12" => parse_table_args!(self, AltYeastNu, table),
            "13" => parse_table_args!(self, AsciMtDna, table),
            "14" => parse_table_args!(self, AltFlatwormMtDna, table),
            "16" => parse_table_args!(self, ChlorMtDna, table),
            "21" => parse_table_args!(self, TrematodeMtDna, table),
            "22" => parse_table_args!(self, ScenedesmusMtDna, table),
            "23" => parse_table_args!(self, ThrausMtDna, table),
            "24" => parse_table_args!(self, RhabdopMtDna, table),
            "25" => parse_table_args!(self, CaDivSR1GraciBac, table),
            "26" => parse_table_args!(self, PachyNu, table),
            "29" => parse_table_args!(self, MesodiniumNu, table),
            "30" => parse_table_args!(self, PeritrichNu, table),
            "33" => parse_table_args!(self, CephalodiscidaeMtDna, table),
            _ => unimplemented!("The Genetic Codes is not supported!"),
        }
    }

    fn get_reading_frame(&self) -> Option<usize> {
        let frame = self
            .matches
            .value_of("reading-frame")
            .expect("Failed getting reading frame values");
        match frame {
            "auto" => None,
            _ => Some(
                frame
                    .parse::<usize>()
                    .expect("Failed parsing reading frame values"),
            ),
        }
    }

    fn show_ncbi_tables(&self) {
        println!("{}", Yellow.paint("Supported NCBI Genetic Code Tables"));
        println!(
            "1. The Standard Code\n\
            2. The Vertebrate Mitochondrial Code\n\
            3. The Yeast Mitochondrial Code\n\
            4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code\n\
            5. The Invertebrate Mitochondrial Code\n\
            6. The Ciliate, Dasycladacean and Hexamita Nuclear Code\n\
            9. The Echinoderm and Flatworm Mitochondrial Code\n\
            10. The Euplotid Nuclear Code\n\
            11. The Bacterial, Archaeal and Plant Plastid Code\n\
            12. The Alternative Yeast Nuclear Code\n\
            13. The Ascidian Mitochondrial Code\n\
            14. The Alternative Flatworm Mitochondrial Code\n\
            16. Chlorophycean Mitochondrial Code\n\
            21. Trematode Mitochondrial Code\n\
            22. Scenedesmus obliquus Mitochondrial Code\n\
            23. Thraustochytrium Mitochondrial Code\n\
            24. Rhabdopleuridae Mitochondrial Code\n\
            25. Candidate Division SR1 and Gracilibacteria Code\n\
            26. Pachysolen tannophilus Nuclear Code\n\
            29. Mesodinium Nuclear Code\n\
            30. Peritrich Nuclear Code\n\
            33. Cephalodiscidae Mitochondrial UAA-Tyr Code\n
            "
        );
    }
}

#[macro_export]
macro_rules! parse_table_args {
    ($self:ident, $code:ident, $table:ident) => {{
        $self.trans_table = GeneticCodes::$code;
        log::info!("{:18}: {}", "Translation Table", $table);
    }};
}
