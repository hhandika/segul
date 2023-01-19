use std::path::PathBuf;

use colored::Colorize;

use crate::cli::AlignSeqPrint;
use crate::handler::translate::Translate;
use crate::helper::types::GeneticCodes;

use super::args::SequenceTranslateArgs;
use super::{collect_paths, AlignSeqInput, InputCli, InputPrint, OutputCli};

impl InputCli for TranslateParser<'_> {}
impl InputPrint for TranslateParser<'_> {}
impl OutputCli for TranslateParser<'_> {}
impl AlignSeqInput for TranslateParser<'_> {}

macro_rules! parse_table_args {
    ($self:ident, $code:ident, $table:ident) => {{
        $self.trans_table = GeneticCodes::$code;
        log::info!("{:18}: {}", "Translation Table", $table);
    }};
}

pub(in crate::cli) struct TranslateParser<'a> {
    args: &'a SequenceTranslateArgs,
    input_dir: Option<PathBuf>,
    trans_table: GeneticCodes,
}

impl<'a> TranslateParser<'a> {
    pub(in crate::cli) fn new(args: &'a SequenceTranslateArgs) -> Self {
        Self {
            args,
            input_dir: None,
            trans_table: GeneticCodes::StandardCode,
        }
    }

    pub(in crate::cli) fn translate(&mut self) {
        if self.args.show_tables {
            self.show_ncbi_tables();
        } else {
            self.translate_all();
        }
    }

    fn translate_all(&mut self) {
        let input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let output_fmt = self.parse_output_fmt(&self.args.out_fmt.output_fmt);
        let frame = self.get_reading_frame();
        let task_desc = "Sequence Translation";
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);
        AlignSeqPrint::new(
            &self.input_dir,
            &input_fmt,
            &datatype,
            task_desc,
            files.len(),
        )
        .print();
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        log::info!("{}", "Params".yellow());
        self.parse_trans_table();
        let translate = Translate::new(&self.trans_table, &input_fmt, &datatype, &output_fmt);
        match frame {
            Some(num) => {
                log::info!("{:18}: {}\n", "Reading frame", &num);
                translate.translate_all(&files, num, &self.args.output);
            }
            None => {
                log::info!("{:18}: Auto\n", "Reading frame");
                translate.translate_all_autoframe(&files, &self.args.output);
            }
        }
    }

    fn parse_trans_table(&mut self) {
        let table = &self.args.table;
        match table {
            1 => log::info!("{:18}: {}", "Translation Table", table),
            2 => parse_table_args!(self, VertMtDna, table),
            3 => parse_table_args!(self, YeastMtDna, table),
            4 => parse_table_args!(self, MoldProtCoelMtDna, table),
            5 => parse_table_args!(self, InvertMtDna, table),
            6 => parse_table_args!(self, CilDasHexNu, table),
            9 => parse_table_args!(self, EchiFlatwormMtDna, table),
            10 => parse_table_args!(self, EuplotidNu, table),
            11 => parse_table_args!(self, BacArchPlantPlast, table),
            12 => parse_table_args!(self, AltYeastNu, table),
            13 => parse_table_args!(self, AsciMtDna, table),
            14 => parse_table_args!(self, AltFlatwormMtDna, table),
            16 => parse_table_args!(self, ChlorMtDna, table),
            21 => parse_table_args!(self, TrematodeMtDna, table),
            22 => parse_table_args!(self, ScenedesmusMtDna, table),
            23 => parse_table_args!(self, ThrausMtDna, table),
            24 => parse_table_args!(self, RhabdopMtDna, table),
            25 => parse_table_args!(self, CaDivSR1GraciBac, table),
            26 => parse_table_args!(self, PachyNu, table),
            29 => parse_table_args!(self, MesodiniumNu, table),
            30 => parse_table_args!(self, PeritrichNu, table),
            33 => parse_table_args!(self, CephalodiscidaeMtDna, table),
            _ => unimplemented!("The Genetic Codes is not supported!"),
        }
    }

    fn get_reading_frame(&self) -> Option<usize> {
        if self.args.auto_read {
            None
        } else {
            Some(self.args.reading_frame)
        }
    }

    fn show_ncbi_tables(&self) {
        println!("{}", "Supported NCBI Genetic Code Tables".yellow());
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
