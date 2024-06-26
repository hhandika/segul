use std::path::PathBuf;

use crate::handler::sequence::translate::Translate;
use crate::helper::logger::AlignSeqLogger;
use crate::helper::types::GeneticCodes;
use colored::Colorize;

use super::args::SequenceTranslateArgs;
use super::{collect_paths, AlignSeqInput, InputCli, OutputCli};

impl InputCli for TranslateParser<'_> {}
impl OutputCli for TranslateParser<'_> {}
impl AlignSeqInput for TranslateParser<'_> {}

pub(in crate::cli) struct TranslateParser<'a> {
    args: &'a SequenceTranslateArgs,
    input_dir: Option<PathBuf>,
}

impl<'a> TranslateParser<'a> {
    pub(in crate::cli) fn new(args: &'a SequenceTranslateArgs) -> Self {
        Self {
            args,
            input_dir: None,
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
        let task = "Sequence Translation";
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);
        AlignSeqLogger::new(
            self.input_dir.as_deref(),
            &input_fmt,
            &datatype,
            files.len(),
        )
        .log(task);
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        log::info!("{}", "Params".yellow());
        self.show_trans_table();
        let table = self.parse_table_num();
        let translate = Translate::new(&input_fmt, &table, &datatype, &output_fmt);
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

    fn parse_table_num(&self) -> GeneticCodes {
        self.args.table.parse().expect("Invalid table number")
    }

    fn show_trans_table(&mut self) {
        log::info!("{:18}: {}", "Translation Table", self.args.table);
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
