use std::fs;
use std::path::{Path, PathBuf};

use ansi_term::Colour::Yellow;
use indexmap::IndexMap;
use rayon::prelude::*;

use crate::helper::sequence::{SeqCheck, Sequence};
use crate::helper::translation::NcbiTables;
use crate::helper::types::{
    DataType, GeneticCodes, Header, InputFmt, OutputFmt, PartitionFmt, SeqMatrix,
};
use crate::helper::utils;
use crate::writer::sequences::SeqWriter;

pub struct Translate<'a> {
    input_fmt: &'a InputFmt,
    trans_table: &'a GeneticCodes,
    datatype: &'a DataType,
    frame: usize,
}

impl<'a> Translate<'a> {
    pub fn new(
        trans_table: &'a GeneticCodes,
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        frame: usize,
    ) -> Self {
        Self {
            trans_table,
            input_fmt,
            datatype,
            frame,
        }
    }

    pub fn translate_all(&self, files: &[PathBuf], output: &Path, output_fmt: &OutputFmt) {
        let spin = utils::set_spinner();
        spin.set_message("Translating dna sequences...");
        fs::create_dir_all(output).expect("Failed creating an output directory");
        files.par_iter().for_each(|file| {
            let (mut seq, _) = Sequence::new(file, self.datatype).get(self.input_fmt);
            let trans_mat = self.translate_matrix(&mut seq);
            let header = self.get_header(&trans_mat);
            let outname = self.get_output_names(output, file, output_fmt);
            let mut writer =
                SeqWriter::new(&outname, &trans_mat, &header, None, &PartitionFmt::None);
            writer
                .write_sequence(output_fmt)
                .expect("Failed writing the output files");
        });

        spin.finish_with_message("Finished translating dna sequences!\n");
        log::info!("{}", Yellow.paint("Output"));
        log::info!("{:18}: {}", "Output dir", output.display());
    }

    fn translate_matrix(&self, matrix: &mut SeqMatrix) -> SeqMatrix {
        let mut trans_matrix: SeqMatrix = IndexMap::new();
        matrix.iter().for_each(|(id, seq)| {
            let sequences = self.translate_seq(seq);
            trans_matrix.insert(id.to_string(), sequences);
        });

        matrix.clear();
        trans_matrix
    }

    fn translate_seq(&self, seq: &str) -> String {
        let table = match self.trans_table {
            GeneticCodes::StandardCode => NcbiTables::new().standard_code(),
            GeneticCodes::VertMtDna => NcbiTables::new().vert_mtdna(),
            GeneticCodes::YeastMtDna => NcbiTables::new().yeast_mtdna(),
            GeneticCodes::MoldProtCoelMtDna => NcbiTables::new().moldprotocoe_mtdna(),
            _ => unimplemented!(),
        };
        let mut translation = String::new();
        seq.to_uppercase()
            .chars()
            .skip(self.frame - 1)
            .collect::<Vec<char>>()
            .chunks(3)
            .map(|c| c.iter().collect::<String>())
            .for_each(|c| {
                let aa = table.get(&c);
                match aa {
                    Some(c) => translation.push_str(c),
                    None => translation.push('X'),
                }
            });

        translation
    }

    fn get_output_names(&self, dir: &Path, file: &Path, output_fmt: &OutputFmt) -> PathBuf {
        let path = dir.join(
            file.file_name()
                .expect("Failed parsing filename for output file"),
        );
        match output_fmt {
            OutputFmt::Fasta | OutputFmt::FastaInt => path.with_extension("fas"),
            OutputFmt::Nexus | OutputFmt::NexusInt => path.with_extension("nex"),
            OutputFmt::Phylip | OutputFmt::PhylipInt => path.with_extension("phy"),
        }
    }

    fn get_header(&self, matrix: &SeqMatrix) -> Header {
        let mut seq_info = SeqCheck::new();
        seq_info.check(matrix);
        let mut header = Header::new();
        header.aligned = seq_info.is_alignment;
        header.nchar = seq_info.longest;
        header.ntax = matrix.len();
        header.datatype = String::from("protein");
        header
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_translate;

    #[macro_export]
    macro_rules! test_translate {
        ($input:expr, $frame:expr, $result:expr, $code:ident) => {
            let trans = Translate::new(
                &GeneticCodes::$code,
                &InputFmt::Fasta,
                &DataType::Dna,
                $frame,
            );
            assert_eq!($result, trans.translate_seq($input));
        };
    }

    #[test]
    fn test_translation_simple() {
        let dna = "AAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("KGDLVRX"), StandardCode);
    }

    #[test]
    fn test_translating_stdcode() {
        let standard_code = "TTTTTCTTATTGCTTCTCCTA\
        CTGATTATCATAATGGTTGTCGTAGTGTCTTC\
        CTCATCGCCTCCCCCACCGACTACCACAACGGC\
        TGCCGCAGCGTATTACTAATAGCATCACCAACAG\
        AATAACAAAAAGGATGACGAAGAGTGTTGCTGAT\
        GGCGTCGCCGACGGAGTAGCAGAAGGGGTGGCGGAGGG";
        let frame = 1;
        test_translate!(
            standard_code,
            frame,
            String::from(
                "FFLLLLLLIIIMVVVVSSSSPPPPTTTT\
            AAAAYY**HHQQNNKKDDEECC*WRRRRSSRRGGGG"
            ),
            StandardCode
        );
    }

    #[test]
    fn test_translating_with_gaps() {
        let dna = "AAAGGGGATTTAGTTAGAA-----";
        let frame = 1;
        test_translate!(dna, frame, String::from("KGDLVRX-"), StandardCode);
    }

    #[test]
    fn test_translation_vertmtdna_simple() {
        let dna = "AAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("KGDLV*X"), VertMtDna);
    }

    #[test]
    fn test_translation_yestmtdna_simple() {
        let dna = "CTTATAAAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("TMKGDLVRX"), YeastMtDna);
    }

    #[test]
    fn test_translating_mold_etal_simple() {
        let dna = "TGAAAAGGGGATTTAGTTAGAA-----";
        let frame = 1;
        test_translate!(dna, frame, String::from("WKGDLVRX-"), MoldProtCoelMtDna);
    }
}
