use std::fs;
use std::path::{Path, PathBuf};

use indexmap::IndexMap;
use rayon::prelude::*;

use crate::helper::sequence::{SeqCheck, Sequence};
use crate::helper::translation::NcbiTables;
use crate::helper::types::{
    DataType, GeneticCodes, Header, InputFmt, OutputFmt, PartitionFmt, SeqMatrix,
};
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

    pub fn translate_sequences(&self, files: &[PathBuf], output: &Path, output_fmt: &OutputFmt) {
        fs::create_dir_all(output).expect("Failed creating an output directory");
        files.par_iter().for_each(|file| {
            let (mut seq, _) = Sequence::new(file, self.datatype).get(self.input_fmt);
            let trans_mat = self.translate(&mut seq);
            let header = self.get_header(&trans_mat);
            let outname = self.get_output_names(output, file, output_fmt);
            let mut writer =
                SeqWriter::new(&outname, &trans_mat, &header, None, &PartitionFmt::None);
            writer
                .write_sequence(output_fmt)
                .expect("Failed writing the output files");
        });
    }

    fn translate(&self, matrix: &mut SeqMatrix) -> SeqMatrix {
        let mut trans_matrix: SeqMatrix = IndexMap::new();
        matrix.iter().for_each(|(id, seq)| {
            let sequences = self.match_translation(seq);
            trans_matrix.insert(id.to_string(), sequences);
        });

        matrix.clear();
        trans_matrix
    }

    fn match_translation(&self, seq: &str) -> String {
        let table = match self.trans_table {
            GeneticCodes::StandardCode => NcbiTables::new().standard_code(),
            GeneticCodes::VertMtDna => NcbiTables::new().vert_mtdna(),
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
            assert_eq!($result, trans.match_translation($input));
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
    fn test_translation_mtdna_simple() {
        let dna = "AAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("KGDLV*X"), VertMtDna);
    }
}
