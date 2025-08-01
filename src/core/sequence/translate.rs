use std::fs;
use std::path::{Path, PathBuf};

use ahash::AHashMap as HashMap;
use colored::Colorize;
use indexmap::IndexMap;
use rayon::prelude::*;

use crate::core::OutputPrint;
use crate::helper::sequence::{SeqCheck, SeqParser};
use crate::helper::translation::NcbiTables;
use crate::helper::types::{DataType, GeneticCodes, Header, InputFmt, OutputFmt, SeqMatrix};
use crate::helper::{files, utils};
use crate::writer::sequences::SeqWriter;

impl OutputPrint for SequenceTranslation<'_> {}

pub struct SequenceTranslation<'a> {
    input_fmt: &'a InputFmt,
    trans_table: &'a GeneticCodes,
    datatype: &'a DataType,
    output_fmt: &'a OutputFmt,
}

impl<'a> SequenceTranslation<'a> {
    pub fn new(
        input_fmt: &'a InputFmt,
        trans_table: &'a GeneticCodes,
        datatype: &'a DataType,
        output_fmt: &'a OutputFmt,
    ) -> Self {
        Self {
            trans_table,
            input_fmt,
            datatype,
            output_fmt,
        }
    }

    pub fn translate(&self, files: &[PathBuf], frame: usize, output: &Path) {
        let spin = utils::set_spinner();
        spin.set_message("Translating dna sequences...");
        fs::create_dir_all(output).expect("Failed creating an output directory");
        files.par_iter().for_each(|file| {
            let (mut seq, _) = SeqParser::new(file, self.datatype).parse(self.input_fmt);
            let (trans_mat, header) = self.translate_matrix(&mut seq, frame);
            let output_fname = files::create_output_fname(output, file, self.output_fmt);
            let mut writer = SeqWriter::new(&output_fname, &trans_mat, &header);
            writer
                .write_sequence(self.output_fmt)
                .expect("Failed writing the output files");
        });

        spin.finish_with_message("Finished translating dna sequences!\n");
        self.print_output_info(output);
    }

    pub fn translate_all_autoframe(&self, files: &[PathBuf], output: &Path) {
        let spin = utils::set_spinner();
        spin.set_message("Translating dna sequences...");
        files.par_iter().for_each(|file| {
            let (mut seq, _) = SeqParser::new(file, self.datatype).parse(self.input_fmt);
            let mut frame = 1;
            self.get_reading_frame(file, &seq, &mut frame);
            let (trans_mat, header) = self.translate_matrix(&mut seq, frame);
            let output_dir = output.join(format!("RF-{}", frame));
            fs::create_dir_all(output).expect("Failed creating an output directory");
            let outname = files::create_output_fname(&output_dir, file, self.output_fmt);
            let mut writer = SeqWriter::new(&outname, &trans_mat, &header);
            writer
                .write_sequence(self.output_fmt)
                .expect("Failed writing the output files");
        });

        spin.finish_with_message("Finished translating dna sequences!\n");
        self.print_output_info(output);
    }

    fn get_reading_frame(&self, file: &Path, matrix: &SeqMatrix, frame: &mut usize) {
        let seq = matrix
            .values()
            .next()
            .expect("Failed getting the first sequence");
        let trans = self.translate_seq(seq, *frame);
        if trans.contains('*') && *frame < 3 {
            *frame += 1;
            self.get_reading_frame(file, matrix, frame);
        } else if trans.contains('*') && *frame == 3 {
            panic!(
                "The alignment {} still contains stop codons \
            after testing all possible reading frames",
                file.display()
            )
        }
    }

    fn translate_matrix(&self, matrix: &mut SeqMatrix, frame: usize) -> (SeqMatrix, Header) {
        let mut trans_matrix: SeqMatrix = IndexMap::new();
        matrix.iter().for_each(|(id, seq)| {
            let sequences = self.translate_seq(seq, frame);
            trans_matrix.insert(id.to_string(), sequences);
        });
        matrix.clear();
        let header = self.get_header(&trans_matrix);
        (trans_matrix, header)
    }

    // We iter over the sequence and translate each codon to its corresponding amino acid
    // We collect minus the first frame - 1 characters and then collect the rest in chunks of 3
    // This allow us to adjust the reading frame when it doesn't start from the first character
    fn translate_seq(&self, seq: &str, frame: usize) -> String {
        let table = self.get_ncbi_tables();
        let mut translation = String::new();
        seq.to_uppercase()
            .chars()
            .skip(frame - 1)
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

    /// Match user input with supported genetic codes
    /// and return the corresponding translation table
    fn get_ncbi_tables(&self) -> HashMap<String, String> {
        match self.trans_table {
            GeneticCodes::StandardCode => NcbiTables::new().standard_code(), // Table 1
            GeneticCodes::VertMtDna => NcbiTables::new().vert_mtdna(),       // Table 2
            GeneticCodes::YeastMtDna => NcbiTables::new().yeast_mtdna(),     // Table 3
            GeneticCodes::MoldProtCoelMtDna => NcbiTables::new().moldprotocoe_mtdna(), // Table 4
            GeneticCodes::InvertMtDna => NcbiTables::new().invert_mtdna(),   // Table 5
            GeneticCodes::CilDasHexNu => NcbiTables::new().cildashex_nudna(), // Table 6
            GeneticCodes::EchiFlatwormMtDna => NcbiTables::new().echiflatworm_mtdna(), // Table 9
            GeneticCodes::EuplotidNu => NcbiTables::new().euplotid_nudna(),  // Table 10
            GeneticCodes::BacArchPlantPlast => NcbiTables::new().standard_code(), // Table 11
            GeneticCodes::AltYeastNu => NcbiTables::new().alt_yeast_nu(),    // Table 12
            GeneticCodes::AsciMtDna => NcbiTables::new().ascidian_mtdna(),   // Table 13
            GeneticCodes::AltFlatwormMtDna => NcbiTables::new().alt_flatworm_mtdna(), // Table 14
            GeneticCodes::ChlorMtDna => NcbiTables::new().chlorophycean_mtdna(), // Table 16
            GeneticCodes::TrematodeMtDna => NcbiTables::new().trematode_mtdna(), // Table 21
            GeneticCodes::ScenedesmusMtDna => NcbiTables::new().scenedesmus_mtdna(), // Table 22
            GeneticCodes::ThrausMtDna => NcbiTables::new().thraustochytrium_mtdna(), // Table 23
            GeneticCodes::RhabdopMtDna => NcbiTables::new().rhabdopleuridae_mtdna(), // Table 24
            GeneticCodes::CaDivSR1GraciBac => NcbiTables::new().candid_div_sr1_gracil(), // Table 25
            GeneticCodes::PachyNu => NcbiTables::new().pachysolen_tanno_nu(), // Table 26
            GeneticCodes::MesodiniumNu => NcbiTables::new().mesodinium_nu(), // Table 29
            GeneticCodes::PeritrichNu => NcbiTables::new().peritrich_nu(),   // Table 30
            GeneticCodes::CephalodiscidaeMtDna => NcbiTables::new().cephalodiscidae_mtdna(), // Table 33
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

    fn print_output_info(&self, output: &Path) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Output dir", output.display());
        self.print_output_fmt(self.output_fmt);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! test_translate {
        ($input:expr_2021, $frame:expr_2021, $result:expr_2021, $code:ident) => {
            let trans = SequenceTranslation::new(
                &InputFmt::Fasta,
                &GeneticCodes::$code,
                &DataType::Dna,
                &OutputFmt::Fasta,
            );
            assert_eq!($result, trans.translate_seq($input, $frame));
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
    fn test_translate_frame_2() {
        let dna = "TAAAGGGGATTTAGTTAGAA";
        let frame = 2;
        test_translate!(dna, frame, String::from("KGDLVRX"), StandardCode);
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
    fn test_translation_yestmtdna() {
        // Dna from https://www.ncbi.nlm.nih.gov/nuccore/NM_001184304.1.
        // Note that the genebank translation is using the standard code.
        let dna = "ATGTCCGGAAAACCACCAGTTTATAGATTACCCCCTCTTCCCAGACTAAAAGTGAAGAAACCTATCATTA\
        GGCAAGAAGCGAACAAATGTCTTGTTTTAATGTCAAACTTATTACAATGCTGGTCTTCATATGGTCATAT\
        GAGCCCTAAATGTGCTGGCTTAGTAACGGAATTGAAAAGTTGCACTAGTGAAAGTGCCCTCGGCAAAAGA\
        AACAACGTTCAAAAGAGCAATATTAACTATCACGCCGCTAGGCTGTATGATAGGATCAACGGTAAACCCC\
        ATGATTAA";
        let frame = 1;
        test_translate!(
            dna,
            frame,
            String::from(
                "MSGKPPVYRLPPTPRTKVKKPIIRQEANKC\
        TVLMSNLLQCWSSYGHMSPKCAGLVTELKSCTSESATGKRNNVQKSNINYHAARTYDRINGKPHD*"
            ),
            YeastMtDna
        );
    }

    #[test]
    fn test_translating_mold_etal_simple() {
        let dna = "TGAAAAGGGGATTTAGTTAGAA-----";
        let frame = 1;
        test_translate!(dna, frame, String::from("WKGDLVRX-"), MoldProtCoelMtDna);
    }

    #[test]
    fn test_translation_invertmtdna_simple() {
        let dna = "AAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("KGDLVSX"), InvertMtDna);
    }

    #[test]
    fn test_invert_mtdna_real_sequence() {
        // Sequence from https://www.ncbi.nlm.nih.gov/nuccore/MT602627.1
        let dna = "TAGGTACATCTTTGAGATGATTAATTCGTACAGAATTAGGTAACCCAGGATCTTTAATTGGAGATGATCA\
        AATTTATAATACTATTGTTACAGCTCATGCTTTTATTATAATTTTTTTTATAGTTATACCTATTATAATT\
        GGAGGATTTGGCAATTGACTAGTACCATTAATATTGGGAGCCCCTGATATAGCTTGTCCTCGAATAAGTA\
        ACATGAGCTTTTGACTTTTACCACCTTCATTAACTTTATTAATTTCTAGAAGAATTGTAGAAAATGGTGC\
        CGGAACTGGCTGAACTGTTTACCCCCCTCTTTCTTCTCATATTGCTCATCATGGAGCTTCTGTTGATTTA\
        GCTATTTTTTCCCTTCACTTAGCAGGAATTTCTTCAATTTTAGGAGCTATTAATTTTATTACTACAATTA\
        TTAATATGCGAATTATAAATTTATCTTTTGACCAATTACCTTTATTTGTATGATCTGGTGGGATTACAGC\
        TTTACTATTATTACTCTCTTTACCTGTCTTAGCAGGAGCTATTACTATACTTTTAACAGATCGAAATTTA\
        AACACATCTGTTTTTGATCCAGCGGGAGGAGGAGATCCTATTTTATATCAACATTTATTATTCGATTT";

        let res = String::from(
            "GTSLSWLIRTELGNPGSLIGDDQIYNTIVTAHAFIMIFFMVMPI\
                    MIGGFGNWLVPLMLGAPDMACPRMSNMSFWLLPPSLTLLISSSIVENGAGTGWTVYPP\
                    LSSHIAHHGASVDLAIFSLHLAGISSILGAINFITTIINMRIMNLSFDQLPLFVWSGG\
                    ITALLLLLSLPVLAGAITMLLTDRNLNTSVFDPAGGGDPILYQHLLFDX",
        );

        let frame = 3;

        test_translate!(dna, frame, res, InvertMtDna);
    }

    #[test]
    fn test_translating_cildashex_simple() {
        let dna = "TAGTAAAAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("QQKGDLVRX"), CilDasHexNu);
    }

    #[test]
    fn test_translating_echiflatworms_simple() {
        let dna = "TGAAGGAAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("WSNGDLVSX"), EchiFlatwormMtDna);
    }

    #[test]
    fn test_translating_euplotid_simple() {
        let dna = "TGAAAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("CKGDLVRX"), EuplotidNu);
    }

    #[test]
    fn test_translating_bacteria_simple() {
        let dna = "AAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("KGDLVRX"), BacArchPlantPlast);
    }

    #[test]
    fn test_translating_altyeast_nu_simple() {
        let dna = "CTGAAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("SKGDLVRX"), AltYeastNu);
    }

    #[test]
    fn test_translating_ascidian_mtdna_simple() {
        let dna = "TGAATAAGGAGAAAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("WMGGKGDLVGX"), AsciMtDna);
    }

    #[test]
    fn test_alt_flatworm_mtdna_simple() {
        let dna = "TAATGAAGGAAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("YWSNGDLVSX"), AltFlatwormMtDna);
    }

    #[test]
    fn test_chlorophycean_mtdna_simple() {
        let dna = "TAGAAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("LKGDLVRX"), ChlorMtDna);
    }

    #[test]
    fn test_trematode_mtdna_simple() {
        let dna = "TGAATAAGGAGAAAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("WMSSNGDLVSX"), TrematodeMtDna);
    }

    #[test]
    fn test_scenedesmus_mtdna_simple() {
        let dna = "TCATAGTGAATAAGGAGAAAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("*L*IRRKGDLVRX"), ScenedesmusMtDna);
    }

    #[test]
    fn test_translating_thrasutochytrium_mtdna_simple() {
        let dna = "TTAAAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("*KGD*VRX"), ThrausMtDna);
    }

    #[test]
    fn test_translating_rhadop_mtdna_simple() {
        let dna = "TGAAGGAGATTAAAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("WKSLKGDLVSX"), RhabdopMtDna);
    }

    #[test]
    fn test_translating_candid_sr1_gracil() {
        let dna = "TGAAAAGGGGATTTAGTTAGAA-----";
        let frame = 1;
        test_translate!(dna, frame, String::from("GKGDLVRX-"), CaDivSR1GraciBac);
    }

    #[test]
    fn test_translating_pachy_simple() {
        let dna = "CTGAAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("AKGDLVRX"), PachyNu);
    }

    #[test]
    fn test_translating_mesodinium_simple() {
        let dna = "TAGTAAAAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("YYKGDLVRX"), MesodiniumNu);
    }

    #[test]
    fn test_translating_peritrich_simple() {
        let dna = "TAGTAACTGAAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(dna, frame, String::from("EELKGDLVRX"), PeritrichNu);
    }

    #[test]
    fn test_translating_cepha_simple() {
        let dna = "TAATGAAGAAGGCTGAAAGGGGATTTAGTTAGAA";
        let frame = 1;
        test_translate!(
            dna,
            frame,
            String::from("YWSKLKGDLVSX"),
            CephalodiscidaeMtDna
        );
    }

    #[test]
    fn test_translating_with_missing() {
        let dna = "NNNAAAGGGGATTTAGTTAGAA?????";
        let frame = 1;
        test_translate!(dna, frame, String::from("XKGDLVRX?"), StandardCode);
    }

    #[test]
    fn test_translating_lowercase() {
        let dna = "nnnaaaggggatttagttagaa?????";
        let frame = 1;
        test_translate!(dna, frame, String::from("XKGDLVRX?"), StandardCode);
    }
}
