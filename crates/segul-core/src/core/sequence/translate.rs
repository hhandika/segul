use std::fs;
use std::path::{Path, PathBuf};

use colored::Colorize;
use rayon::prelude::*;

use crate::core::OutputPrint;
use crate::helper::sequence::{SeqCheck, SeqParser};
use crate::helper::translation::{CodonTable, codon_table, translate_bytes};
use crate::helper::types::{DataType, GeneticCodes, Header, InputFmt, OutputFmt, SeqMatrix};
use crate::helper::{files, utils};
use crate::report::spec::{SequenceTranslateParams, TaskKind, canonical};
use crate::report::stats::TaskCounts;
use crate::report::{DirSnapshot, ReportContext, RunReport, RunReportBuilder, report_path};
use crate::writer::sequences::SeqWriter;

impl OutputPrint for SequenceTranslation<'_> {}

pub struct SequenceTranslation<'a> {
    input_fmt: &'a InputFmt,
    trans_table: &'a GeneticCodes,
    datatype: &'a DataType,
    output_fmt: &'a OutputFmt,
    report_context: ReportContext,
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
            report_context: ReportContext::default(),
        }
    }

    /// Attach frontend provenance and SHA-256 checksum settings for `report.json`.
    pub fn with_report_context(mut self, context: ReportContext) -> Self {
        self.report_context = context;
        self
    }

    /// Translate sequences using a fixed reading frame and write `report.json`.
    pub fn translate_with_report(
        &self,
        files: &[PathBuf],
        frame: usize,
        output: &Path,
    ) -> anyhow::Result<RunReport> {
        let mut sorted = files.to_vec();
        alphanumeric_sort::sort_path_slice(&mut sorted);
        let params = SequenceTranslateParams {
            input_format: canonical::input_fmt(self.input_fmt).to_string(),
            datatype: canonical::datatype(self.datatype).to_string(),
            output_format: canonical::output_fmt(self.output_fmt).to_string(),
            table: canonical::genetic_code(self.trans_table).to_string(),
            reading_frame: frame,
            auto_frame: false,
        };
        let snapshot = DirSnapshot::capture(output);
        let builder = RunReportBuilder::new(TaskKind::SequenceTranslate, &self.report_context)
            .primary_inputs(&sorted)
            .parameters(&params)?;
        self.translate(&sorted, frame, output);
        let written = snapshot.written_files(output);
        let counts = TaskCounts {
            files_read: sorted.len(),
            files_written: written.len(),
            items_processed: Some(written.len()),
            items_skipped: None,
        };
        let report_file = report_path(output, None);
        builder
            .outputs(output, written.into_iter().map(|p| ("sequences", p)))
            .statistics(&counts)?
            .write_to(&report_file)
    }

    /// Translate sequences with automatic frame detection and write `report.json`.
    pub fn translate_all_autoframe_with_report(
        &self,
        files: &[PathBuf],
        output: &Path,
    ) -> anyhow::Result<RunReport> {
        let mut sorted = files.to_vec();
        alphanumeric_sort::sort_path_slice(&mut sorted);
        let params = SequenceTranslateParams {
            input_format: canonical::input_fmt(self.input_fmt).to_string(),
            datatype: canonical::datatype(self.datatype).to_string(),
            output_format: canonical::output_fmt(self.output_fmt).to_string(),
            table: canonical::genetic_code(self.trans_table).to_string(),
            reading_frame: 1,
            auto_frame: true,
        };
        let snapshot = DirSnapshot::capture(output);
        let builder = RunReportBuilder::new(TaskKind::SequenceTranslate, &self.report_context)
            .primary_inputs(&sorted)
            .parameters(&params)?;
        self.translate_all_autoframe(&sorted, output);
        let written = snapshot.written_files(output);
        let counts = TaskCounts {
            files_read: sorted.len(),
            files_written: written.len(),
            items_processed: Some(written.len()),
            items_skipped: None,
        };
        let report_file = report_path(output, None);
        builder
            .outputs(output, written.into_iter().map(|p| ("sequences", p)))
            .statistics(&counts)?
            .write_to(&report_file)
    }

    pub fn translate(&self, files: &[PathBuf], frame: usize, output: &Path) {
        let spin = utils::set_spinner();
        spin.set_message("Translating dna sequences...");
        fs::create_dir_all(output).expect("Failed creating an output directory");
        let table = codon_table(*self.trans_table);
        files.par_iter().for_each(|file| {
            let (mut seq, _) = SeqParser::new(file, self.datatype).parse(self.input_fmt);
            let (trans_mat, header) = self.translate_matrix(&mut seq, frame, &table);
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
        let table = codon_table(*self.trans_table);
        files.par_iter().for_each(|file| {
            let (mut seq, _) = SeqParser::new(file, self.datatype).parse(self.input_fmt);
            let mut frame = 1;
            self.get_reading_frame(file, &seq, &mut frame, &table);
            let (trans_mat, header) = self.translate_matrix(&mut seq, frame, &table);
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

    fn get_reading_frame(
        &self,
        file: &Path,
        matrix: &SeqMatrix,
        frame: &mut usize,
        table: &CodonTable,
    ) {
        let seq = matrix
            .values()
            .next()
            .expect("Failed getting the first sequence");
        let has_internal_stop = self.has_internal_stop(seq, *frame, table);
        if has_internal_stop && *frame < 3 {
            *frame += 1;
            self.get_reading_frame(file, matrix, frame, table);
        } else if has_internal_stop && *frame == 3 {
            panic!(
                "The alignment {} still contains stop codons \
            after testing all possible reading frames",
                file.display()
            )
        }
    }

    // A complete CDS ends in a stop, so a stop is internal only if a full codon follows it.
    // We count bases because X in the translation can be a partial, missing, or ambiguous codon.
    fn has_internal_stop(&self, seq: &[u8], frame: usize, table: &CodonTable) -> bool {
        let mut translated = Vec::with_capacity(seq.len().div_ceil(3));
        translate_bytes(seq, frame, table, &mut translated);
        translated
            .iter()
            .position(|&aa| aa == b'*')
            .is_some_and(|pos| {
                let start = frame - 1 + (pos + 1) * 3;
                let bases_after_stop = seq
                    .get(start..)
                    .unwrap_or(&[])
                    .iter()
                    .filter(|&&b| !matches!(b, b'-' | b'?' | b'.' | b'N' | b'n'))
                    .count();
                bases_after_stop >= 3
            })
    }

    fn translate_matrix(
        &self,
        matrix: &mut SeqMatrix,
        frame: usize,
        table: &CodonTable,
    ) -> (SeqMatrix, Header) {
        let mut trans_matrix =
            SeqMatrix::with_capacity(matrix.len(), matrix.as_bytes().len() / 3 + matrix.len());
        for (id, seq) in matrix.iter() {
            trans_matrix.push_with(id.clone(), |buf| {
                translate_bytes(seq, frame, table, buf);
            });
        }
        matrix.clear();
        let header = self.get_header(&trans_matrix);
        (trans_matrix, header)
    }

    #[cfg(test)]
    fn translate_seq(&self, seq: &str, frame: usize) -> String {
        let table = codon_table(*self.trans_table);
        let mut out = Vec::new();
        translate_bytes(seq.as_bytes(), frame, &table, &mut out);
        String::from_utf8(out).expect("Translation output is ASCII")
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

    macro_rules! test_reading_frame {
        ($input:expr_2021, $result:expr_2021, $code:ident) => {
            let trans = SequenceTranslation::new(
                &InputFmt::Fasta,
                &GeneticCodes::$code,
                &DataType::Dna,
                &OutputFmt::Fasta,
            );
            let mut matrix = SeqMatrix::new();
            matrix.insert(String::from("seq"), $input.as_bytes());
            let mut frame = 1;
            let table = codon_table(GeneticCodes::$code);
            trans.get_reading_frame(Path::new("cds.fas"), &matrix, &mut frame, &table);
            assert_eq!($result, frame, "{}", $input);
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

    #[test]
    fn test_reading_frame_complete_cds() {
        // Frame 2 has no stop codons, so it must not win over the complete CDS in frame 1.
        test_reading_frame!("ATGGCCTAA", 1, StandardCode);
    }

    #[test]
    fn test_reading_frame_complete_cds_stops_in_other_frames() {
        test_reading_frame!("ATGATAGTTAGCTAA", 1, StandardCode);
    }

    #[test]
    fn test_reading_frame_stop_with_partial_codon() {
        test_reading_frame!("ATGGCCTAAG", 1, StandardCode);
        test_reading_frame!("ATGGCCTAAGC", 1, StandardCode);
    }

    #[test]
    fn test_reading_frame_stop_with_gaps() {
        test_reading_frame!("ATGGCCTAA------", 1, StandardCode);
        test_reading_frame!("ATGGCCTAA--", 1, StandardCode);
        test_reading_frame!("ATGGCCTAAG-----", 1, StandardCode);
        test_reading_frame!("ATGGCCTAA--G--C", 1, StandardCode);
    }

    #[test]
    fn test_reading_frame_stop_with_missing() {
        test_reading_frame!("ATGGCCTAA??????", 1, StandardCode);
        test_reading_frame!("ATGGCCTAANNNNNN", 1, StandardCode);
        test_reading_frame!("ATGGCCTAA......", 1, StandardCode);
        test_reading_frame!("atggcctaannn---", 1, StandardCode);
    }

    #[test]
    fn test_reading_frame_complete_cds_offset() {
        test_reading_frame!("TAATGATAGTTAGCTAA", 3, StandardCode);
        test_reading_frame!("TATGATAGTTAGCTAA", 2, StandardCode);
    }

    #[test]
    fn test_reading_frame_complete_cds_vertmtdna() {
        test_reading_frame!("ATGGCCAGA", 1, VertMtDna);
    }

    #[test]
    fn test_reading_frame_internal_stop() {
        test_reading_frame!("ATGTAAGCC", 2, StandardCode);
        test_reading_frame!("ATGGCCTAATAA", 2, StandardCode);
        test_reading_frame!("ATGTAA---GCC", 2, StandardCode);
        test_reading_frame!("ATGGCCTAARAC", 2, StandardCode);
    }

    #[test]
    #[should_panic(expected = "still contains stop codons")]
    fn test_reading_frame_internal_stop_all_frames() {
        test_reading_frame!("TAAATAAATAAATAA", 1, StandardCode);
    }
}
