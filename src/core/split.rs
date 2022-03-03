use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use ansi_term::Colour::Yellow;
use indexmap::IndexMap;
use rayon::prelude::*;

use crate::core::OutputPrint;
use crate::helper::filenames;
use crate::helper::sequence::Sequence;
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, PartitionFmt, SeqMatrix};
use crate::helper::utils;
use crate::parser::partition::PartitionParser;
use crate::writer::sequences::SeqWriter;

impl OutputPrint for Splitter<'_> {}

pub struct Splitter<'a> {
    input: &'a Path,
    datatype: &'a DataType,
    input_fmt: &'a InputFmt,
    output: &'a Path,
    output_fmt: &'a OutputFmt,
}

impl<'a> Splitter<'a> {
    pub fn new(
        input: &'a Path,
        datatype: &'a DataType,
        input_fmt: &'a InputFmt,
        output: &'a Path,
        output_fmt: &'a OutputFmt,
    ) -> Self {
        Self {
            input,
            datatype,
            input_fmt,
            output,
            output_fmt,
        }
    }

    pub fn split_alignment(&self, part_path: &Path, partition_fmt: &PartitionFmt) {
        let partitions = PartitionParser::new(part_path, partition_fmt).parse();
        self.print_partition_info(part_path, &partitions.len());
        let spin = utils::set_spinner();
        spin.set_message("Parsing input sequence file...");
        let aln_matrix = self.parse_sequence();
        spin.set_message("Splitting alignment...");
        let file_counts = AtomicUsize::new(0);
        partitions.par_iter().for_each(|part| {
            let start_pos = part.start - 1;
            let end_pos = part.end;
            let matrix = self.generate_new_matrix(&aln_matrix, start_pos, end_pos);
            let mut header = Header::new();
            if let DataType::Aa = self.datatype {
                header.datatype = String::from("protein");
            }
            header.nchar = end_pos - start_pos;
            header.ntax = matrix.len();
            header.aligned = true;
            let filename = self.parse_filename(&part.gene);
            let output_path =
                filenames::create_output_fname(self.output, &filename, self.output_fmt);
            let mut out = SeqWriter::new(&output_path, &matrix, &header, None, &PartitionFmt::None);
            out.write_sequence(self.output_fmt)
                .expect("Failed writing the output file");
            file_counts.fetch_add(1, Ordering::Relaxed);
        });

        spin.finish_with_message("Finished splitting alignment!\n");

        self.print_output_info(file_counts.load(Ordering::Relaxed));
    }

    // Generate a filename for each locus based on the locus name
    // We get rid of any characters that are not alphanumeric, underscore, or dash
    fn parse_filename(&self, gene_name: &str) -> PathBuf {
        let mut filename = String::from(gene_name);
        filename.retain(|c| !r#"()/\,"';:?!"#.contains(c));
        if filename.contains('.') {
            filename = filename.replace('.', "_");
        }

        PathBuf::from(filename)
    }

    fn generate_new_matrix(
        &self,
        matrix: &SeqMatrix,
        start_pos: usize,
        end_pos: usize,
    ) -> SeqMatrix {
        let mut seq_matrix: SeqMatrix = IndexMap::new();
        matrix.iter().for_each(|(taxon, seq)| {
            // We substract the start position to match rust indexing.
            let part_seq = self.slice_sequence(seq, start_pos, end_pos);
            if !self.is_contain_seq(&part_seq) {
                seq_matrix.insert(taxon.clone(), part_seq);
            }
        });
        seq_matrix
    }

    fn is_contain_seq(&self, seq: &str) -> bool {
        let empty_seq = b"-?";
        seq.bytes().all(|char| empty_seq.contains(&char))
    }

    fn slice_sequence(&self, seq: &str, start: usize, end: usize) -> String {
        seq.get(start..end)
            .expect("Errors in slicing the alignment")
            .to_string()
    }

    fn parse_sequence(&self) -> SeqMatrix {
        let aln = Sequence::new(self.input, self.datatype);
        let (matrix, _) = aln.get_alignment(self.input_fmt);
        matrix
    }

    fn print_partition_info(&self, part_path: &Path, part_counts: &usize) {
        log::info!("{}", Yellow.paint("Partitions"));
        log::info!("{:18}: {}", "File path", part_path.display());
        log::info!(
            "{:18}: {}\n",
            "Partition counts",
            utils::fmt_num(part_counts)
        );
    }

    fn print_output_info(&self, file_counts: usize) {
        log::info!("{}", Yellow.paint("Output"));
        log::info!("{:18}: {}", "File counts", file_counts);
        log::info!("{:18}: {}", "Output dir", self.output.display());
        self.print_output_fmt(self.output_fmt);
    }
}

#[cfg(test)]
mod test {
    use super::*;

    macro_rules! input_split {
        ($var:ident, $input:expr) => {
            let $var = Splitter::new(
                &Path::new($input),
                &DataType::Dna,
                &InputFmt::Fasta,
                &Path::new("test_files/"),
                &OutputFmt::Fasta,
            );
        };
    }

    #[test]
    fn test_parse_filename() {
        input_split!(splitter, "test_files/test.fasta");
        assert_eq!(
            splitter.parse_filename(r#"'test!?'"#),
            PathBuf::from("test")
        );
    }

    #[test]
    fn test_sequence_slicing() {
        input_split!(split, "test_files/test_aln.fasta");
        let seq = "AAAAAAAAGGGGTTTTCCCC";

        let slice = split.slice_sequence(seq, 0, 10);
        let slice_2 = split.slice_sequence(seq, 10, 15);
        assert_eq!(slice, "AAAAAAAAGG");
        assert_eq!(slice_2, "GGTTT");
    }

    #[test]
    fn test_generate_new_matrix() {
        input_split!(split, "test_files/partition/concat_part.fas");
        let matrix = split.parse_sequence();
        let new_matrix = split.generate_new_matrix(&matrix, 0, 10);
        let new_matrix_2 = split.generate_new_matrix(&matrix, 10, 15);
        assert_eq!(new_matrix.len(), 4);
        assert_eq!(new_matrix.get("ABCD").unwrap(), "aaaaaggggg");
        assert_eq!(new_matrix_2.get("ABCD").unwrap(), "ttttt");
    }
}
