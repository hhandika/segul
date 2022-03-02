#[allow(dead_code, unused_imports)]
use std::ffi::OsStr;
use std::path::Path;

#[allow(dead_code, unused_imports)]
use crate::helper::finder::IDs;
use crate::helper::sequence::Sequence;
use crate::helper::types::{DataType, InputFmt, OutputFmt, PartitionFmt, SeqMatrix};
// use crate::helper::utils;
use crate::parser::partition::PartitionParser;
// use crate::writer::sequences::SeqWriter;

#[allow(dead_code)]
pub struct Splitter<'a> {
    input: &'a Path,
    datatype: &'a DataType,
    input_fmt: &'a InputFmt,
    output: &'a Path,
    output_fmt: &'a OutputFmt,
}

#[allow(dead_code, unused_variables)]
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

    fn split_alignment(&self, partition: &Path, partition_fmt: &PartitionFmt) {
        let partitions = PartitionParser::new(partition, partition_fmt).parse();
        let aln_matrix = self.parse_sequence();
    }

    fn slice_sequence(&self, seq: &str, start: &usize, end: &usize) -> String {
        seq.get(*start..*end)
            .expect("Errors in slicing the alignment")
            .to_string()
    }

    fn parse_sequence(&self) -> SeqMatrix {
        let aln = Sequence::new(self.input, self.datatype);
        let (matrix, _) = aln.get_alignment(self.input_fmt);
        matrix
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_sequence_slicing() {
        let split = Splitter::new(
            &Path::new("tests/data/alignment.fasta"),
            &DataType::Dna,
            &InputFmt::Fasta,
            &Path::new("tests/data/alignment.fasta"),
            &OutputFmt::Fasta,
        );
        let seq = "AAAAAAAAGGGGTTTTCCCC";

        let slice = split.slice_sequence(seq, &0, &10);
        let slice_2 = split.slice_sequence(seq, &10, &15);
        assert_eq!(slice, "AAAAAAAAGG");
        assert_eq!(slice_2, "GGTTT");
    }
}
