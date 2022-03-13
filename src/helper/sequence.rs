use std::ffi::OsStr;
use std::path::Path;

use crate::helper::types::{DataType, Header, InputFmt, SeqMatrix};
use crate::parser::fasta::Fasta;
use crate::parser::nexus::Nexus;
use crate::parser::phylip::Phylip;

macro_rules! parse_sequence {
    ($self:ident, $format:ident) => {{
        let mut seq = $format::new($self.file, $self.datatype);
        seq.parse();
        (seq.matrix, seq.header)
    }};
}

pub fn infer_input_auto(input: &Path) -> InputFmt {
    let ext: &str = input
        .extension()
        .and_then(OsStr::to_str)
        .expect("Failed parsing extension");
    match ext {
        "fas" | "fa" | "fasta" => InputFmt::Fasta,
        "nex" | "nexus" => InputFmt::Nexus,
        "phy" | "phylip" => InputFmt::Phylip,
        _ => panic!(
            "Ups... The program cannot recognize the file extension. \
        Maybe try to specify the input format using the -f or --input-format option."
        ),
    }
}

pub struct Sequence<'a> {
    file: &'a Path,
    datatype: &'a DataType,
}

impl<'a> Sequence<'a> {
    pub fn new(file: &'a Path, datatype: &'a DataType) -> Self {
        Self { file, datatype }
    }

    pub fn get_alignment(&self, input_fmt: &'a InputFmt) -> (SeqMatrix, Header) {
        let (matrix, header) = self.get(input_fmt);
        assert!(
            header.aligned,
            "Ups. Something is wrong. {} is not an alignment",
            self.file.display()
        );
        (matrix, header)
    }

    pub fn get(&self, input_fmt: &'a InputFmt) -> (SeqMatrix, Header) {
        match input_fmt {
            InputFmt::Fasta => parse_sequence!(self, Fasta),
            InputFmt::Nexus => parse_sequence!(self, Nexus),
            InputFmt::Phylip => parse_sequence!(self, Phylip),
            InputFmt::Auto => {
                let input_fmt = infer_input_auto(self.file);
                self.get(&input_fmt)
            }
        }
    }
}

pub struct SeqCheck {
    pub shortest: usize,
    pub longest: usize,
    pub is_alignment: bool,
}

impl SeqCheck {
    pub fn new() -> Self {
        Self {
            shortest: 0,
            longest: 0,
            is_alignment: false,
        }
    }

    pub fn check(&mut self, matrix: &SeqMatrix) {
        assert!(
            !matrix.is_empty(),
            "The data matrix is empty. \
        Make user the input format is correct."
        );
        self.get_shortest_seq_len(matrix);
        self.get_longest_seq_len(matrix);
        self.check_is_alignment();
    }

    fn check_is_alignment(&mut self) {
        self.is_alignment = self.shortest == self.longest;
    }

    fn get_shortest_seq_len(&mut self, matrix: &SeqMatrix) {
        self.shortest = matrix
            .values()
            .map(|s| s.len())
            .min_by(|a, b| a.cmp(b))
            .expect("Failed getting the shortest sequence length");
    }

    fn get_longest_seq_len(&mut self, matrix: &SeqMatrix) {
        self.longest = matrix
            .values()
            .map(|s| s.len())
            .max_by(|a, b| a.cmp(b))
            .expect("Failed getting the longest sequence length");
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_alignment_simple() {
        let file = Path::new("test_files/simple.nex");
        let datatype = DataType::Dna;
        let input_fmt = InputFmt::Nexus;
        let aln = Sequence::new(&file, &datatype);
        let (matrix, header) = aln.get_alignment(&input_fmt);
        assert_eq!(1, header.ntax);
        assert_eq!(6, header.nchar);
        assert_eq!(1, matrix.len());
    }

    #[test]
    fn test_parsing_input_fmt() {
        let file = Path::new("test_files/simple.nex");
        let input_fmt = infer_input_auto(&file);
        assert_eq!(InputFmt::Nexus, input_fmt);
    }
}
