use std::path::Path;

use indexmap::IndexMap;

use crate::helper::common::{self, DataType, Header, InputFmt};
use crate::parse_sequence;
use crate::parser::fasta::Fasta;
use crate::parser::nexus::Nexus;
use crate::parser::phylip::Phylip;

pub struct Sequence<'a> {
    file: &'a Path,
    datatype: &'a DataType,
}

impl<'a> Sequence<'a> {
    pub fn new(file: &'a Path, datatype: &'a DataType) -> Self {
        Self {
            file,
            // input_fmt,
            datatype,
        }
    }

    #[inline]
    pub fn get_alignment(&self, input_fmt: &'a InputFmt) -> (IndexMap<String, String>, Header) {
        let (matrix, header) = self.get(input_fmt);
        self.check_is_alignment(self.file, &header);
        (matrix, header)
    }

    pub fn get(&self, input_fmt: &'a InputFmt) -> (IndexMap<String, String>, Header) {
        match input_fmt {
            InputFmt::Nexus => parse_sequence!(self, file, datatype, Nexus),
            InputFmt::Phylip => parse_sequence!(self, file, datatype, Phylip),
            InputFmt::Fasta => parse_sequence!(self, file, datatype, Fasta),
            InputFmt::Auto => {
                let input_fmt = common::infer_input_auto(self.file);
                self.get(&input_fmt)
            }
        }
    }

    fn check_is_alignment(&self, file: &Path, header: &Header) {
        if !header.aligned {
            panic!(
                "Ups. Something is wrong. {} is not an alignment",
                file.display()
            );
        }
    }
}

#[macro_export]
macro_rules! parse_sequence {
    ($self:ident, $file:ident, $datatype:ident, $format:ident) => {{
        let mut seq = $format::new($self.$file, $self.$datatype);
        seq.parse();
        (seq.matrix, seq.header)
    }};
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn alignment_simple_test() {
        let file = Path::new("test_files/simple.nex");
        let datatype = DataType::Dna;
        let input_fmt = InputFmt::Nexus;
        let aln = Sequence::new(&file, &datatype);
        let (matrix, header) = aln.get_alignment(&input_fmt);
        assert_eq!(1, header.ntax);
        assert_eq!(6, header.nchar);
        assert_eq!(1, matrix.len());
    }
}
