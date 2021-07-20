use std::ffi::OsStr;
use std::path::Path;

use indexmap::IndexMap;

use crate::helper::common::{self, DataType, Header, InputFmt};
use crate::parser::fasta::Fasta;
use crate::parser::nexus::Nexus;
use crate::parser::phylip::Phylip;

pub struct Sequence {
    pub matrix: IndexMap<String, String>,
    pub header: Header,
    pub name: String,
}

impl Sequence {
    pub fn new() -> Self {
        Self {
            matrix: IndexMap::new(),
            header: Header::new(),
            name: String::new(),
        }
    }

    #[inline]
    pub fn get_alignment(&mut self, file: &Path, input_fmt: &InputFmt, datatype: &DataType) {
        self.get(file, input_fmt, datatype);
        self.check_is_alignment(file);
    }

    pub fn get(&mut self, file: &Path, input_fmt: &InputFmt, datatype: &DataType) {
        self.name.push_str(
            &file
                .file_stem()
                .and_then(OsStr::to_str)
                .expect("Failed getting alignment name from the file"),
        );

        match input_fmt {
            InputFmt::Nexus => self.from_nexus(file, datatype),
            InputFmt::Phylip => self.from_phylip(file, datatype),
            InputFmt::Fasta => self.from_fasta(file, datatype),
            InputFmt::Auto => {
                let input_fmt = common::infer_input_auto(file);
                self.get(file, &input_fmt, datatype);
            }
        }
    }

    fn from_nexus(&mut self, file: &Path, datatype: &DataType) {
        let mut nex = Nexus::new(file, datatype);
        nex.parse().expect("Failed reading a nexus file");
        self.get_sequence(nex.matrix, nex.header)
    }

    fn from_phylip(&mut self, file: &Path, datatype: &DataType) {
        let mut phy = Phylip::new(file, datatype);
        phy.parse().expect("Failed reading a phylip file");
        self.get_sequence(phy.matrix, phy.header);
    }

    fn from_fasta(&mut self, file: &Path, datatype: &DataType) {
        let mut fas = Fasta::new(file, datatype);
        fas.parse();
        self.get_sequence(fas.matrix, fas.header);
    }

    #[inline]
    fn get_sequence(&mut self, matrix: IndexMap<String, String>, header: Header) {
        self.matrix = matrix;
        self.header = header;
    }

    fn check_is_alignment(&self, file: &Path) {
        if !self.header.aligned {
            panic!(
                "Ups. Something is wrong. {} is not an alignment",
                file.display()
            );
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn alignment_simple_test() {
        let mut aln = Sequence::new();
        let mut file = Path::new("test_files/simple.nex");
        let datatype = DataType::Dna;
        let input_fmt = InputFmt::Nexus;
        aln.get_alignment(&mut file, &input_fmt, &datatype);
        assert_eq!(String::from("simple"), aln.name);
        assert_eq!(1, aln.header.ntax);
        assert_eq!(6, aln.header.nchar);
        assert_eq!(1, aln.matrix.len());
    }
}
