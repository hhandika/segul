use std::ffi::OsStr;
use std::path::Path;

use indexmap::IndexMap;

use crate::helper::common::{self, DataType, Header, InputFmt};
use crate::parser::fasta::Fasta;
use crate::parser::nexus::Nexus;
use crate::parser::phylip::Phylip;

pub struct Sequence {
    pub alignment: IndexMap<String, String>,
    pub header: Header,
    pub name: String,
}

impl Sequence {
    pub fn new() -> Self {
        Self {
            alignment: IndexMap::new(),
            header: Header::new(),
            name: String::new(),
        }
    }

    pub fn get_aln_any(&mut self, file: &Path, input_fmt: &InputFmt, datatype: &DataType) {
        self.name.push_str(
            file.file_stem()
                .and_then(OsStr::to_str)
                .expect("Failed getting alignment name from the file"),
        );

        match input_fmt {
            InputFmt::Nexus => self.from_nexus(file, datatype),
            InputFmt::Phylip => self.from_phylip(file, datatype),
            InputFmt::Fasta => self.from_fasta(file, datatype),
            InputFmt::Auto => {
                let input_fmt = common::infer_input_auto(file);
                self.get_aln_any(file, &input_fmt, datatype);
            }
        }

        self.check_is_alignment(file);
    }

    fn from_nexus(&mut self, file: &Path, datatype: &DataType) {
        let mut nex = Nexus::new(file, datatype);
        nex.parse().expect("Failed reading a nexus file");
        self.get(nex.matrix, nex.header)
    }

    fn from_phylip(&mut self, file: &Path, datatype: &DataType) {
        let mut phy = Phylip::new(file, datatype);
        phy.parse().expect("Failed reading a phylip file");
        self.get(phy.matrix, phy.header);
    }

    fn from_fasta(&mut self, file: &Path, datatype: &DataType) {
        let mut fas = Fasta::new(file, datatype);
        fas.parse();
        self.get(fas.matrix, fas.header);
    }

    fn get(&mut self, alignment: IndexMap<String, String>, header: Header) {
        self.alignment = alignment;
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
        aln.get_aln_any(&mut file, &input_fmt, &datatype);
        assert_eq!(String::from("simple"), aln.name);
        assert_eq!(1, aln.header.ntax);
        assert_eq!(6, aln.header.nchar);
        assert_eq!(1, aln.alignment.len());
    }
}
