use std::ffi::OsStr;
use std::path::Path;

use indexmap::IndexMap;

use crate::helper::common::{self, DataType, Header, InputFmt};
use crate::parser::fasta::Fasta;
use crate::parser::nexus::Nexus;
use crate::parser::phylip::Phylip;

pub struct Alignment {
    pub alignment: IndexMap<String, String>,
    pub header: Header,
    pub name: String,
}

impl Alignment {
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
            InputFmt::Nexus => self.get_aln_from_nexus(file, datatype),
            InputFmt::Phylip => self.get_aln_from_phylip(file, datatype),
            InputFmt::Fasta => self.get_aln_from_fasta(file, datatype),
            InputFmt::Auto => {
                let input_fmt = common::infer_input_auto(file);
                self.get_aln_any(file, &input_fmt, datatype);
            }
        }
    }

    fn get_aln_from_nexus(&mut self, file: &Path, datatype: &DataType) {
        let mut nex = Nexus::new(file, datatype);
        nex.parse().expect("Failed reading a nexus file");
        self.check_is_alignment(&file, nex.is_alignment);
        self.get_alignment(nex.matrix, nex.header)
    }

    fn get_aln_from_phylip(&mut self, file: &Path, datatype: &DataType) {
        let mut phy = Phylip::new(file, datatype);
        phy.parse().expect("Failed reading a phylip file");
        self.check_is_alignment(file, phy.is_alignment);
        self.get_alignment(phy.matrix, phy.header);
    }

    fn get_aln_from_fasta(&mut self, file: &Path, datatype: &DataType) {
        let mut fas = Fasta::new(file, datatype);
        fas.parse();
        self.check_is_alignment(file, fas.is_alignment);
        self.get_alignment(fas.matrix, fas.header);
    }

    fn get_alignment(&mut self, alignment: IndexMap<String, String>, header: Header) {
        self.alignment = alignment;
        self.header = header;
    }

    fn check_is_alignment(&self, file: &Path, aligned: bool) {
        if !aligned {
            panic!(
                "Ups. Something is wrong. {} is not an alignment",
                file.display()
            );
        }
    }
}
