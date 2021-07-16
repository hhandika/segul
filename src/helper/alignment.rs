use std::path::Path;

use indexmap::IndexMap;

use crate::helper::common::{self, Header, InputFmt};
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

    pub fn get_aln_any(&mut self, file: &Path, input_fmt: &InputFmt) {
        self.name
            .push_str(&file.file_stem().unwrap().to_string_lossy());
        match input_fmt {
            InputFmt::Nexus => self.get_aln_from_nexus(file),
            InputFmt::Phylip => self.get_aln_from_phylip(file, false),
            InputFmt::PhylipInt => self.get_aln_from_phylip(file, true),
            InputFmt::Fasta => self.get_aln_from_fasta(file),
            InputFmt::Auto => {
                let input_fmt = common::infer_input_auto(file);
                self.get_aln_any(file, &input_fmt);
            }
        }
    }

    fn get_aln_from_nexus(&mut self, file: &Path) {
        let mut nex = Nexus::new(file);
        nex.parse().expect("Failed reading a nexus file");
        self.check_is_alignment(&file, nex.is_alignment);
        self.get_alignment(nex.matrix, nex.header)
    }

    fn get_aln_from_phylip(&mut self, file: &Path, interleave: bool) {
        let mut phy = Phylip::new(file, interleave);
        phy.parse().expect("Failed reading a phylip file");
        self.check_is_alignment(file, phy.is_alignment);
        self.get_alignment(phy.matrix, phy.header);
    }

    fn get_aln_from_fasta(&mut self, file: &Path) {
        let mut fas = Fasta::new(file);
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
