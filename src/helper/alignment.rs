use std::path::Path;

use indexmap::IndexMap;

use crate::helper::common::{self, Header, SeqFormat};
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

    pub fn get_aln_any(&mut self, file: &Path, input_fmt: &SeqFormat) {
        self.name
            .push_str(&file.file_stem().unwrap().to_string_lossy());
        match input_fmt {
            SeqFormat::Nexus => self.get_aln_from_nexus(file),
            SeqFormat::Phylip => self.get_aln_from_phylip(file, false),
            SeqFormat::PhylipInt => self.get_aln_from_phylip(file, true),
            SeqFormat::Fasta => self.get_aln_from_fasta(file),
            SeqFormat::Auto => {
                let input_fmt = common::infer_input_auto(file);
                self.get_aln_any(file, &input_fmt);
            }
            _ => (),
        }
    }

    fn get_aln_from_nexus(&mut self, file: &Path) {
        let mut nex = Nexus::new(file);
        nex.read().expect("CANNOT READ A NEXUS FILE");
        self.check_is_alignment(&file, nex.is_alignment);
        self.get_alignment(nex.matrix, nex.header)
    }

    fn get_aln_from_phylip(&mut self, file: &Path, interleave: bool) {
        let mut phy = Phylip::new(file, interleave);
        phy.read().expect("CANNOT READ A PHYLIP FILE");
        self.check_is_alignment(file, phy.is_alignment);
        self.get_alignment(phy.matrix, phy.header);
    }

    fn get_aln_from_fasta(&mut self, file: &Path) {
        let mut fas = Fasta::new(file);
        fas.read();
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
                "INVALID INPUT FILES. {} IS NOT AN ALIGNMENT",
                file.display()
            );
        }
    }
}
