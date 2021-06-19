use std::path::Path;

use indexmap::IndexMap;

use crate::common::InputFormat;
use crate::fasta::Fasta;
use crate::nexus::Nexus;
use crate::phylip::Phylip;

pub struct Alignment {
    pub alignment: IndexMap<String, String>,
    pub nchar: usize,
}

impl Alignment {
    pub fn new() -> Self {
        Self {
            alignment: IndexMap::new(),
            nchar: 0,
        }
    }

    pub fn get_aln_any(&mut self, file: &Path, input_format: &InputFormat) {
        match input_format {
            InputFormat::Nexus => self.get_aln_from_nexus(file),
            InputFormat::Phylip => self.get_aln_from_phylip(file),
            InputFormat::Fasta => self.get_aln_from_fasta(file),
        }
    }

    fn get_aln_from_nexus(&mut self, file: &Path) {
        let mut nex = Nexus::new(file);
        nex.read().expect("CANNOT READ A NEXUS FILE");
        self.check_is_alignment(&file, nex.is_alignment);
        self.get_alignment(nex.matrix, nex.nchar)
    }

    fn get_aln_from_phylip(&mut self, file: &Path) {
        let mut phy = Phylip::new(file);
        phy.read().expect("CANNOT READ A PHYLIP FILE");
        self.check_is_alignment(file, phy.is_alignment);
        self.get_alignment(phy.matrix, phy.nchar);
    }

    fn get_aln_from_fasta(&mut self, file: &Path) {
        let mut fas = Fasta::new(file);
        fas.read();
        self.check_is_alignment(file, fas.is_alignment);
        self.get_alignment(fas.matrix, fas.nchar);
    }

    fn get_alignment(&mut self, alignment: IndexMap<String, String>, nchar: usize) {
        self.alignment = alignment;
        self.nchar = nchar;
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
