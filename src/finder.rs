//! Module for files finding and IDs indexing.

use std::path::PathBuf;

use glob::glob;
// use indexmap::IndexMap;
use indexmap::IndexSet;
// use rayon::prelude::*;

use crate::common::SeqFormat;
use crate::fasta;
use crate::nexus::Nexus;
use crate::phylip::Phylip;

pub struct Files<'a> {
    dir: &'a str,
    input_format: &'a SeqFormat,
    pattern: String,
}

impl<'a> Files<'a> {
    pub fn new(dir: &'a str, input_format: &'a SeqFormat) -> Self {
        Self {
            dir,
            input_format,
            pattern: String::new(),
        }
    }

    pub fn get_files(&mut self) -> Vec<PathBuf> {
        self.get_pattern();
        let files = glob(&self.pattern)
            .expect("COULD NOT FIND FILES")
            .filter_map(|ok| ok.ok())
            .collect::<Vec<PathBuf>>();
        self.check_glob_results(&files);

        files
    }

    fn check_glob_results(&self, files: &[PathBuf]) {
        if files.is_empty() {
            panic!("NO VALID ALIGNMENT FILES FOUND");
        }
    }

    fn get_pattern(&mut self) {
        self.pattern = match self.input_format {
            SeqFormat::Nexus => format!("{}/*.nex*", self.dir),
            SeqFormat::Phylip => format!("{}/*.phy*", self.dir),
            SeqFormat::Fasta => format!("{}/*.fa*", self.dir),
        };
    }
}

pub struct IDs<'a> {
    files: &'a [PathBuf],
    input_format: &'a SeqFormat,
}

impl<'a> IDs<'a> {
    pub fn new(files: &'a [PathBuf], input_format: &'a SeqFormat) -> Self {
        Self {
            files,
            input_format,
        }
    }

    pub fn get_id_all(&self, interleave: bool) -> IndexSet<String> {
        let mut id = IndexSet::new();
        match self.input_format {
            SeqFormat::Nexus => self.get_id_from_nexus(&mut id),
            SeqFormat::Phylip => self.get_id_from_phylip(&mut id, interleave),
            SeqFormat::Fasta => self.get_id_from_fasta(&mut id),
        };
        id
    }

    fn get_id_from_phylip(&self, id: &mut IndexSet<String>, interleave: bool) {
        self.files.iter().for_each(|file| {
            let mut phy = Phylip::new(file, interleave);
            let ids = phy.read_only_id();
            self.get_id(&ids, id);
        });
    }

    fn get_id_from_nexus(&self, id: &mut IndexSet<String>) {
        self.files.iter().for_each(|file| {
            let mut nex = Nexus::new(file);
            let ids = nex.read_only_id();
            self.get_id(&ids, id);
        });
    }

    fn get_id_from_fasta(&self, id: &mut IndexSet<String>) {
        self.files.iter().for_each(|file| {
            let ids = fasta::read_only_id(file);
            self.get_id(&ids, id)
        });
    }

    fn get_id(&self, ids: &IndexSet<String>, id: &mut IndexSet<String>) {
        ids.iter().for_each(|val| {
            if !id.contains(val) {
                id.insert(val.to_string());
            }
        });
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn get_files_test() {
        let path = "test_files/concat/";
        let mut finder = Files::new(path, &SeqFormat::Nexus);
        let files = finder.get_files();
        assert_eq!(4, files.len());
    }

    #[test]
    #[should_panic]
    fn check_empty_files_test() {
        let path = "test_files/empty/";
        let mut finder = Files::new(path, &SeqFormat::Nexus);
        let files = finder.get_files();
        finder.check_glob_results(&files);
    }
}
