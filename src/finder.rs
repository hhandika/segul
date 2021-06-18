//! Module for files finding and IDs indexing.

use std::path::PathBuf;

use glob::glob;
use indexmap::IndexMap;
use indexmap::IndexSet;

use crate::common::InputFormat;
use crate::fasta::Fasta;
use crate::nexus::Nexus;
use crate::phylip::Phylip;

pub struct Files<'a> {
    dir: &'a str,
    input_format: &'a InputFormat,
    pattern: String,
}

impl<'a> Files<'a> {
    pub fn new(dir: &'a str, input_format: &'a InputFormat) -> Self {
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
            InputFormat::Nexus => format!("{}/*.nex*", self.dir),
            InputFormat::Phylip => format!("{}/*.phy*", self.dir),
            InputFormat::Fasta => format!("{}/*.fa*", self.dir),
        };
    }
}

pub struct IDs<'a> {
    files: &'a [PathBuf],
    input_format: &'a InputFormat,
}

impl<'a> IDs<'a> {
    pub fn new(files: &'a [PathBuf], input_format: &'a InputFormat) -> Self {
        Self {
            files,
            input_format,
        }
    }

    pub fn get_id_all(&self) -> IndexSet<String> {
        let mut id = IndexSet::new();
        match self.input_format {
            InputFormat::Nexus => self.get_id_from_nexus(&mut id),
            InputFormat::Phylip => self.get_id_from_phylip(&mut id),
            InputFormat::Fasta => self.get_id_from_fasta(&mut id),
        };
        id
    }

    fn get_id_from_phylip(&self, id: &mut IndexSet<String>) {
        self.files.iter().for_each(|file| {
            let mut phy = Phylip::new(file);
            phy.read().expect("CANNOT READ A PHYLIP FILE");
            self.get_id(&phy.matrix, id);
        });
    }

    fn get_id_from_nexus(&self, id: &mut IndexSet<String>) {
        self.files.iter().for_each(|file| {
            let mut nex = Nexus::new(file);
            nex.read().expect("CANNOT READ A NEXUS FILE");
            self.get_id(&nex.matrix, id);
        });
    }

    fn get_id_from_fasta(&self, id: &mut IndexSet<String>) {
        self.files.iter().for_each(|file| {
            let mut fas = Fasta::new(file);
            fas.read();
            self.get_id(&fas.matrix, id);
        });
    }

    fn get_id(&self, matrix: &IndexMap<String, String>, id: &mut IndexSet<String>) {
        matrix.keys().for_each(|key| {
            if !id.contains(key) {
                id.insert(key.to_string());
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
        let mut finder = Files::new(path, &InputFormat::Nexus);
        let files = finder.get_files();
        assert_eq!(4, files.len());
    }

    #[test]
    #[should_panic]
    fn check_empty_files_test() {
        let path = "test_files/empty/";
        let mut finder = Files::new(path, &InputFormat::Nexus);
        let files = finder.get_files();
        finder.check_glob_results(&files);
    }
}
