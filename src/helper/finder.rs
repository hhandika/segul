//! Module for files finding and IDs indexing.

use std::path::PathBuf;
use std::sync::mpsc::channel;

use glob::glob;
use indexmap::IndexSet;
use rayon::prelude::*;

use crate::helper::common::SeqFormat;
use crate::parser::fasta;
use crate::parser::nexus::Nexus;
use crate::parser::phylip::Phylip;

pub struct Files<'a> {
    dir: &'a str,
    input_fmt: &'a SeqFormat,
    pattern: String,
}

impl<'a> Files<'a> {
    pub fn new(dir: &'a str, input_fmt: &'a SeqFormat) -> Self {
        Self {
            dir,
            input_fmt,
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
            panic!("NO ALIGNMENT FILES FOUND THAT MATCH {}", self.pattern);
        }
    }

    fn get_pattern(&mut self) {
        self.pattern = match self.input_fmt {
            SeqFormat::Fasta => format!("{}/*.fa*", self.dir),
            SeqFormat::Nexus => format!("{}/*.nex*", self.dir),
            SeqFormat::Phylip | SeqFormat::PhylipInt => format!("{}/*.phy*", self.dir),
            SeqFormat::Auto => panic!(
                "YOUR INPUT FORMAT IS THE DEFAULT AUTO. \
            THE PROGRAM CANNOT USE AUTO FOR DIR INPUT. PLEASE, \
            SPECIFY INPUT FORMAT USING THE OPTION -f or --format OR USE USE WILDCARD."
            ),
            _ => panic!("UNSUPPORTED FORMAT! USE FASTA, NEXUS, OR PHYLIP ONLY"),
        };
    }
}

pub struct IDs<'a> {
    files: &'a [PathBuf],
    input_fmt: &'a SeqFormat,
}

impl<'a> IDs<'a> {
    pub fn new(files: &'a [PathBuf], input_fmt: &'a SeqFormat) -> Self {
        Self { files, input_fmt }
    }

    pub fn get_id_all(&self) -> IndexSet<String> {
        let all_ids = match self.input_fmt {
            SeqFormat::Nexus => self.get_id_from_nexus(),
            SeqFormat::Phylip => self.get_id_from_phylip(false),
            SeqFormat::PhylipInt => self.get_id_from_phylip(true),
            SeqFormat::Fasta => self.get_id_from_fasta(),
            _ => panic!("USE FASTA, NEXUS, OR PHYLIP ONLY"),
        };
        self.get_id(&all_ids)
    }

    fn get_id_from_phylip(&self, interleave: bool) -> Vec<IndexSet<String>> {
        let (sender, receiver) = channel();
        self.files.par_iter().for_each_with(sender, |s, file| {
            s.send(Phylip::new(file, interleave).read_only_id())
                .unwrap();
        });
        receiver.iter().collect()
    }

    fn get_id_from_nexus(&self) -> Vec<IndexSet<String>> {
        let (sender, receiver) = channel();
        self.files.par_iter().for_each_with(sender, |s, file| {
            s.send(Nexus::new(file).read_only_id()).unwrap();
        });
        receiver.iter().collect()
    }

    fn get_id_from_fasta(&self) -> Vec<IndexSet<String>> {
        let (sender, receiver) = channel();
        self.files.par_iter().for_each_with(sender, |s, file| {
            s.send(fasta::read_only_id(file)).unwrap();
        });
        receiver.iter().collect()
    }

    fn get_id(&self, all_ids: &[IndexSet<String>]) -> IndexSet<String> {
        let mut id = IndexSet::new();
        all_ids.iter().for_each(|ids| {
            ids.iter().for_each(|val| {
                if !id.contains(val) {
                    id.insert(val.to_string());
                }
            });
        });

        id
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
