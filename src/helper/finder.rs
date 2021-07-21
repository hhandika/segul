//! Module for files finding and IDs indexing.

use std::path::PathBuf;
use std::sync::mpsc::channel;

use glob::glob;
use indexmap::IndexSet;
use rayon::prelude::*;

use crate::helper::sequence;
use crate::helper::types::{DataType, InputFmt};
use crate::parser::fasta;
use crate::parser::nexus::Nexus;
use crate::parser::phylip::Phylip;

pub struct Files<'a> {
    dir: &'a str,
    input_fmt: &'a InputFmt,
    pattern: String,
}

impl<'a> Files<'a> {
    pub fn new(dir: &'a str, input_fmt: &'a InputFmt) -> Self {
        Self {
            dir,
            input_fmt,
            pattern: String::new(),
        }
    }

    pub fn get_files(&mut self) -> Vec<PathBuf> {
        self.get_pattern();
        let files = glob(&self.pattern)
            .expect("Failed globbing files")
            .filter_map(|ok| ok.ok())
            .collect::<Vec<PathBuf>>();
        self.check_glob_results(&files);

        files
    }

    fn check_glob_results(&self, files: &[PathBuf]) {
        if files.is_empty() {
            panic!(
                "Failed finding files that match {}. \
            Maybe try using wildcard option -c or --wcard",
                self.pattern
            );
        }
    }

    fn get_pattern(&mut self) {
        self.pattern = match self.input_fmt {
            InputFmt::Fasta => format!("{}/*.fa*", self.dir),
            InputFmt::Nexus => format!("{}/*.nex*", self.dir),
            InputFmt::Phylip => format!("{}/*.phy*", self.dir),
            InputFmt::Auto => panic!(
                "The input format is the default auto. \
            The program cannot use auto for dir input. \
            Try to specify input format using the option -f or --format \
            or use the wildcard option -c or --wcard."
            ),
        };
    }
}

pub struct IDs<'a> {
    files: &'a [PathBuf],
    input_fmt: &'a InputFmt,
    datatype: &'a DataType,
}

impl<'a> IDs<'a> {
    pub fn new(files: &'a [PathBuf], input_fmt: &'a InputFmt, datatype: &'a DataType) -> Self {
        Self {
            files,
            input_fmt,
            datatype,
        }
    }

    pub fn get_id_all(&self) -> IndexSet<String> {
        let all_ids = match self.input_fmt {
            InputFmt::Nexus => self.get_id_from_nexus(),
            InputFmt::Phylip => self.get_id_from_phylip(),
            InputFmt::Fasta => self.get_id_from_fasta(),
            InputFmt::Auto => self.get_id_auto(),
        };
        self.get_id(&all_ids)
    }

    fn get_id_auto(&self) -> Vec<IndexSet<String>> {
        let (sender, receiver) = channel();
        self.files.par_iter().for_each_with(sender, |s, file| {
            let input_fmt = sequence::infer_input_auto(file);
            match input_fmt {
                InputFmt::Fasta => s.send(fasta::parse_only_id(file)).unwrap(),
                InputFmt::Nexus => s
                    .send(Nexus::new(file, self.datatype).parse_only_id())
                    .unwrap(),
                InputFmt::Phylip => s
                    .send(Phylip::new(file, self.datatype).parse_only_id())
                    .unwrap(),
                _ => unreachable!(),
            }
        });
        receiver.iter().collect()
    }

    fn get_id_from_phylip(&self) -> Vec<IndexSet<String>> {
        let (sender, receiver) = channel();
        self.files.par_iter().for_each_with(sender, |s, file| {
            s.send(Phylip::new(file, self.datatype).parse_only_id())
                .unwrap();
        });
        receiver.iter().collect()
    }

    fn get_id_from_nexus(&self) -> Vec<IndexSet<String>> {
        let (sender, receiver) = channel();
        self.files.par_iter().for_each_with(sender, |s, file| {
            s.send(Nexus::new(file, self.datatype).parse_only_id())
                .unwrap();
        });
        receiver.iter().collect()
    }

    fn get_id_from_fasta(&self) -> Vec<IndexSet<String>> {
        let (sender, receiver) = channel();
        self.files.par_iter().for_each_with(sender, |s, file| {
            s.send(fasta::parse_only_id(file)).unwrap();
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
        let mut finder = Files::new(path, &InputFmt::Nexus);
        let files = finder.get_files();
        assert_eq!(4, files.len());
    }

    #[test]
    #[should_panic]
    fn check_empty_files_test() {
        let path = "test_files/empty/";
        let mut finder = Files::new(path, &InputFmt::Nexus);
        let files = finder.get_files();
        finder.check_glob_results(&files);
    }
}
