//! Find input files and parse IDs from input files.

use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;

use glob::glob;
use indexmap::IndexSet;
use rayon::prelude::*;

use crate::helper::sequence;
use crate::helper::types::RawReadFmt;
use crate::helper::types::{DataType, InputFmt};
use crate::parser::fasta;
use crate::parser::nexus::Nexus;
use crate::parser::phylip::Phylip;

macro_rules! id_non_fasta {
    ($self:ident,  $type: ident, $datatype:ident) => {{
        let (sender, receiver) = channel();
        $self.files.par_iter().for_each_with(sender, |s, file| {
            s.send($type::new(file, $self.$datatype).parse_only_id())
                .expect("Failed parallel processing IDs");
        });
        receiver.iter().collect()
    }};
}

/// Find input files from a directory.
pub struct Files<'a> {
    /// Input directory.
    dir: &'a Path,
    /// Glob pattern.
    pattern: String,
}

impl<'a> Files<'a> {
    /// Create a new `Files` instance.
    pub fn new(dir: &'a Path) -> Self {
        Self {
            dir,
            pattern: String::new(),
        }
    }

    /// Find input files for sequence and alignment.
    /// Return a vector of input files.
    /// # Example
    /// ```
    /// use std::path::Path;
    /// use segul::helper::types::InputFmt;
    /// use segul::helper::finder::Files;
    ///
    /// let dir = Path::new("tests/files/concat");
    /// let input_fmt = InputFmt::Nexus;
    /// let files = Files::new(&dir).find(&input_fmt);
    /// assert_eq!(files.len(), 4);
    /// ```
    pub fn find(&mut self, input_fmt: &'a InputFmt) -> Vec<PathBuf> {
        self.pattern(input_fmt);
        let files = self.glob_files();
        self.check_glob_results(&files);

        files
    }

    /// Find input files for raw reads.
    /// Return a vector of input files.
    /// # Example
    /// ```
    /// use std::path::Path;
    /// use segul::helper::types::RawReadFmt;
    /// use segul::helper::finder::Files;
    ///
    /// let dir = Path::new("tests/files/raw");
    /// let input_fmt = RawReadFmt::Fastq;
    /// let files = Files::new(&dir).find_raw_read(&input_fmt);
    /// assert_eq!(files.len(), 4);
    pub fn find_raw_read(&mut self, input_fmt: &'a RawReadFmt) -> Vec<PathBuf> {
        self.raw_pattern(&input_fmt);
        let files = self.glob_files();
        self.check_glob_results(&files);

        files
    }

    fn glob_files(&self) -> Vec<PathBuf> {
        glob(&self.pattern)
            .expect("Failed globbing files")
            .filter_map(|ok| ok.ok())
            .collect::<Vec<PathBuf>>()
    }

    fn check_glob_results(&self, files: &[PathBuf]) {
        if files.is_empty() {
            panic!(
                "Failed finding input files using {}. \
                Check the input directory and the input format.",
                self.pattern
            );
        }
    }

    fn raw_pattern(&mut self, input_fmt: &'a RawReadFmt) {
        self.pattern = match input_fmt {
            RawReadFmt::Fastq => format!("{}/*.f*q", self.dir.display()),
            RawReadFmt::Gzip => format!("{}/*.f*q.gz*", self.dir.display()),
            RawReadFmt::Auto => panic!(
                "The input format is the default auto. \
            The program cannot use auto for dir input. \
            Try to specify input format."
            ),
        };
    }

    fn pattern(&mut self, input_fmt: &'a InputFmt) {
        self.pattern = match input_fmt {
            InputFmt::Fasta => format!("{}/*.fa*", self.dir.display()),
            InputFmt::Nexus => format!("{}/*.nex*", self.dir.display()),
            InputFmt::Phylip => format!("{}/*.phy*", self.dir.display()),
            InputFmt::Auto => panic!(
                "The input format is the default auto. \
            The program cannot use auto for dir input. \
            Try to specify input format."
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

    pub fn id_unique(&self) -> IndexSet<String> {
        let all_ids = self.parse_id();
        self.filter_unique(&all_ids)
    }

    fn parse_id(&self) -> Vec<IndexSet<String>> {
        match self.input_fmt {
            InputFmt::Nexus => id_non_fasta!(self, Nexus, datatype),
            InputFmt::Phylip => id_non_fasta!(self, Phylip, datatype),
            InputFmt::Fasta => self.id_from_fasta(),
            InputFmt::Auto => self.id_auto(),
        }
    }

    fn id_auto(&self) -> Vec<IndexSet<String>> {
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

    fn id_from_fasta(&self) -> Vec<IndexSet<String>> {
        let (sender, receiver) = channel();
        self.files.par_iter().for_each_with(sender, |s, file| {
            s.send(fasta::parse_only_id(file)).unwrap();
        });
        receiver.iter().collect()
    }

    fn filter_unique(&self, all_ids: &[IndexSet<String>]) -> IndexSet<String> {
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

    macro_rules! input {
        ($files: ident) => {
            let path = Path::new("tests/files/concat");

            let mut $files = Files::new(path);
        };
    }

    #[test]
    fn test_files() {
        input!(finder);
        let fmt = InputFmt::Nexus;
        let files = finder.find(&fmt);
        assert_eq!(4, files.len());
    }

    #[test]
    fn test_raw_pattern() {
        let path = Path::new("tests/files/raw");
        let fmt = RawReadFmt::Fastq;
        let mut files = Files::new(path);
        let found_files = files.find_raw_read(&fmt);
        assert_eq!(4, found_files.len());
        assert_eq!("tests/files/raw/*.f*q", files.pattern);
    }

    #[test]
    fn test_pattern() {
        input!(files);
        let fmt = InputFmt::Nexus;
        files.pattern(&fmt);
        assert_eq!("tests/files/concat/*.nex*", files.pattern);
    }

    #[test]
    #[should_panic]
    fn test_check_empty_files() {
        let path = Path::new("tests/files/empty/");
        let mut finder = Files::new(path);
        let files = finder.find(&InputFmt::Nexus);
        finder.check_glob_results(&files);
    }

    #[test]
    fn test_id() {
        input!(finder);
        let input_fmt = InputFmt::Nexus;
        let datatype = DataType::Dna;
        let files = finder.find(&input_fmt);
        let id = IDs::new(&files, &input_fmt, &datatype);
        let ids = id.id_unique();
        assert_eq!(3, ids.len());
    }
}
