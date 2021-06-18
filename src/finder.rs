use std::path::PathBuf;

use glob::glob;

use crate::common::InputFormat;

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
