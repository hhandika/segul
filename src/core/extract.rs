use std::path::PathBuf;

use regex::Regex;

pub struct Extract<'a> {
    re: &'a Option<String>,
}

impl<'a> Extract<'a> {
    pub fn new(re: &'a Option<String>) -> Self {
        Self { re }
    }

    #[allow(unused_variables)]
    pub fn extract_sequences(&self, files: &[PathBuf]) {
        match &self.re {
            Some(regex) => {
                self.match_id("Sequence 1", regex);
            }
            None => println!("No regex provided"),
        }
    }

    fn match_id(&self, id: &str, re: &str) -> bool {
        let re = Regex::new(re).expect("Failed capturing nexus commands");
        re.is_match(id)
    }
}
