use std::path::PathBuf;

use regex::Regex;

pub enum Params {
    Regex(String),
    File(PathBuf),
    None,
}

pub struct Extract<'a> {
    params: &'a Params,
}

impl<'a> Extract<'a> {
    pub fn new(params: &'a Params) -> Self {
        Self { params }
    }

    #[allow(unused_variables)]
    pub fn extract_sequences(&self, files: &[PathBuf]) {
        match self.params {
            Params::Regex(re) => {
                let re = self.match_id("Sequence 1", re);
                println!("Match: {}\n", re);
            }
            Params::File(path) => println!("Path: {}\n", path.display()),
            Params::None => unreachable!("Please, specify a matching parameter!"),
        }
    }

    fn match_id(&self, id: &str, re: &str) -> bool {
        let re = Regex::new(re).expect("Failed capturing nexus commands");
        re.is_match(id)
    }
}
