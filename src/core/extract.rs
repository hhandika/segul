// use lazy_static::lazy_static;
use regex::Regex;

pub struct Extract {
    re: Option<String>,
}

impl Extract {
    pub fn new(re: Option<String>) -> Self {
        Self { re }
    }

    pub fn print_input_info(&self) {
        let text = "Sequence_1 ATGTGTG";
        println!("{:18}: {}", "Regex", self.re.as_ref().unwrap());

        match &self.re {
            Some(re) => {
                let res = self.match_id(text, re);
                println!("{:18}: {}", "Result", res);
            }
            None => println!("None"),
        }
    }

    fn match_id(&self, id: &str, re: &str) -> String {
        let re = Regex::new(re).unwrap();
        match re.captures(id) {
            Some(word) => word[0].to_lowercase(),
            None => String::from("No match"),
        }
    }
}
