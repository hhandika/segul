use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

pub fn parse_text_file(input: &Path) -> Vec<String> {
    let file = File::open(input).expect("Failed opening a text file");
    let buff = BufReader::new(file);
    let mut contents = Vec::new();
    buff.lines()
        .map_while(Result::ok)
        .for_each(|line| contents.push(line));
    contents
}
