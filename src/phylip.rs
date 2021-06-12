use std::collections::BTreeMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;

use nom::{character::complete, sequence, IResult};

pub fn read_phylip(path: &str) {
    let mut phylip = Phylip::new();
    phylip.read(path);

    // Test print to fasta
    phylip.matrix.iter().for_each(|(tax, seq)| {
        println!(">{}", tax);
        println!("{}", seq);
    })
}

struct Phylip {
    matrix: BTreeMap<String, String>,
    num_tax: usize,
    num_chars: usize,
}

impl Phylip {
    fn new() -> Self {
        Self {
            matrix: BTreeMap::new(),
            num_tax: 0,
            num_chars: 0,
        }
    }

    fn read(&mut self, path: &str) {
        let file = File::open(path).expect("CANNOT OPEN THE INPUT FILE.");
        let mut buff = BufReader::new(file);

        let mut header_line = String::new();
        buff.read_line(&mut header_line).unwrap();
        self.parse_header(&header_line.trim());

        buff.lines().filter_map(|ok| ok.ok()).for_each(|line| {
            let seq: Vec<&str> = line.split_whitespace().collect();
            match seq.len() {
                2 => self.matrix.insert(seq[0].to_string(), seq[1].to_string()),
                _ => {
                    panic!(
                        "UNSUPPORTED PHYLIP. \
                    THE PROGRAM ONLY WORK WITH NON-INTERLEAVED PHYLIP."
                    )
                }
            };
        })
    }

    fn parse_header<'a>(&'a mut self, header_line: &'a str) {
        let header: IResult<&str, (&str, &str)> =
            sequence::separated_pair(complete::digit0, complete::space0, complete::digit0)(
                header_line,
            );

        match header {
            Ok((_, (tax, chars))) => self.parse_num(tax, chars),
            Err(_) => eprintln!("UNKNOWN HEADER! FAILED TO PARSE"),
        };
    }

    fn parse_num(&mut self, tax: &str, seq: &str) {
        self.num_tax = tax
            .parse::<usize>()
            .expect("HEADER TAXA NUMBER IS NOT A NUMBER");
        self.num_chars = seq
            .parse::<usize>()
            .expect("HEADER CHARS LENGTH IS NOT A NUMBER");
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn read_phylip_simple_test() {
        let path = "test_files/simple.phy";
        let mut phylip = Phylip::new();
        phylip.read(path);

        assert_eq!(2, phylip.num_tax);
        assert_eq!(4, phylip.num_chars);
        assert_eq!(2, phylip.matrix.len());
    }
}
