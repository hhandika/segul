use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Result};
use std::path::Path;

use indexmap::IndexMap;
use nom::{character::complete, sequence, IResult};

use crate::common::{self, Header, SeqCheck};

pub struct Phylip<'a> {
    input: &'a Path,
    pub matrix: IndexMap<String, String>,
    pub ntax: usize,
    pub nchar: usize,
    pub is_alignment: bool,
}

impl SeqCheck for Phylip<'_> {}

impl<'a> Phylip<'a> {
    pub fn new(input: &'a Path) -> Self {
        Self {
            input,
            matrix: IndexMap::new(),
            ntax: 0,
            nchar: 0,
            is_alignment: false,
        }
    }

    pub fn read(&mut self) -> Result<()> {
        let file = File::open(self.input).expect("CANNOT OPEN THE INPUT FILE.");
        let mut buff = BufReader::new(file);

        let mut header_line = String::new();
        buff.read_line(&mut header_line)?;
        self.parse_header(&header_line.trim());

        buff.lines().filter_map(|ok| ok.ok()).for_each(|line| {
            let (id, dna) = self.parse_sequence(line.trim());
            self.insert_matrix(id, dna);
        });

        let (shortest, longest) = self.get_sequence_len(&self.matrix);
        self.is_alignment = self.check_is_alignment(&shortest, &longest);
        self.check_ntax_matches();
        self.check_nchar_matches(longest);
        Ok(())
    }

    pub fn get_header(&self) -> Header {
        let mut header = Header::new();
        header.ntax = Some(self.ntax);
        header.nchar = Some(self.nchar);
        header
    }

    fn parse_sequence(&self, line: &str) -> (String, String) {
        let seq: Vec<&str> = line.split_whitespace().collect();
        self.check_seq_len(seq.len());
        let id = seq[0].to_string();
        let dna = seq[1].to_string().to_lowercase();
        common::check_valid_dna(&self.input, &id, &dna);
        (id, dna)
    }

    fn insert_matrix(&mut self, id: String, dna: String) {
        #[allow(clippy::all)]
        if self.matrix.contains_key(&id) {
            panic!(
                "DUPLICATE SAMPLES FOR FILE {}. FIRST DUPLICATE FOUND: {}",
                self.input.display(),
                id
            );
        } else {
            self.matrix.insert(id, dna);
        }
    }

    fn check_seq_len(&self, len: usize) {
        if len != 2 {
            panic!(
                "THE FILE {} IS UNSUPPORTED PHYLIP FORMAT. \
            MAKE SURE THERE IS NO SPACE IN THE SAMPLE IDs",
                self.input.display()
            );
        }
    }

    fn parse_header(&mut self, header_line: &str) {
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
        self.ntax = tax
            .parse::<usize>()
            .expect("HEADER TAXA NUMBER IS NOT A NUMBER");
        self.nchar = seq
            .parse::<usize>()
            .expect("HEADER CHARS LENGTH IS NOT A NUMBER");
    }

    fn check_ntax_matches(&self) {
        if self.matrix.len() != self.ntax {
            panic!(
                "ERROR READING PHYLIP FILE: {}. \
            THE NUMBER OF TAXA DOES NOT MATCH THE INFORMATION IN THE HEADER.\
            IN THE HEADER: {} \
            AND TAXA FOUND: {}",
                self.input.display(),
                self.ntax,
                self.matrix.len()
            );
        }
    }

    fn check_nchar_matches(&self, longest: usize) {
        if self.nchar != longest {
            panic!(
                "ERROR READING PHYLIP FILE {}, \
            THE NCHAR VALUE IN THE HEADER DOES NOT MATCH THE SEQUENCE LENGTH. \
            THE VALUE IN THE HEADER {}. \
            THE SEQUENCE LENGTH {}.",
                self.input.display(),
                self.nchar,
                longest
            );
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn read_phylip_simple_test() {
        let path = Path::new("test_files/simple.phy");
        let mut phylip = Phylip::new(path);
        phylip.read().unwrap();

        assert_eq!(2, phylip.ntax);
        assert_eq!(4, phylip.nchar);
        assert_eq!(2, phylip.matrix.len());
    }

    #[test]
    #[should_panic]
    fn read_phylip_invalid_test() {
        let path = Path::new("test_files/invalid.phy");
        let mut phylip = Phylip::new(path);
        phylip.read().unwrap();
    }

    #[test]
    fn read_phylip_whitespace_test() {
        let path = Path::new("test_files/whitespaces.phy");
        let mut phylip = Phylip::new(path);
        phylip.read().unwrap();
        assert_eq!(2, phylip.ntax);
        assert_eq!(4, phylip.nchar);
        assert_eq!(2, phylip.matrix.len());
    }

    #[test]
    fn parse_phylip_header_test() {
        let header = "2 24";
        let mut phy = Phylip::new(Path::new("."));
        phy.parse_header(header);

        assert_eq!(2, phy.ntax);
        assert_eq!(24, phy.nchar);
    }
}
