use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Result};
use std::path::Path;

use indexmap::{IndexMap, IndexSet};
use nom::{character::complete, sequence, IResult};

use crate::helper::common::{self, Header, SeqCheck};

pub struct Phylip<'a> {
    input: &'a Path,
    interleave: bool,
    pub matrix: IndexMap<String, String>,
    pub header: Header,
    pub is_alignment: bool,
}

impl<'a> Phylip<'a> {
    pub fn new(input: &'a Path, interleave: bool) -> Self {
        Self {
            input,
            interleave,
            matrix: IndexMap::new(),
            header: Header::new(),
            is_alignment: false,
        }
    }

    pub fn read(&mut self) -> Result<()> {
        self.read_file()?;
        let mut seq_info = SeqCheck::new();
        seq_info.get_sequence_info(&self.matrix);
        self.is_alignment = seq_info.is_alignment;
        self.check_ntax_matches();
        self.check_nchar_matches(seq_info.longest);
        Ok(())
    }

    pub fn read_only_id(&mut self) -> IndexSet<String> {
        let file = File::open(self.input).expect("CANNOT READ THE FILE");
        let mut buff = BufReader::new(file);
        let mut header_line = String::new();
        buff.read_line(&mut header_line).unwrap();
        self.parse_header(&header_line.trim());
        let mut ids = IndexSet::new();
        buff.lines().filter_map(|ok| ok.ok()).for_each(|line| {
            let seq: Vec<&str> = line.split_whitespace().collect();
            if seq.len() == 2 {
                ids.insert(seq[0].to_string());
            }
        });
        assert!(
            ids.len() == self.header.ntax,
            "FAILED PARSING {}. \
        THE NUMBER OF TAXA DOES NOT MATCH THE INFORMATION IN THE HEADER.\
            IN THE HEADER: {} \
            AND TAXA FOUND: {}
        ",
            self.input.display(),
            self.header.ntax,
            ids.len()
        );
        ids
    }

    fn read_file(&mut self) -> Result<()> {
        let file = File::open(self.input)?;

        if self.interleave {
            self.read_interleave(file)?;
        } else {
            self.read_sequential(file)?;
        }

        Ok(())
    }

    fn read_sequential<R: Read>(&mut self, file: R) -> Result<()> {
        let mut buff = BufReader::new(file);
        let mut header_line = String::new();
        buff.read_line(&mut header_line)?;
        self.parse_header(&header_line.trim());

        buff.lines()
            .filter_map(|ok| ok.ok())
            .filter(|l| !l.is_empty())
            .for_each(|line| {
                let (id, dna) = self.parse_sequence(line.trim());
                self.insert_matrix(id, dna);
            });

        Ok(())
    }

    fn read_interleave<R: Read>(&mut self, file: R) -> Result<()> {
        let mut buff = BufReader::new(file);
        let mut header_line = String::new();
        buff.read_line(&mut header_line)?;
        self.parse_header(&header_line.trim());
        let mut pos: usize = 1;
        let mut ids: IndexMap<usize, String> = IndexMap::new();
        let mut seq = false;
        buff.lines()
            .filter_map(|ok| ok.ok())
            .filter(|l| !l.is_empty())
            .for_each(|line| {
                if !seq {
                    // then, the line contains the sequence id.
                    let (id, dna) = self.parse_sequence(line.trim());
                    ids.insert(pos, id.clone());
                    self.insert_matrix(id, dna);
                    pos += 1;
                } else if let Some(id) = ids.get(&pos) {
                    if let Some(value) = self.matrix.get_mut(id) {
                        value.push_str(line.trim());
                        pos += 1;
                    }
                }

                // We reset the pos position after reaching
                // the end of the id lines.
                if pos == self.header.ntax + 1 {
                    pos = 1;
                    seq = true;
                }
            });

        Ok(())
    }

    fn parse_sequence(&self, line: &str) -> (String, String) {
        let seq: Vec<&str> = line.split_whitespace().collect();
        self.check_seq_len(seq.len());
        let id = seq[0].to_string();
        let dna = seq[1].to_string();
        common::check_valid_dna(&self.input, &id, &dna);
        (id, dna)
    }

    fn insert_matrix(&mut self, id: String, dna: String) {
        match self.matrix.get(&id) {
            Some(_) => panic!(
                "DUPLICATE SAMPLES FOR FILE {}. FIRST DUPLICATE FOUND: {}",
                self.input.display(),
                id
            ),
            None => {
                self.matrix.insert(id, dna);
            }
        }
    }

    fn check_seq_len(&self, len: usize) {
        if len != 2 {
            panic!(
                "CAN'T PARSE {}. \
            MAKE SURE THERE IS NO SPACE IN THE SAMPLE IDs. \
            OR YOUR FILES MAY BE IN INTERLEAVED FORMAT. \
            TRY USING --interleave FLAG. ",
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
        self.header.ntax = tax
            .parse::<usize>()
            .expect("HEADER TAXA NUMBER IS NOT A NUMBER");
        self.header.nchar = seq
            .parse::<usize>()
            .expect("HEADER CHARS LENGTH IS NOT A NUMBER");
    }

    fn check_ntax_matches(&self) {
        if self.matrix.len() != self.header.ntax {
            panic!(
                "ERROR READING PHYLIP FILE: {}. \
            THE NUMBER OF TAXA DOES NOT MATCH THE INFORMATION IN THE HEADER.\
            IN THE HEADER: {} \
            AND TAXA FOUND: {}",
                self.input.display(),
                self.header.ntax,
                self.matrix.len()
            );
        }
    }

    fn check_nchar_matches(&self, longest: usize) {
        if self.header.nchar != longest {
            panic!(
                "ERROR READING PHYLIP FILE {}, \
            THE NCHAR VALUE IN THE HEADER DOES NOT MATCH THE SEQUENCE LENGTH. \
            THE VALUE IN THE HEADER {}. \
            THE SEQUENCE LENGTH {}.",
                self.input.display(),
                self.header.nchar,
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
        let mut phylip = Phylip::new(path, false);
        phylip.read().unwrap();

        assert_eq!(2, phylip.header.ntax);
        assert_eq!(4, phylip.header.nchar);
        assert_eq!(2, phylip.matrix.len());
    }

    #[test]
    #[should_panic]
    fn read_phylip_invalid_test() {
        let path = Path::new("test_files/invalid.phy");
        let mut phylip = Phylip::new(path, false);
        phylip.read().unwrap();
    }

    #[test]
    fn read_phylip_whitespace_test() {
        let path = Path::new("test_files/whitespaces.phy");
        let mut phylip = Phylip::new(path, false);
        phylip.read().unwrap();
        assert_eq!(2, phylip.header.ntax);
        assert_eq!(4, phylip.header.nchar);
        assert_eq!(2, phylip.matrix.len());
    }

    #[test]
    fn parse_phylip_header_test() {
        let header = "2 24";
        let mut phy = Phylip::new(Path::new("."), false);
        phy.parse_header(header);

        assert_eq!(2, phy.header.ntax);
        assert_eq!(24, phy.header.nchar);
    }

    #[test]
    fn parse_interleave_phylip_test() {
        let path = Path::new("test_files/interleave.phy");
        let file = File::open(path).unwrap();
        let mut phy = Phylip::new(path, true);
        phy.read_interleave(file).unwrap();
        let res = phy.matrix.get("ABCD");
        assert_eq!(Some(&String::from("agccatggaa")), res);
    }

    #[test]
    fn read_interleave_phylip_test() {
        let path = Path::new("test_files/interleave.phy");
        let mut phy = Phylip::new(path, true);
        phy.read().unwrap();
        let res = phy.matrix.get("ABCD");
        assert_eq!(Some(&String::from("agccatggaa")), res);
    }

    #[test]
    fn reade_int_phylip_whitespaces_test() {
        let path = Path::new("test_files/interleave_whitespaces.phy");
        let mut phy = Phylip::new(path, true);
        phy.read().unwrap();
        let res = phy.matrix.get("ABCD");
        assert_eq!(Some(&String::from("agccatggcc")), res);
    }

    #[test]
    #[should_panic]
    fn read_interleave_phylip_panic_test() {
        // The header does not matches the character length.
        let path = Path::new("test_files/invalid_interleave.phy");
        let mut phy = Phylip::new(path, true);
        phy.read().unwrap();
    }
}