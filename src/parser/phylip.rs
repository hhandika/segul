use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Lines, Result};
use std::path::Path;

use indexmap::IndexMap;
use nom::{character::complete, sequence, IResult};

use crate::helper::common::{self, Header, SeqCheck};

pub struct Phylip<'a> {
    input: &'a Path,
    pub matrix: IndexMap<String, String>,
    pub header: Header,
    pub is_alignment: bool,
}

impl<'a> Phylip<'a> {
    pub fn new(input: &'a Path) -> Self {
        Self {
            input,
            matrix: IndexMap::new(),
            header: Header::new(),
            is_alignment: false,
        }
    }

    pub fn parse(&mut self) -> Result<()> {
        self.parse_matrix()?;
        let mut seq_info = SeqCheck::new();
        seq_info.get_sequence_info(&self.matrix);
        self.is_alignment = seq_info.is_alignment;
        self.check_ntax_matches();
        self.check_nchar_matches(seq_info.longest);
        Ok(())
    }

    pub fn parse_only_id(&mut self) -> Vec<String> {
        let file = File::open(self.input).expect("CANNOT READ THE FILE");
        let mut buff = BufReader::new(file);
        let mut header_line = String::new();
        buff.read_line(&mut header_line).unwrap();
        self.parse_header(&header_line.trim());
        let mut ids = Vec::new();
        buff.lines().filter_map(|ok| ok.ok()).for_each(|line| {
            let seq: Vec<&str> = line.split_whitespace().collect();
            if seq.len() == 2 {
                ids.push(seq[0].to_string());
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

    fn parse_matrix(&mut self) -> Result<()> {
        let file = File::open(self.input)?;
        let mut buff = BufReader::new(file);
        let mut header_line = String::new();
        buff.read_line(&mut header_line)?;
        self.parse_header(&header_line.trim());
        let records = Reader::new(buff, self.header.ntax);
        let mut ids: IndexMap<usize, String> = IndexMap::new();
        records.into_iter().for_each(|rec| match rec.id {
            Some(id) => {
                common::check_valid_dna(&self.input, &id, &rec.seq);
                ids.insert(rec.pos, id.clone());
                self.insert_matrix(id, rec.seq);
            }
            None => {
                if let Some(id) = ids.get(&rec.pos) {
                    if let Some(value) = self.matrix.get_mut(id) {
                        common::check_valid_dna(&self.input, &id, &rec.seq);
                        value.push_str(&rec.seq);
                    }
                }
            }
        });

        Ok(())
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

struct Records {
    id: Option<String>,
    seq: String,
    pos: usize,
}

impl Records {
    fn new() -> Self {
        Self {
            id: None,
            seq: String::new(),
            pos: 0,
        }
    }
}

struct Reader<R> {
    reader: Lines<BufReader<R>>,
    interleave: bool,
    pos: usize,
    ntax: usize,
}

impl<R: Read> Reader<R> {
    fn new(file: R, ntax: usize) -> Self {
        Self {
            reader: BufReader::new(file).lines(),
            interleave: false,
            pos: 1,
            ntax,
        }
    }

    fn next_seq(&mut self) -> Option<Records> {
        while let Some(Ok(lines)) = self.reader.by_ref().next() {
            let line = lines.trim();
            if !line.is_empty() {
                let mut records = Records::new();
                if !self.interleave {
                    let seq: Vec<&str> = line.split_whitespace().collect();
                    if seq.len() == 2 {
                        records.id = Some(seq[0].trim().to_string());
                        records.seq = seq[1].trim().to_string();
                    }
                } else {
                    records.id = None;
                    records.seq = line.to_string();
                }

                records.pos = self.pos;
                self.pos += 1;

                if self.pos == self.ntax + 1 {
                    self.pos = 1;
                    self.interleave = true;
                }

                return Some(records);
            }
            continue;
        }
        None
    }
}

impl<R: Read> Iterator for Reader<R> {
    type Item = Records;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_seq()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn read_phylip_simple_test() {
        let path = Path::new("test_files/simple.phy");
        let mut phylip = Phylip::new(path);
        phylip.parse().unwrap();

        assert_eq!(2, phylip.header.ntax);
        assert_eq!(4, phylip.header.nchar);
        assert_eq!(2, phylip.matrix.len());
    }

    #[test]
    #[should_panic]
    fn read_phylip_invalid_test() {
        let path = Path::new("test_files/invalid.phy");
        let mut phylip = Phylip::new(path);
        phylip.parse().unwrap();
    }

    #[test]
    fn read_phylip_whitespace_test() {
        let path = Path::new("test_files/whitespaces.phy");
        let mut phylip = Phylip::new(path);
        phylip.parse().unwrap();
        assert_eq!(2, phylip.header.ntax);
        assert_eq!(4, phylip.header.nchar);
        assert_eq!(2, phylip.matrix.len());
    }

    #[test]
    fn parse_phylip_header_test() {
        let header = "2 24";
        let mut phy = Phylip::new(Path::new("."));
        phy.parse_header(header);

        assert_eq!(2, phy.header.ntax);
        assert_eq!(24, phy.header.nchar);
    }

    #[test]
    fn parse_matrix_interleave_phylip_test() {
        let path = Path::new("test_files/interleave.phy");
        let mut phy = Phylip::new(path);
        phy.parse_matrix().unwrap();
        let res = phy.matrix.get("ABCD");
        assert_eq!(Some(&String::from("agccatggaa")), res);
    }

    #[test]
    fn read_interleave_phylip_test() {
        let path = Path::new("test_files/interleave.phy");
        let mut phy = Phylip::new(path);
        phy.parse().unwrap();
        let res = phy.matrix.get("ABCD");
        assert_eq!(Some(&String::from("agccatggaa")), res);
    }

    #[test]
    fn reade_int_phylip_whitespaces_test() {
        let path = Path::new("test_files/interleave_whitespaces.phy");
        let mut phy = Phylip::new(path);
        phy.parse().unwrap();
        let res = phy.matrix.get("ABCD");
        assert_eq!(Some(&String::from("agccatggcc")), res);
    }

    #[test]
    #[should_panic]
    fn read_interleave_phylip_panic_test() {
        // The header does not matches the character length.
        let path = Path::new("test_files/invalid_interleave.phy");
        let mut phy = Phylip::new(path);
        phy.parse().unwrap();
    }
}
