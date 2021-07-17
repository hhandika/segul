use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Read, Result};
use std::path::Path;

use indexmap::IndexMap;
use nom::{bytes::complete, character, sequence, IResult};

use crate::helper::common::{self, Header, SeqCheck};

pub struct Nexus<'a> {
    input: &'a Path,
    pub matrix: IndexMap<String, String>,
    pub header: Header,
    pub interleave: bool,
    pub is_alignment: bool,
}

impl<'a> Nexus<'a> {
    pub fn new(input: &'a Path) -> Self {
        Self {
            input,
            matrix: IndexMap::new(),
            header: Header::new(),
            interleave: false,
            is_alignment: false,
        }
    }

    pub fn parse(&mut self) -> Result<()> {
        let blocks = self.get_blocks();
        blocks.iter().for_each(|block| match block {
            Block::Dimensions(dimensions) => self.parse_dimensions(&dimensions),
            Block::Format(format) => self.parse_format(&format),
            Block::Matrix(matrix) => self.parse_matrix(&matrix),
            Block::Undetermined => (),
        });
        let mut seq_info = SeqCheck::new();
        seq_info.get_sequence_info(&self.matrix);
        self.is_alignment = seq_info.is_alignment;
        self.check_ntax_matches();
        self.check_nchar_matches(seq_info.longest);
        Ok(())
    }

    pub fn parse_only_id(&mut self) -> Vec<String> {
        let blocks = self.get_blocks();
        let mut ids = Vec::new();
        blocks.iter().for_each(|block| match block {
            Block::Dimensions(dimensions) => self.parse_dimensions(dimensions),
            Block::Matrix(matrix) => {
                matrix.iter().for_each(|(id, _)| {
                    ids.push(id.to_string());
                });
            }
            _ => (),
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

    fn get_blocks(&mut self) -> Vec<Block> {
        let input = File::open(self.input).expect("CANNOT OPEN THE INPUT FILE");
        let mut buff = BufReader::new(input);
        let mut header = String::new();
        buff.read_line(&mut header)
            .expect("CANNOT READ THE HEADER FILE");
        self.check_nexus(&header.trim());
        let reader = NexusReader::new(buff);
        reader.into_iter().collect()
    }

    fn parse_dimensions(&mut self, blocks: &[String]) {
        blocks.iter().for_each(|dimension| match dimension {
            tag if tag.starts_with("ntax") => self.header.ntax = self.parse_ntax(&dimension),
            tag if tag.starts_with("nchar") => {
                self.header.nchar = self.parse_characters(&dimension)
            }
            _ => (),
        });
    }

    fn parse_format(&mut self, blocks: &[String]) {
        blocks.iter().for_each(|format| match format {
            token if token.starts_with("datatype") => {
                self.header.datatype = self.parse_datatype(&format)
            }
            token if token.starts_with("missing") => {
                self.header.missing = self.parse_missing(&format)
            }
            token if token.starts_with("gap") => self.header.gap = self.parse_gap(&format),
            token if token.starts_with("interleave") => self.parse_interleave(&format),
            _ => (),
        });
    }

    fn parse_matrix(&mut self, matrix: &[(String, String)]) {
        matrix.iter().for_each(|(id, seq)| {
            common::check_valid_dna(&self.input, &id, &seq);
            if self.interleave {
                self.insert_matrix_interleave(id.to_string(), seq.to_string());
            } else {
                self.insert_matrix(id.to_string(), seq.to_string());
            }
        });
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

    fn insert_matrix_interleave(&mut self, id: String, dna: String) {
        match self.matrix.get_mut(&id) {
            Some(value) => value.push_str(&dna),
            None => {
                self.matrix.insert(id, dna);
            }
        }
    }

    fn parse_interleave(&mut self, tokens: &str) {
        match tokens {
            "interleave=yes" => self.interleave = true,
            "interleave" => self.interleave = true,
            "interleave=no" => self.interleave = true,
            _ => (),
        }
    }

    fn parse_datatype(&self, input: &str) -> String {
        let tag: IResult<&str, &str> =
            sequence::preceded(complete::tag("datatype="), character::complete::alpha1)(input);
        self.parse_string(tag)
    }

    fn parse_gap(&self, input: &str) -> char {
        let gap = input.replace("gap=", "");
        self.parse_char(&gap)
    }

    fn parse_missing(&self, input: &str) -> char {
        let missing = input.replace("missing=", "");
        self.parse_char(&missing)
    }

    fn parse_string(&self, tag: IResult<&str, &str>) -> String {
        let mut text = String::new();
        self.convert_nomtag_to_string(tag, &mut text);
        text
    }

    fn parse_char(&self, text: &str) -> char {
        text.parse::<char>().expect("CANNOT PARSE TAG TO CHAR")
    }

    fn parse_ntax(&self, input: &str) -> usize {
        let tag: IResult<&str, &str> =
            sequence::preceded(complete::tag("ntax="), character::complete::digit1)(input);
        self.parse_usize(tag)
    }

    fn parse_characters(&self, input: &str) -> usize {
        let tag: IResult<&str, &str> =
            sequence::preceded(complete::tag("nchar="), character::complete::digit1)(input);
        self.parse_usize(tag)
    }

    fn parse_usize(&self, tag: IResult<&str, &str>) -> usize {
        let mut text = String::new();
        self.convert_nomtag_to_string(tag, &mut text);
        text.parse::<usize>()
            .expect("HEADER TAXA NUMBER IS NOT A NUMBER")
    }

    fn convert_nomtag_to_string(&self, tag: IResult<&str, &str>, text: &mut String) {
        match tag {
            Ok((_, out)) => text.push_str(out.trim()),
            Err(_) => eprintln!("CANNOT PARSE NEXUS TAG"),
        }
    }

    fn check_nexus(&self, line: &str) {
        if !line.to_lowercase().starts_with("#nexus") {
            panic!("THE FILE {} IS INVALID NEXUS FORMAT", self.input.display());
        }
    }

    fn check_ntax_matches(&self) {
        if self.matrix.len() != self.header.ntax {
            panic!(
                "ERROR READING NEXUS FILE: {}. \
            THE NUMBER OF TAXA DOES NOT MATCH THE INFORMATION IN THE BLOCK.\
            IN THE BLOCK: {} \
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
                "ERROR READING NEXUS FILE {}, \
            THE NCHAR VALUE IN THE BLOCK DOES NOT MATCH THE SEQUENCE LENGTH. \
            THE VALUE IN THE BLOCK {}. \
            THE SEQUENCE LENGTH {}.",
                self.input.display(),
                self.header.nchar,
                longest
            );
        }
    }
}

enum Block {
    Dimensions(Vec<String>),
    Format(Vec<String>),
    Matrix(Vec<(String, String)>),
    Undetermined,
}

struct NexusReader<R> {
    reader: BufReader<R>,
    buffer: Vec<u8>,
}

impl<R: Read> NexusReader<R> {
    fn new(file: R) -> Self {
        Self {
            reader: BufReader::new(file),
            buffer: Vec::new(),
        }
    }

    fn next_block(&mut self) -> Option<Block> {
        self.buffer.clear();
        let bytes = self
            .reader
            .read_until(b';', &mut self.buffer)
            .expect("ERROR READING FILES");
        if bytes == 0 {
            None
        } else {
            let mut block: String = std::str::from_utf8(&self.buffer)
                .expect("Failed parsing nexus block")
                .trim()
                .to_string();
            block.pop(); // remove terminated semicolon
            match block.to_lowercase() {
                b if b.starts_with("dimensions") => {
                    Some(Block::Dimensions(self.parse_header(&block)))
                }
                b if b.starts_with("format") => Some(Block::Format(self.parse_header(&block))),
                b if b.starts_with("matrix") => Some(Block::Matrix(self.parse_matrix(&block))),
                _ => Some(Block::Undetermined),
            }
        }
    }

    fn parse_header(&self, block: &str) -> Vec<String> {
        let headers: Vec<&str> = block.split_whitespace().collect();
        let mut tokens: Vec<String> = Vec::new();
        headers[1..]
            .iter()
            .filter(|h| !h.is_empty())
            .for_each(|h| tokens.push(h.to_lowercase()));
        tokens
    }

    fn parse_matrix(&self, block: &str) -> Vec<(String, String)> {
        let matrix: Vec<&str> = block.split('\n').collect();
        let mut sequence = Vec::new();
        matrix[1..].iter().filter(|s| !s.is_empty()).for_each(|s| {
            let seq: Vec<&str> = s.split_whitespace().collect();
            if seq.len() == 2 {
                sequence.push((seq[0].to_string(), seq[1].to_string()));
            }
        });

        sequence
    }
}

// Iterate over the file.
// Collect each of the nexus block terminated by semi-colon.
impl<R: Read> Iterator for NexusReader<R> {
    type Item = Block;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_block()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn nexus_reading_simple_test() {
        let sample = Path::new("test_files/simple.nex");
        let mut nex = Nexus::new(sample);
        nex.parse().unwrap();
        assert_eq!(1, nex.matrix.len());
    }

    #[test]
    fn nexus_reading_complete_test() {
        let sample = Path::new("test_files/complete.nex");
        let mut nex = Nexus::new(sample);
        nex.parse().unwrap();
        assert_eq!(5, nex.matrix.len());
    }

    #[test]
    fn nexus_reading_tabulated_test() {
        let sample = Path::new("test_files/tabulated.nex");
        let mut nex = Nexus::new(sample);
        nex.parse().unwrap();
        assert_eq!(2, nex.matrix.len());
    }

    #[test]
    fn nexus_parsing_object_test() {
        let sample = Path::new("test_files/complete.nex");
        let mut nex = Nexus::new(sample);
        nex.parse().unwrap();
        assert_eq!(5, nex.header.ntax);
        assert_eq!(802, nex.header.nchar);
        assert_eq!("dna", nex.header.datatype);
        assert_eq!('-', nex.header.gap);
        assert_eq!('?', nex.header.missing);
    }

    #[test]
    fn nexus_parse_ntax_test() {
        let sample = Path::new(".");
        let tax = "ntax=5";
        let nex = Nexus::new(sample);
        let res = nex.parse_ntax(tax);
        assert_eq!(5, res);
    }

    #[test]
    // #[should_panic]
    fn check_match_ntax_test() {
        let sample = Path::new("test_files/simple.nex");
        let mut nex = Nexus::new(sample);
        nex.parse().unwrap();
    }

    #[test]
    #[should_panic]
    fn check_match_ntax_panic_test() {
        let sample = Path::new("test_files/unmatched_block.nex");
        let mut nex = Nexus::new(sample);
        nex.parse().unwrap();
    }

    #[test]
    #[should_panic]
    fn check_invalid_nexus_test() {
        let sample = Path::new("test_files/simple.fas");
        let mut nex = Nexus::new(sample);
        nex.parse().unwrap();
    }

    #[test]
    #[should_panic]
    fn nexus_duplicate_panic_test() {
        let sample = Path::new("test_files/duplicates.nex");
        let mut nex = Nexus::new(sample);
        nex.parse().unwrap();
    }

    #[test]
    #[should_panic]
    fn nexus_space_panic_test() {
        let sample = Path::new("test_files/idspaces.nex");
        let mut nex = Nexus::new(sample);
        nex.parse().unwrap();
    }

    #[test]
    fn nexus_sequence_test() {
        let sample = Path::new("test_files/tabulated.nex");
        let mut nex = Nexus::new(sample);
        nex.parse().unwrap();
        let key = String::from("ABEF");
        let res = String::from("gatata---");
        assert_eq!(Some(&res), nex.matrix.get(&key));
    }

    #[test]
    fn nexus_parse_interleave() {
        let sample = Path::new("test_files/interleave.nex");
        let mut nex = Nexus::new(sample);
        nex.parse().unwrap();
        assert_eq!(3, nex.matrix.len());
    }

    #[test]
    fn nexus_parse_interleave_res_test() {
        let sample = Path::new("test_files/interleave.nex");
        let mut nex = Nexus::new(sample);
        nex.parse().unwrap();
        let key = String::from("ABCD");
        let res = String::from("gatatagatatt");
        assert_eq!(Some(&res), nex.matrix.get(&key));
    }
}
