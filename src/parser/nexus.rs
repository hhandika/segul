use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Read};
use std::path::Path;

use ansi_term::Colour::Red;
use indexmap::{IndexMap, IndexSet};
use lazy_static::lazy_static;
use nom::{bytes::complete, character, sequence, IResult};
use regex::Regex;

use crate::helper::alphabet;
use crate::helper::sequence::SeqCheck;
use crate::helper::types::{DataType, Header, SeqMatrix};
use crate::parser;

pub struct Nexus<'a> {
    input: &'a Path,
    datatype: &'a DataType,
    pub matrix: SeqMatrix,
    pub header: Header,
    pub interleave: bool,
}

impl<'a> Nexus<'a> {
    pub fn new(input: &'a Path, datatype: &'a DataType) -> Self {
        Self {
            input,
            datatype,
            matrix: IndexMap::new(),
            header: Header::new(),
            interleave: false,
        }
    }

    pub fn parse(&mut self) {
        let blocks = self.get_blocks();
        self.parse_blocks(&blocks);
        let mut seq_info = SeqCheck::new();
        seq_info.check(&self.matrix);
        self.header.aligned = seq_info.is_alignment;
        self.check_ntax_matches();
        self.check_nchar_matches(seq_info.longest);
    }

    pub fn parse_only_id(&mut self) -> IndexSet<String> {
        let blocks = self.get_blocks();
        let mut ids = IndexSet::new();
        blocks.iter().for_each(|block| match block {
            Block::Dimensions(dimensions) => self.parse_dimensions(dimensions),
            Block::Matrix(matrix) => {
                for (id, _) in matrix.iter() {
                    if !ids.contains(id) {
                        ids.insert(id.to_string());
                    } else {
                        break;
                    }
                }
            }
            _ => (),
        });

        parser::warn_duplicate_ids!(self, ids);

        ids
    }

    fn get_blocks(&mut self) -> Vec<Block> {
        let input = File::open(self.input).expect("Failed opening nexus file");
        let mut buff = BufReader::new(input);
        let mut header = String::new();
        buff.read_line(&mut header)
            .unwrap_or_else(|_| panic!("Failed reading nexus header for {}", self.input.display()));
        self.check_nexus(header.trim());
        let reader = NexusReader::new(buff);
        reader.into_iter().collect()
    }

    fn parse_blocks(&mut self, blocks: &[Block]) {
        blocks.iter().for_each(|block| match block {
            Block::Dimensions(dimensions) => self.parse_dimensions(dimensions),
            Block::Format(format) => self.parse_format(format),
            Block::Matrix(matrix) => self.parse_matrix(matrix),
            Block::Undetermined => (),
        });
    }

    fn parse_dimensions(&mut self, blocks: &[String]) {
        blocks.iter().for_each(|dimension| match dimension {
            token if token.starts_with("ntax") => self.header.ntax = self.parse_ntax(dimension),
            token if token.starts_with("nchar") => {
                self.header.nchar = self.parse_characters(dimension)
            }
            _ => (),
        });
    }

    fn parse_format(&mut self, blocks: &[String]) {
        blocks.iter().for_each(|format| match format {
            token if token.starts_with("datatype") => {
                self.header.datatype = self.parse_datatype(format)
            }
            token if token.starts_with("missing") => {
                self.header.missing = self.parse_missing(format)
            }
            token if token.starts_with("gap") => self.header.gap = self.parse_gap(format),
            token if token.starts_with("interleave") => self.parse_interleave(format),
            _ => (),
        });
    }

    fn parse_matrix(&mut self, matrix: &[(String, String)]) {
        self.matrix.reserve(self.header.ntax);
        matrix.iter().for_each(|(id, seq)| {
            alphabet::check_valid_seq(self.input, self.datatype, id, seq);
            if self.interleave {
                self.insert_matrix_interleave(id.to_string(), seq.to_string());
            } else {
                parser::insert_matrix!(self, id, seq);
            }
        });
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
            "interleave=no" => self.interleave = false,
            _ => (),
        }
    }

    #[inline]
    fn parse_datatype(&self, input: &str) -> String {
        let tag: IResult<&str, &str> =
            sequence::preceded(complete::tag("datatype="), character::complete::alpha1)(input);
        self.parse_string(tag)
    }

    #[inline]
    fn parse_gap(&self, input: &str) -> char {
        let gap = input.replace("gap=", "");
        self.parse_char(&gap)
    }

    #[inline]
    fn parse_missing(&self, input: &str) -> char {
        let missing = input.replace("missing=", "");
        self.parse_char(&missing)
    }

    #[inline]
    fn parse_string(&self, tag: IResult<&str, &str>) -> String {
        let mut text = String::new();
        self.convert_nomtag_to_string(tag, &mut text);
        text
    }

    #[inline]
    fn parse_char(&self, text: &str) -> char {
        text.parse::<char>()
            .expect("Gaps or missing tags are not a char")
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
        text.parse::<usize>().expect("Header taxa is not a number")
    }

    fn convert_nomtag_to_string(&self, tag: IResult<&str, &str>, text: &mut String) {
        match tag {
            Ok((_, out)) => text.push_str(out.trim()),
            Err(_) => eprintln!("Failed parsing nexus header"),
        }
    }

    #[inline]
    fn check_nexus(&self, line: &str) {
        if !line.to_lowercase().starts_with("#nexus") {
            panic!("The file {} is invalid nexus format.", self.input.display());
        }
    }

    fn check_ntax_matches(&self) {
        if self.matrix.len() != self.header.ntax {
            panic!(
                "Error reading nexus file: {}. \
            The number of taxa does not match the information in the block.\
            In the block: {} \
            and taxa found: {}",
                self.input.display(),
                self.header.ntax,
                self.matrix.len()
            );
        }
    }

    fn check_nchar_matches(&self, longest: usize) {
        if self.header.nchar != longest {
            panic!(
                "Error reading nexus file {}, \
            the NCHAR value in the header does not match the sequence length. \
            the value in the block {}. \
            the sequence length {}.",
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
            .expect("Failed reading nexus file");
        if bytes == 0 {
            None
        } else {
            let mut block: String = std::str::from_utf8(&self.buffer)
                .expect("Failed parsing nexus block")
                .trim()
                .to_string();
            block.pop(); // remove terminated semicolon
            let commands = get_commands(&block);
            match commands.as_str() {
                "dimensions" => Some(Block::Dimensions(self.parse_header(&block))),
                "format" => Some(Block::Format(self.parse_header(&block))),
                "matrix" => Some(Block::Matrix(self.parse_matrix(&block))),
                _ => Some(Block::Undetermined),
            }
        }
    }

    fn parse_header(&self, block: &str) -> Vec<String> {
        let headers: Vec<&str> = block.split_whitespace().collect();
        let mut tokens: Vec<String> = Vec::new();
        headers[1..] // ignore NEXUS commands
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

fn get_commands(text: &str) -> String {
    lazy_static! { // Match the first word in the block
        static ref RE: Regex = Regex::new(r"^(\w+)").expect("Failed capturing nexus commands");
    }

    match RE.captures(text) {
        Some(word) => word[0].to_lowercase(),
        None => String::from(""),
    }
}

#[cfg(test)]
mod test {
    use super::*;

    const DNA: DataType = DataType::Dna;
    const AA: DataType = DataType::Aa;

    #[test]
    fn test_nexus_reading_simple() {
        let sample = Path::new("tests/files/simple.nex");
        let mut nex = Nexus::new(sample, &DNA);
        nex.parse();
        assert_eq!(1, nex.matrix.len());
    }

    #[test]
    fn test_nexus_reading_complete() {
        let sample = Path::new("tests/files/complete.nex");
        let mut nex = Nexus::new(sample, &DNA);
        nex.parse();
        assert_eq!(5, nex.matrix.len());
    }

    #[test]
    fn test_nexus_reading_tabulated() {
        let sample = Path::new("tests/files/tabulated.nex");
        let mut nex = Nexus::new(sample, &DNA);
        nex.parse();
        assert_eq!(2, nex.matrix.len());
    }

    #[test]
    fn test_nexus_parsing_object() {
        let sample = Path::new("tests/files/complete.nex");
        let mut nex = Nexus::new(sample, &DNA);
        nex.parse();
        assert_eq!(5, nex.header.ntax);
        assert_eq!(802, nex.header.nchar);
        assert_eq!("dna", nex.header.datatype);
        assert_eq!('-', nex.header.gap);
        assert_eq!('?', nex.header.missing);
    }

    #[test]
    fn test_nexus_parse_ntax() {
        let sample = Path::new(".");
        let tax = "ntax=5";
        let nex = Nexus::new(sample, &DNA);
        let res = nex.parse_ntax(tax);
        assert_eq!(5, res);
    }

    #[test]
    fn test_check_match_ntax() {
        let sample = Path::new("tests/files/simple.nex");
        let mut nex = Nexus::new(sample, &DNA);
        nex.parse();
    }

    #[test]
    #[should_panic]
    fn test_check_match_ntax_panic() {
        let sample = Path::new("tests/files/unmatched_block.nex");
        let mut nex = Nexus::new(sample, &DNA);
        nex.parse();
    }

    #[test]
    #[should_panic]
    fn test_check_invalid_nexus() {
        let sample = Path::new("tests/files/simple.fas");
        let mut nex = Nexus::new(sample, &DNA);
        nex.parse();
    }

    #[test]
    #[should_panic]
    fn test_nexus_duplicate_panic() {
        let sample = Path::new("tests/files/duplicates.nex");
        let mut nex = Nexus::new(sample, &DNA);
        nex.parse();
    }

    #[test]
    #[should_panic]
    fn test_nexus_space_panic() {
        let sample = Path::new("tests/files/idspaces.nex");
        let mut nex = Nexus::new(sample, &DNA);
        nex.parse();
    }

    #[test]
    fn test_nexus_sequence() {
        let sample = Path::new("tests/files/tabulated.nex");
        let mut nex = Nexus::new(sample, &DNA);
        nex.parse();
        let key = String::from("ABEF");
        let res = String::from("gatata---");
        assert_eq!(Some(&res), nex.matrix.get(&key));
    }

    #[test]
    fn test_nexus_parse_interleave() {
        let sample = Path::new("tests/files/interleave.nex");
        let mut nex = Nexus::new(sample, &DNA);
        nex.parse();
        assert_eq!(3, nex.matrix.len());
    }

    #[test]
    fn test_nexus_parse_interleave_res() {
        let sample = Path::new("tests/files/interleave.nex");
        let mut nex = Nexus::new(sample, &DNA);
        nex.parse();
        let key = String::from("ABCD");
        let res = String::from("gatatagatatt");
        assert_eq!(Some(&res), nex.matrix.get(&key));
    }

    #[test]
    fn test_nexus_simple_aa() {
        let sample = Path::new("tests/files/simple_aa.nex");
        let mut nex = Nexus::new(sample, &AA);
        nex.parse();
        let key = String::from("ABCE");
        let res = String::from("MAYPMQLGFQDATSPI");
        assert_eq!(Some(&res), nex.matrix.get(&key));
    }

    #[test]
    fn test_regex_command() {
        let text = "Matrix\n ABCD AGTC";
        assert_eq!(String::from("matrix"), get_commands(text));
    }
}
