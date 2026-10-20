//! Nexus parser, support sequential and interleaved format.
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use std::path::Path;

use colored::Colorize;
use indexmap::{IndexMap, IndexSet};
use lazy_static::lazy_static;
use nom::{IResult, Parser, bytes::complete, character, sequence};
use regex::Regex;

use crate::helper::alphabet;
use crate::helper::sequence::SeqCheck;
use crate::helper::types::{DataType, Header, SeqMatrix};
use crate::parser;

/// Parse a nexus file and return a sequence matrix.
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
        self.read_matrix(false);
        let mut seq_info = SeqCheck::new();
        seq_info.check(&self.matrix);
        self.header.aligned = seq_info.is_alignment;
        self.check_ntax_matches();
        self.check_nchar_matches(seq_info.longest);
    }

    pub fn parse_only_id(&mut self) -> IndexSet<String> {
        let ids = self.read_matrix(true);
        parser::warn_duplicate_ids!(self, ids);
        ids
    }

    /// Stream the file line by line. Header commands are buffered
    /// until their terminating semicolon (they are small),
    /// matrix lines are inserted directly into the matrix
    /// without keeping a copy of the whole block in memory.
    /// When `only_id` is true, sequences are never stored
    /// and the IDs are returned instead.
    fn read_matrix(&mut self, only_id: bool) -> IndexSet<String> {
        let file = File::open(self.input).expect("Failed opening nexus file");
        let mut reader = BufReader::new(file);
        let mut line = Vec::new();
        reader
            .read_until(b'\n', &mut line)
            .unwrap_or_else(|_| panic!("Failed reading nexus header for {}", self.input.display()));
        self.check_nexus(self.line_to_str(&line).trim());

        let mut ids = IndexSet::new();
        let mut command = String::new();
        let mut in_matrix = false;
        loop {
            line.clear();
            let bytes = reader
                .read_until(b'\n', &mut line)
                .expect("Failed reading nexus file");
            if bytes == 0 {
                break;
            }
            let text = self.line_to_str(&line);
            let done = if in_matrix {
                self.read_matrix_line(text, only_id, &mut ids)
            } else {
                in_matrix = self.read_header_line(text, &mut command);
                // Data that follows the matrix command on the same line.
                in_matrix && self.read_matrix_line(&command, only_id, &mut ids)
            };
            if done {
                break;
            }
        }

        ids
    }

    /// Accumulate header commands until a semicolon is found.
    /// Returns true when the matrix command starts.
    fn read_header_line(&mut self, text: &str, command: &mut String) -> bool {
        let mut rest = text;
        while let Some(pos) = rest.find(';') {
            command.push_str(&rest[..pos]);
            if get_commands(command.trim()) == "matrix" {
                // Matrix terminated on the same line, e.g. `matrix ABCD ACGT;`
                self.strip_matrix_command(command);
                command.push(';');
                return true;
            }
            self.parse_command(command.trim());
            command.clear();
            rest = &rest[pos + 1..];
        }
        command.push_str(rest);
        if get_commands(command.trim()) == "matrix" {
            self.strip_matrix_command(command);
            return true;
        }
        false
    }

    /// Keep only the text after the matrix command.
    fn strip_matrix_command(&mut self, command: &mut String) {
        let trimmed = command.trim_start();
        let data = trimmed["matrix".len()..].to_string();
        *command = data;
        self.matrix.reserve(self.header.ntax);
    }

    fn parse_command(&mut self, command: &str) {
        match get_commands(command).as_str() {
            "dimensions" => {
                let tokens = parse_header_tokens(command);
                self.parse_dimensions(&tokens)
            }
            "format" => {
                let tokens = parse_header_tokens(command);
                self.parse_format(&tokens)
            }
            _ => (),
        }
    }

    /// Parse a matrix line. Returns true when the matrix ends
    /// or, for ID-only parsing, when IDs start repeating.
    fn read_matrix_line(&mut self, text: &str, only_id: bool, ids: &mut IndexSet<String>) -> bool {
        let (data, end) = match text.find(';') {
            Some(pos) => (&text[..pos], true),
            None => (text, false),
        };
        if let Some((id, seq)) = parser::split_id_seq(data) {
            if only_id {
                if ids.contains(id) {
                    return true;
                }
                ids.insert(id.to_string());
            } else {
                alphabet::check_valid_seq(self.input, self.datatype, id, seq);
                if self.interleave {
                    self.insert_matrix_interleave(id, seq);
                } else {
                    parser::insert_matrix!(self, id, seq);
                }
            }
        }
        end
    }

    fn line_to_str<'b>(&self, line: &'b [u8]) -> &'b str {
        std::str::from_utf8(line)
            .unwrap_or_else(|_| panic!("Failed parsing nexus file {}", self.input.display()))
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

    fn insert_matrix_interleave(&mut self, id: &str, seq: &str) {
        match self.matrix.get_mut(id) {
            Some(value) => value.push_str(seq),
            None => {
                // Reserve the full length upfront to avoid
                // growing (and over-allocating) the string.
                let mut value = String::with_capacity(self.header.nchar.min(parser::MAX_RESERVE));
                value.push_str(seq);
                self.matrix.insert(id.to_string(), value);
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
            sequence::preceded(complete::tag("datatype="), character::complete::alpha1)
                .parse(input);
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
            sequence::preceded(complete::tag("ntax="), character::complete::digit1).parse(input);
        self.parse_usize(tag)
    }

    fn parse_characters(&self, input: &str) -> usize {
        let tag: IResult<&str, &str> =
            sequence::preceded(complete::tag("nchar="), character::complete::digit1).parse(input);
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

/// Split a header command into lowercase tokens,
/// ignoring the command itself.
fn parse_header_tokens(command: &str) -> Vec<String> {
    command
        .split_whitespace()
        .skip(1)
        .map(|token| token.to_lowercase())
        .collect()
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
    fn test_nexus_parse_only_id_interleave() {
        let sample = Path::new("tests/files/interleave.nex");
        let mut nex = Nexus::new(sample, &DNA);
        let ids = nex.parse_only_id();
        assert_eq!(
            vec!["ABCD", "ABCE", "ABCF"],
            ids.into_iter().collect::<Vec<_>>()
        );
        assert!(nex.matrix.is_empty());
    }

    #[test]
    fn test_nexus_parse_only_id_terminated_line() {
        let sample = Path::new("tests/files/tabulated.nex");
        let mut nex = Nexus::new(sample, &DNA);
        let ids = nex.parse_only_id();
        assert_eq!(2, ids.len());
        assert_eq!(2, nex.header.ntax);
    }

    #[test]
    fn test_regex_command() {
        let text = "Matrix\n ABCD AGTC";
        assert_eq!(String::from("matrix"), get_commands(text));
    }
}
