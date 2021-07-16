use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Read, Result};
use std::path::Path;

use indexmap::{IndexMap, IndexSet};
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

    pub fn read(&mut self) -> Result<()> {
        let mut commands = self.get_commands();
        self.parse_dimensions(&mut commands.dimensions);
        self.parse_format(&mut commands.format);
        self.parse_matrix(&mut commands.matrix);
        let mut seq_info = SeqCheck::new();
        seq_info.get_sequence_info(&self.matrix);
        self.is_alignment = seq_info.is_alignment;
        self.check_ntax_matches();
        self.check_nchar_matches(seq_info.longest);
        Ok(())
    }

    pub fn read_only_id(&mut self) -> IndexSet<String> {
        let mut commands = self.get_commands();
        self.parse_dimensions(&mut commands.dimensions);
        let ids = self.parse_matrix_id(&mut commands.matrix);
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

    fn get_commands(&mut self) -> Commands {
        let input = File::open(self.input).expect("CANNOT OPEN THE INPUT FILE");
        let mut buff = BufReader::new(input);
        let mut header = String::new();
        buff.read_line(&mut header)
            .expect("CANNOT READ THE HEADER FILE");
        self.check_nexus(&header.trim());
        self.parse_blocks(buff)
    }

    fn parse_blocks<R: Read>(&self, buff: R) -> Commands {
        let reader = NexusReader::new(buff);
        let mut commands = Commands::new();
        reader.into_iter().for_each(|read| {
            match read.to_lowercase() {
                command if command.starts_with("dimensions") => {
                    commands.dimensions.push_str(&read.trim().to_lowercase())
                }
                command if command.starts_with("format") => {
                    commands.format.push_str(&read.trim().to_lowercase())
                }
                command if command.starts_with("matrix") => commands.matrix.push_str(&read.trim()),
                _ => (),
            };
        });

        commands
    }

    fn parse_dimensions(&mut self, input: &mut String) {
        input.pop();
        let dimensions: Vec<&str> = input.split_whitespace().collect();
        dimensions
            .iter()
            .map(|d| d.trim())
            .for_each(|dimension| match dimension {
                tag if tag.starts_with("ntax") => self.header.ntax = self.parse_ntax(&dimension),
                tag if tag.starts_with("nchar") => {
                    self.header.nchar = self.parse_characters(&dimension)
                }
                _ => (),
            });
    }

    fn parse_format(&mut self, input: &mut String) {
        input.pop();
        let formats: Vec<&str> = input.split_whitespace().collect();
        formats
            .iter()
            .map(|f| f.trim())
            .filter(|f| !f.is_empty())
            .for_each(|format| match format {
                tag if tag.starts_with("datatype") => {
                    self.header.datatype = self.parse_datatype(&format)
                }
                tag if tag.starts_with("missing") => {
                    self.header.missing = self.parse_missing(&format)
                }
                tag if tag.starts_with("gap") => self.header.gap = self.parse_gap(&format),
                "interleave=yes" => self.interleave = true,
                "interleave" => self.interleave = true,
                _ => (),
            });
    }

    fn parse_matrix_id(&mut self, read: &mut String) -> IndexSet<String> {
        read.pop(); // remove terminated semicolon.
        let matrix: Vec<&str> = read.split('\n').collect();
        let mut ids = IndexSet::new();
        matrix[1..]
            .iter()
            .map(|l| l.trim())
            .filter(|l| !l.is_empty())
            .for_each(|line| {
                let seq: Vec<&str> = line.split_whitespace().collect();
                if seq.len() == 2 && !ids.contains(seq[0]) {
                    ids.insert(seq[0].to_string());
                }
            });
        read.clear();
        ids
    }

    fn parse_matrix(&mut self, read: &mut String) {
        read.pop(); // remove terminated semicolon.
        let matrix: Vec<&str> = read.split('\n').collect();
        matrix[1..]
            .iter()
            .map(|l| l.trim())
            .filter(|l| !l.is_empty())
            .for_each(|line| {
                let (id, dna) = self.parse_sequence(line);
                if self.interleave {
                    self.insert_matrix_interleave(id, dna);
                } else {
                    self.insert_matrix(id, dna);
                }
            });
        read.clear();
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

    fn parse_sequence(&self, line: &str) -> (String, String) {
        let seq: Vec<&str> = line.split_whitespace().collect();
        self.check_seq_len(seq.len());
        let id = seq[0].to_string();
        let dna = seq[1].to_string();
        common::check_valid_dna(&self.input, &id, &dna);
        (id, dna)
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

    fn check_seq_len(&self, len: usize) {
        if len != 2 {
            panic!(
                "THE FILE {} IS UNSUPPORTED NEXUS FORMAT. \
            MAKE SURE THERE IS NO SPACE IN THE SAMPLE IDs",
                self.input.display()
            );
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

struct Commands {
    matrix: String,
    dimensions: String,
    format: String,
}

impl Commands {
    fn new() -> Self {
        Self {
            matrix: String::new(),
            dimensions: String::new(),
            format: String::new(),
        }
    }
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

    fn next_block(&mut self) -> Option<String> {
        self.buffer.clear();
        let bytes = self
            .reader
            .read_until(b';', &mut self.buffer)
            .expect("ERROR READING FILES");
        if bytes == 0 {
            None
        } else {
            let mut token = String::new();
            if let Ok(tok) = std::str::from_utf8(&self.buffer) {
                token.push_str(tok.trim());
            }
            Some(token)
        }
    }
}

// Iterate over the file.
// Collect each of the nexus block terminated by semi-colon.
impl<R: Read> Iterator for NexusReader<R> {
    type Item = String;

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
        nex.read().unwrap();
        assert_eq!(1, nex.matrix.len());
    }

    #[test]
    fn nexus_reading_complete_test() {
        let sample = Path::new("test_files/complete.nex");
        let mut nex = Nexus::new(sample);
        nex.read().unwrap();
        assert_eq!(5, nex.matrix.len());
    }

    #[test]
    fn nexus_reading_tabulated_test() {
        let sample = Path::new("test_files/tabulated.nex");
        let mut nex = Nexus::new(sample);
        nex.read().unwrap();
        assert_eq!(2, nex.matrix.len());
    }

    #[test]
    fn nexus_parsing_object_test() {
        let sample = Path::new("test_files/complete.nex");
        let mut nex = Nexus::new(sample);
        nex.read().unwrap();
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
        nex.read().unwrap();
    }

    #[test]
    #[should_panic]
    fn check_match_ntax_panic_test() {
        let sample = Path::new("test_files/unmatched_block.nex");
        let mut nex = Nexus::new(sample);
        nex.read().unwrap();
    }

    #[test]
    #[should_panic]
    fn check_invalid_nexus_test() {
        let sample = Path::new("test_files/simple.fas");
        let mut nex = Nexus::new(sample);
        nex.read().unwrap();
    }

    #[test]
    #[should_panic]
    fn nexus_duplicate_panic_test() {
        let sample = Path::new("test_files/duplicates.nex");
        let mut nex = Nexus::new(sample);
        nex.read().unwrap();
    }

    #[test]
    #[should_panic]
    fn nexus_space_panic_test() {
        let sample = Path::new("test_files/idspaces.nex");
        let mut nex = Nexus::new(sample);
        nex.read().unwrap();
    }

    #[test]
    fn nexus_sequence_test() {
        let sample = Path::new("test_files/tabulated.nex");
        let mut nex = Nexus::new(sample);
        nex.read().unwrap();
        let key = String::from("ABEF");
        let res = String::from("gatata---");
        assert_eq!(Some(&res), nex.matrix.get(&key));
    }

    #[test]
    fn nexus_parse_interleave() {
        let sample = Path::new("test_files/interleave.nex");
        let mut nex = Nexus::new(sample);
        nex.read().unwrap();
        assert_eq!(3, nex.matrix.len());
    }

    #[test]
    fn nexus_parse_interleave_res_test() {
        let sample = Path::new("test_files/interleave.nex");
        let mut nex = Nexus::new(sample);
        nex.read().unwrap();
        let key = String::from("ABCD");
        let res = String::from("gatatagatatt");
        assert_eq!(Some(&res), nex.matrix.get(&key));
    }
}
