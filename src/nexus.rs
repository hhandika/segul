use std::collections::BTreeMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Lines, Read, Result};

use nom::{bytes::complete, character, sequence, IResult};

use crate::common::{self, SeqFormat};
use crate::writer::SeqWriter;

pub fn convert_nexus(path: &str, filetype: SeqFormat) {
    let mut nex = Nexus::new();
    nex.read(path).expect("CANNOT READ NEXUS FILES");
    let mut convert = SeqWriter::new(
        path,
        &nex.matrix,
        Some(nex.ntax),
        Some(nex.nchar),
        Some(nex.datatype),
        Some(nex.missing),
        Some(nex.gap),
    );

    match filetype {
        SeqFormat::Phylip => convert.write_sequence(&filetype),
        SeqFormat::Fasta => convert.write_fasta(),
        _ => (),
    }
}

pub struct Nexus {
    pub matrix: BTreeMap<String, String>,
    pub ntax: usize,
    pub nchar: usize,
    pub datatype: String,
    pub missing: char,
    pub gap: char,
}

impl Nexus {
    pub fn new() -> Self {
        Self {
            matrix: BTreeMap::new(),
            ntax: 0,
            nchar: 0,
            datatype: String::new(),
            missing: '?',
            gap: '-',
        }
    }

    pub fn read(&mut self, path: &str) -> Result<()> {
        let input = File::open(path).expect("CANNOT OPEN THE INPUT FILE");
        let mut buff = BufReader::new(input);
        let mut header = String::new();
        buff.read_line(&mut header)?;
        self.check_nexus(&header.trim());
        let mut commands = self.parse_blocks(buff);
        self.parse_matrix(&mut commands.matrix);
        self.parse_dimensions(&mut commands.dimensions);
        self.parse_format(&mut commands.format);
        self.check_ntax_matches();
        Ok(())
    }

    fn check_nexus(&self, line: &str) {
        if !line.to_lowercase().starts_with("#nexus") {
            panic!("INVALID NEXUS FORMAT");
        }
    }

    fn check_ntax_matches(&self) {
        if self.matrix.len() != self.ntax {
            panic!(
                "ERROR READING NEXUS FILES. \
            THE NUMBER OF TAXA IS NOT MATCH THE INFORMATION IN THE BLOCK.\
            IN THE BLOCK: {} \
            AND TAXA FOUND: {}",
                self.ntax,
                self.matrix.len()
            );
        }
    }

    fn parse_blocks<R: Read>(&self, buff: R) -> Commands {
        let reader = Reader::new(buff);
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
                tag if tag.starts_with("ntax") => self.ntax = self.parse_ntax(&dimension),
                tag if tag.starts_with("nchar") => self.nchar = self.parse_characters(&dimension),
                _ => (),
            });
    }

    fn parse_format(&mut self, input: &mut String) {
        input.pop();
        let formats: Vec<&str> = input.split_ascii_whitespace().collect();
        formats
            .iter()
            .map(|f| f.trim())
            .filter(|f| !f.is_empty())
            .for_each(|format| match format {
                tag if tag.starts_with("datatype") => self.datatype = self.parse_datatype(&format),
                tag if tag.starts_with("missing") => self.missing = self.parse_missing(&format),
                tag if tag.starts_with("gap") => self.gap = self.parse_gap(&format),
                _ => (),
            });
    }

    // Iterate each matrix. Convert dna to lowercase for
    // consistency with the other blocks lettercase.
    // This also matches biopython format.
    // The ID text is left intact.
    fn parse_matrix(&mut self, read: &mut String) {
        read.pop(); // remove terminated semicolon.
        let matrix: Vec<&str> = read.split('\n').collect();
        matrix[1..]
            .iter()
            .map(|l| l.trim())
            .filter(|l| !l.is_empty())
            .for_each(|line| {
                let seq: Vec<&str> = line.split_whitespace().collect();
                self.check_seq_len(seq.len());
                let id = seq[0].to_string();
                let dna = seq[1].to_string().to_lowercase();
                self.check_valid_dna(&id, &dna);
                #[allow(clippy::all)]
                if self.matrix.contains_key(&id) {
                    panic!("DUPLICATE SAMPLES. FIRST DUPLICATE FOUND: {}", id);
                } else {
                    self.matrix.insert(id, dna);
                }
            });
        read.clear();
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

    fn check_valid_dna(&self, id: &str, dna: &str) {
        if !common::is_valid_dna(dna) {
            panic!("INVALID DNA SEQUENCE FOUND FOR {}", id);
        }
    }

    fn check_seq_len(&self, len: usize) {
        if len != 2 {
            panic!(
                "UNSUPPORTED NEXUS FORMAT. \
            MAKE SURE THERE IS NO SPACE IN THE SAMPLE IDs"
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

struct Reader<R> {
    reader: Lines<BufReader<R>>,
    buffer: String,
    content: String,
}

impl<R: Read> Reader<R> {
    fn new(file: R) -> Self {
        Self {
            reader: BufReader::new(file).lines(),
            buffer: String::new(),
            content: String::new(),
        }
    }
}

// Iterate over the file.
// Collect each of the nexus block terminated by semi-colon.
impl<R: Read> Iterator for Reader<R> {
    type Item = String;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(Ok(line)) = self.reader.next() {
            self.buffer.push_str(&line);
            if !line.is_empty() {
                self.buffer.push('\n');
            }
            if line.ends_with(';') {
                self.content.push_str(&self.buffer);
                self.buffer.clear();
            }
            let token = self.content.trim().to_string();
            self.content.clear();
            Some(token)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn nexus_reading_simple_test() {
        let sample = "test_files/simple.nex";
        let mut nex = Nexus::new();
        nex.read(sample).unwrap();
        assert_eq!(1, nex.matrix.len());
    }

    #[test]
    fn nexus_reading_complete_test() {
        let sample = "test_files/complete.nex";
        let mut nex = Nexus::new();
        nex.read(sample).unwrap();
        assert_eq!(5, nex.matrix.len());
    }

    #[test]
    fn nexus_reading_tabulated_test() {
        let sample = "test_files/tabulated.nex";
        let mut nex = Nexus::new();
        nex.read(sample).unwrap();
        assert_eq!(2, nex.matrix.len());
    }

    #[test]
    fn nexus_parsing_object_test() {
        let sample = "test_files/complete.nex";
        let mut nex = Nexus::new();
        nex.read(sample).unwrap();
        assert_eq!(5, nex.ntax);
        assert_eq!(802, nex.nchar);
        assert_eq!("dna", nex.datatype);
        assert_eq!('-', nex.gap);
        assert_eq!('?', nex.missing);
    }

    #[test]
    fn nexus_parse_ntax_test() {
        let tax = "ntax=5";
        let nex = Nexus::new();
        let res = nex.parse_ntax(tax);
        assert_eq!(5, res);
    }

    #[test]
    #[should_panic]
    fn check_invalid_dna_panic_test() {
        let nex = Nexus::new();
        let id = "ABCD";
        let dna = String::from("AGTC?-Z");
        nex.check_valid_dna(id, &dna);
    }

    #[test]
    // #[should_panic]
    fn check_match_ntax_test() {
        let sample = "test_files/simple.nex";
        let mut nex = Nexus::new();
        nex.read(sample).unwrap();
    }

    #[test]
    #[should_panic]
    fn check_match_ntax_panic_test() {
        let sample = "test_files/unmatched_block.nex";
        let mut nex = Nexus::new();
        nex.read(sample).unwrap();
    }

    #[test]
    #[should_panic]
    fn check_invalid_nexus_test() {
        let sample = "test_files/simple.fas";
        let mut nex = Nexus::new();
        nex.read(sample).unwrap();
    }

    #[test]
    #[should_panic]
    fn nexus_duplicate_panic_test() {
        let sample = "test_files/duplicates.nex";
        let mut nex = Nexus::new();
        nex.read(sample).unwrap();
    }

    #[test]
    #[should_panic]
    fn nexus_space_panic_test() {
        let sample = "test_files/idspaces.nex";
        let mut nex = Nexus::new();
        nex.read(sample).unwrap();
    }

    #[test]
    fn nexus_sequence_test() {
        let sample = "test_files/tabulated.nex";
        let mut nex = Nexus::new();
        nex.read(sample).unwrap();
        let key = String::from("ABEF");
        let res = String::from("gatata---");
        assert_eq!(Some(&res), nex.matrix.get(&key));
    }
}
