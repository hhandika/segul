use std::collections::BTreeMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Lines, Read, Result};

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
        SeqFormat::Phylip => convert.write_phylip(),
        SeqFormat::Fasta => convert.write_fasta(),
    }
}

struct Nexus {
    matrix: BTreeMap<String, String>,
    ntax: usize,
    nchar: usize,
    datatype: String,
    missing: char,
    gap: char,
}

impl Nexus {
    fn new() -> Self {
        Self {
            matrix: BTreeMap::new(),
            ntax: 0,
            nchar: 0,
            datatype: String::new(),
            missing: '?',
            gap: '-',
        }
    }

    fn read(&mut self, path: &str) -> Result<()> {
        let input = File::open(path).expect("CANNOT OPEN THE INPUT FILE");
        let mut buff = BufReader::new(input);
        let mut header = String::new();
        buff.read_line(&mut header)?;
        self.check_nexus(&header.trim());
        let mut matrix = self.parse_nexus(buff);
        self.parse_matrix(&mut matrix);
        Ok(())
    }

    fn check_nexus(&self, line: &str) {
        if !line.to_lowercase().starts_with("#nexus") {
            panic!("INVALID NEXUS FORMAT");
        }
    }

    fn parse_nexus<R: Read>(&self, buff: R) -> String {
        let reader = Reader::new(buff);
        let mut matrix = String::new();
        reader.into_iter().for_each(|read| {
            match read.to_lowercase() {
                command if command.starts_with("matrix") => matrix.push_str(read.trim()),
                _ => (),
            };
        });

        matrix
    }

    fn parse_matrix(&mut self, matrix: &mut String) {
        matrix.pop(); // remove terminated semicolon.
        let content: Vec<&str> = matrix.split('\n').collect();
        content[1..]
            .iter()
            .map(|l| l.trim())
            .filter(|l| !l.is_empty())
            .for_each(|line| {
                let seq: Vec<&str> = line.split_whitespace().collect();
                self.check_seq_len(seq.len());
                let id = seq[0].to_string();
                let dna = seq[1].to_string();
                self.check_valid_dna(&id, &dna);
                #[allow(clippy::all)]
                if self.matrix.contains_key(&id) {
                    panic!("DUPLICATE SAMPLES. FIRST DUPLICATE FOUND: {}", id);
                } else {
                    self.matrix.insert(id, dna);
                }
            });
        matrix.clear();
    }

    fn check_valid_dna(&self, id: &str, dna: &String) {
        if !self.is_valid_dna(dna) {
            panic!("INVALID DNA SEQUENCE FOUND FOR {}", id);
        }
    }

    fn is_valid_dna(&self, dna: &String) -> bool {
        dna.chars().all(|char| common::valid_dna().contains(char))
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
    fn check_valid_dna_test() {
        let nex = Nexus::new();
        let dna = String::from("AGTC?-");
        assert_eq!(true, nex.is_valid_dna(&dna));
    }

    #[test]
    fn check_invalid_dna_test() {
        let nex = Nexus::new();
        let dna = String::from("AGTC?-Z");
        assert_eq!(false, nex.is_valid_dna(&dna));
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
        let res = String::from("GATATA---");
        assert_eq!(Some(&res), nex.matrix.get(&key));
    }
}
