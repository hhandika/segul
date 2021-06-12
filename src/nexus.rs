use std::collections::BTreeMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Lines, Read, Result};

use crate::converter::Converter;

pub fn convert_to_fasta(path: &str) {
    let mut nex = Nexus::new(path);
    nex.read().expect("CANNOT READ NEXUS FILES");
    let matrix = nex.parse_matrix();
    let mut convert = Converter::new(path, &matrix);
    convert.write_fasta();
}

pub fn convert_to_phylip(path: &str) {
    let mut nex = Nexus::new(path);
    nex.read().expect("CANNOT READ NEXUS FILES");
    let matrix = nex.parse_matrix();
    let mut convert = Converter::new(path, &matrix);
    convert.write_phylip();
}

struct Nexus {
    input: String,
    matrix: String,
}

impl Nexus {
    fn new(path: &str) -> Self {
        Self {
            input: String::from(path),
            matrix: String::new(),
        }
    }

    fn read(&mut self) -> Result<()> {
        let input = File::open(&self.input).expect("CANNOT OPEN THE INPUT FILE");
        let mut buff = BufReader::new(input);
        let mut header = String::new();
        buff.read_line(&mut header).unwrap();
        self.check_nexus(&header.trim());
        let reader = Reader::new(buff);
        reader.into_iter().for_each(|r| {
            if r.to_lowercase().starts_with("matrix") {
                self.matrix = r.trim().to_string();
            }
        });

        Ok(())
    }

    fn check_nexus(&self, line: &str) {
        if !line.to_lowercase().starts_with("#nexus") {
            panic!("INVALID NEXUS FORMAT");
        }
    }

    fn parse_matrix(&mut self) -> BTreeMap<String, String> {
        self.matrix.pop(); // remove terminated semicolon.
        let matrix: Vec<&str> = self.matrix.split('\n').collect();
        let mut seqs = BTreeMap::new();
        matrix[1..]
            .iter()
            .map(|l| l.trim())
            .filter(|l| !l.is_empty())
            .for_each(|l| {
                let seq: Vec<&str> = l.split_whitespace().collect();
                if seq.len() != 2 {
                    panic!(
                        "UNSUPPORTED NEXUS FORMAT. MAKE SURE THERE IS NO SPACE IN THE SAMPLE IDs"
                    );
                }
                let id = seq[0].to_string();
                let dna = seq[1].to_string();
                #[allow(clippy::all)]
                if seqs.contains_key(&id) {
                    panic!("DUPLICATE SAMPLES. FIRST DUPLICATE FOUND: {}", id);
                } else {
                    seqs.insert(id, dna);
                }
            });
        seqs
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
        let mut nex = Nexus::new(sample);
        nex.read().unwrap();
        let read = nex.parse_matrix();
        assert_eq!(1, read.len());
    }

    #[test]
    fn nexus_reading_complete_test() {
        let sample = "test_files/complete.nex";
        let mut nex = Nexus::new(sample);
        nex.read().unwrap();
        let read = nex.parse_matrix();
        assert_eq!(5, read.len());
    }

    #[test]
    fn nexus_reading_tabulated_test() {
        let sample = "test_files/tabulated.nex";
        let mut nex = Nexus::new(sample);
        nex.read().unwrap();
        let read = nex.parse_matrix();
        assert_eq!(2, read.len());
    }

    #[test]
    #[should_panic]
    fn check_invalid_nexus_test() {
        let sample = "test_files/simple.fas";
        let mut nex = Nexus::new(sample);
        nex.read().unwrap();
    }

    #[test]
    #[should_panic]
    fn nexus_duplicate_panic_test() {
        let sample = "test_files/duplicates.nex";
        let mut nex = Nexus::new(sample);
        nex.read().unwrap();
        nex.parse_matrix();
    }

    #[test]
    #[should_panic]
    fn nexus_space_panic_test() {
        let sample = "test_files/idspaces.nex";
        let mut nex = Nexus::new(sample);
        nex.read().unwrap();
        nex.parse_matrix();
    }

    #[test]
    fn nexus_sequence_test() {
        let sample = "test_files/tabulated.nex";
        let mut nex = Nexus::new(sample);
        nex.read().unwrap();
        let read = nex.parse_matrix();
        let key = String::from("ABEF");
        let res = String::from("GATATA---");
        assert_eq!(Some(&res), read.get(&key));
    }
}
