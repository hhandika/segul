use std::collections::BTreeMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Lines, Read, Result};
use std::path::Path;

use crate::converter::Converter;

pub fn read_nexus<P: AsRef<Path>>(path: &P) {
    let input = File::open(path).unwrap();
    let buff = BufReader::new(input);
    let mut nex = Nexus::new();
    nex.read(buff).unwrap();

    let matrix = nex.parse_matrix();

    matrix.iter().for_each(|(id, seq)| {
        println!(">{}", id);
        println!("{}", seq);
    });
}

pub fn convert_to_fasta(path: &str) {
    let input = File::open(path).unwrap();
    let buff = BufReader::new(input);
    let mut nex = Nexus::new();

    nex.read(buff).expect("CANNOT READ NEXUS FILES");
    let matrix = nex.parse_matrix();
    let mut convert = Converter::new(path, &matrix);
    convert.write_fasta();
}

pub fn convert_to_phylip(path: &str) {
    let input = File::open(path).unwrap();
    let buff = BufReader::new(input);
    let mut nex = Nexus::new();

    nex.read(buff).expect("CANNOT READ NEXUS FILES");
    let matrix = nex.parse_matrix();
    let mut convert = Converter::new(path, &matrix);
    convert.write_phylip();
}

struct Nexus {
    matrix: String,
}

impl Nexus {
    fn new() -> Self {
        Self {
            matrix: String::new(),
        }
    }

    fn read<R: Read>(&mut self, reader: R) -> Result<()> {
        let reader = Reader::new(reader);
        reader.into_iter().for_each(|r| {
            if r.to_lowercase().contains("matrix") {
                self.matrix = r.trim().to_string();
            }
        });

        Ok(())
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
        let input = File::open(sample).unwrap();
        let buff = BufReader::new(input);
        let mut nex = Nexus::new();
        nex.read(buff).unwrap();
        let read = nex.parse_matrix();
        assert_eq!(1, read.len());
    }

    #[test]
    fn nexus_reading_complete_test() {
        let sample = "test_files/complete.nex";
        let input = File::open(sample).unwrap();
        let buff = BufReader::new(input);
        let mut nex = Nexus::new();
        nex.read(buff).unwrap();
        let read = nex.parse_matrix();
        assert_eq!(5, read.len());
    }

    #[test]
    fn nexus_reading_tabulated_test() {
        let sample = "test_files/tabulated.nex";
        let input = File::open(sample).unwrap();
        let buff = BufReader::new(input);
        let mut nex = Nexus::new();
        nex.read(buff).unwrap();
        let read = nex.parse_matrix();
        assert_eq!(2, read.len());
    }

    #[test]
    #[should_panic]
    fn nexus_duplicate_panic_test() {
        let sample = "test_files/duplicates.nex";
        read_nexus(&sample);
    }

    #[test]
    #[should_panic]
    fn nexus_space_panic_test() {
        let sample = "test_files/idspaces.nex";
        read_nexus(&sample);
    }

    #[test]
    fn nexus_sequence_test() {
        let sample = "test_files/tabulated.nex";
        let input = File::open(sample).unwrap();
        let buff = BufReader::new(input);
        let mut nex = Nexus::new();
        nex.read(buff).unwrap();
        let read = nex.parse_matrix();
        let key = String::from("ABEF");
        let res = String::from("GATATA---");
        assert_eq!(Some(&res), read.get(&key));
    }
}
