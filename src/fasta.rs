use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Lines};
use std::path::Path;

use indexmap::IndexMap;

use crate::common::{Header, SeqCheck};

pub struct Fasta<'a> {
    input: &'a Path,
    pub matrix: IndexMap<String, String>,
    pub is_alignment: bool,
    pub ntax: usize,
    pub nchar: usize,
}

impl SeqCheck for Fasta<'_> {}

impl<'a> Fasta<'a> {
    pub fn new(input: &'a Path) -> Self {
        Self {
            input,
            matrix: IndexMap::new(),
            is_alignment: false,
            ntax: 0,
            nchar: 0,
        }
    }

    pub fn read(&mut self) {
        let file = File::open(self.input).expect("CANNOT OPEN THE FILE");
        let buff = BufReader::new(file);
        self.parse_fasta(buff);
        let (shortest, longest) = self.get_sequence_len(&self.matrix);
        self.is_alignment = self.check_is_alignment(&shortest, &longest);
        self.nchar = longest;
        self.ntax = self.matrix.len();
    }

    pub fn get_header(&self) -> Header {
        Header::new()
    }

    fn parse_fasta<R: Read>(&mut self, buff: R) {
        let fasta = FastaReader::new(buff);
        fasta.into_iter().for_each(|fas| {
            #[allow(clippy::all)]
            if self.matrix.contains_key(&fas.id) {
                panic!("DUPLICATE SAMPLES. FIRST DUPLICATE FOUND: {}", fas.id);
            } else {
                self.matrix.insert(fas.id, fas.seq);
            }
        });
    }
}

pub struct Records {
    id: String,
    seq: String,
}

impl Records {
    fn new() -> Self {
        Self {
            id: String::new(),
            seq: String::new(),
        }
    }
}

pub struct FastaReader<R> {
    reader: Lines<BufReader<R>>,
    pub id: Option<String>,
    pub seq: String,
}

impl<R: Read> FastaReader<R> {
    fn new(file: R) -> Self {
        Self {
            reader: BufReader::new(file).lines(),
            id: None,
            seq: String::new(),
        }
    }
}

impl<R: Read> Iterator for FastaReader<R> {
    type Item = Records;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(Ok(line)) = self.reader.next() {
            if let Some(id) = line.strip_prefix('>') {
                if self.id.is_some() {
                    let mut res = Records::new();
                    res.id.push_str(&self.id.as_ref().unwrap());
                    res.seq.push_str(&self.seq);
                    self.id = Some(String::from(id));
                    self.seq.clear();
                    return Some(res);
                } else {
                    self.id = Some(String::from(id));
                    self.seq.clear();
                }
                continue;
            }
            self.seq.push_str(line.trim());
        }
        if self.id.is_some() {
            let mut res = Records::new();
            res.id.push_str(&self.id.as_ref().unwrap());
            res.seq.push_str(&self.seq);
            self.id = None;
            self.seq.clear();
            self.seq.shrink_to_fit();
            return Some(res);
        }
        None
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn read_fasta_simple_test() {
        let path = Path::new("test_files/simple.fas");
        let mut fasta = Fasta::new(path);
        fasta.read();

        assert_eq!(2, fasta.matrix.len());
    }

    #[test]
    fn check_is_alignment() {
        let path = Path::new("test_files/simple.fas");
        let mut fasta = Fasta::new(path);
        fasta.read();

        assert_eq!(true, fasta.is_alignment);
    }

    #[test]
    fn check_isnot_alignment() {
        let path = Path::new("test_files/unaligned.fas");
        let mut fasta = Fasta::new(path);
        fasta.read();

        assert_eq!(false, fasta.is_alignment);
    }
}
