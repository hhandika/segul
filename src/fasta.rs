//! A module for parsing fasta files.

use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

use indexmap::IndexMap;
use indexmap::IndexSet;

use crate::common::{Header, SeqCheck};

pub struct Fasta<'a> {
    input: &'a Path,
    pub matrix: IndexMap<String, String>,
    pub is_alignment: bool,
    pub header: Header,
}

pub fn read_only_id(input: &Path) -> IndexSet<String> {
    let file = File::open(input).expect("CANNOT READ A FASTA FILE");
    let buff = BufReader::new(file);
    let mut ids: IndexSet<String> = IndexSet::new();
    buff.lines()
        .filter_map(|ok| ok.ok())
        .filter(|line| line.starts_with('>'))
        .for_each(|line| {
            if let Some(id) = line.strip_prefix('>') {
                ids.insert(id.trim().to_string());
            }
        });
    ids
}

impl<'a> Fasta<'a> {
    pub fn new(input: &'a Path) -> Self {
        Self {
            input,
            matrix: IndexMap::new(),
            is_alignment: false,
            header: Header::new(),
        }
    }

    pub fn read(&mut self) {
        let file = File::open(self.input).expect("CANNOT OPEN THE FILE");
        let buff = BufReader::new(file);
        self.parse_fasta(buff);
        let mut seq_info = SeqCheck::new();
        seq_info.get_sequence_info(&self.matrix);
        self.is_alignment = seq_info.is_alignment;
        self.header.nchar = seq_info.longest;
        self.header.ntax = self.matrix.len();
    }

    fn parse_fasta<R: Read>(&mut self, buff: R) {
        let fasta = FastaReader::new(buff);
        fasta
            .into_iter()
            .for_each(|fas| match self.matrix.get(&fas.id) {
                Some(_) => panic!(
                    "DUPLICATE SAMPLES FOR FILE {}. FIRST DUPLICATE FOUND: {}",
                    self.input.display(),
                    fas.id
                ),
                None => {
                    self.matrix.insert(fas.id, fas.seq);
                }
            });
    }
}

struct Records {
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

struct FastaReader<R> {
    reader: BufReader<R>,
    id: String,
    seq: String,
    found_rec: bool,
}

impl<R: Read> FastaReader<R> {
    fn new(file: R) -> Self {
        Self {
            reader: BufReader::new(file),
            id: String::new(),
            seq: String::new(),
            found_rec: false,
        }
    }

    fn next_read(&mut self) -> Option<Records> {
        while let Some(Ok(line)) = self.reader.by_ref().lines().next() {
            if let Some(id) = line.strip_prefix('>') {
                if self.found_rec {
                    let mut res = Records::new();
                    res.id.push_str(&self.id);
                    res.seq.push_str(&self.seq);
                    self.id = String::from(id);
                    self.seq.clear();
                    return Some(res);
                } else {
                    self.id = String::from(id);
                    self.found_rec = true;
                    self.seq.clear();
                }
                continue;
            }
            self.seq.push_str(line.trim());
        }
        if self.found_rec {
            let mut res = Records::new();
            res.id.push_str(&self.id);
            res.seq.push_str(&self.seq);
            self.id.clear();
            self.found_rec = false;
            self.seq.clear();
            return Some(res);
        } else {
            None
        }
    }
}

impl<R: Read> Iterator for FastaReader<R> {
    type Item = Records;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_read()
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
