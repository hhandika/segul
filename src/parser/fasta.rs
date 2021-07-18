//! A module for parsing fasta files.

use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

use indexmap::IndexMap;

use crate::helper::common::{self, Header, SeqCheck};

pub struct Fasta<'a> {
    input: &'a Path,
    pub matrix: IndexMap<String, String>,
    pub is_alignment: bool,
    pub header: Header,
}

pub fn parse_only_id(input: &Path) -> Vec<String> {
    let file = File::open(input).expect("CANNOT READ A FASTA FILE");
    let buff = BufReader::new(file);
    let mut ids: Vec<String> = Vec::new();
    buff.lines()
        .filter_map(|ok| ok.ok())
        .filter(|line| line.starts_with('>'))
        .for_each(|line| {
            if let Some(id) = line.strip_prefix('>') {
                ids.push(id.trim().to_string());
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

    pub fn parse(&mut self) {
        let file = File::open(self.input).expect("CANNOT OPEN THE FILE");
        let buff = BufReader::new(file);
        self.parse_matrix(buff);
        let mut seq_info = SeqCheck::new();
        seq_info.get_sequence_info(&self.matrix);
        self.is_alignment = seq_info.is_alignment;
        self.header.nchar = seq_info.longest;
        self.header.ntax = self.matrix.len();
    }

    fn parse_matrix<R: Read>(&mut self, buff: R) {
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
                    common::check_valid_dna(&self.input, &fas.id, &fas.seq);
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
    fn new(id: &str, seq: &str) -> Self {
        Self {
            id: String::from(id),
            seq: String::from(seq),
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

    fn next_seq(&mut self) -> Option<Records> {
        while let Some(Ok(line)) = self.reader.by_ref().lines().next() {
            if let Some(id) = line.strip_prefix('>') {
                if !self.found_rec {
                    self.id = String::from(id);
                    self.found_rec = true;
                    self.seq.clear();
                } else {
                    let recs = self.get_recs(&self.id, &self.seq);
                    self.id = String::from(id);
                    self.seq.clear();
                    return Some(recs);
                }
            } else {
                self.seq.push_str(line.trim());
            }
        }
        if self.found_rec {
            let recs = self.get_recs(&self.id, &self.seq);
            self.found_rec = false;
            self.id.clear();
            self.seq.clear();
            Some(recs)
        } else {
            None
        }
    }

    fn get_recs(&self, id: &str, seq: &str) -> Records {
        Records::new(id, seq)
    }
}

impl<R: Read> Iterator for FastaReader<R> {
    type Item = Records;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_seq()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn read_fasta_simple_test() {
        let path = Path::new("test_files/simple.fas");
        let mut fasta = Fasta::new(path);
        fasta.parse();

        assert_eq!(2, fasta.matrix.len());
    }

    #[test]
    fn check_is_alignment_test() {
        let path = Path::new("test_files/simple.fas");
        let mut fasta = Fasta::new(path);
        fasta.parse();

        assert_eq!(true, fasta.is_alignment);
    }

    #[test]
    fn check_isnot_alignment_test() {
        let path = Path::new("test_files/unaligned.fas");
        let mut fasta = Fasta::new(path);
        fasta.parse();

        assert_eq!(false, fasta.is_alignment);
    }

    #[test]
    fn interleaved_fas_test() {
        let path = Path::new("test_files/interleave.fas");
        let file = File::open(path).unwrap();
        let rec = FastaReader::new(file);
        let mut seq = IndexMap::new();
        rec.into_iter().for_each(|r| {
            seq.insert(r.id, r.seq);
        });

        let res = String::from("AGTATGATGTATATGTAT");
        let res_2 = String::from("AGTATGATGTATAAAAAA");
        assert_eq!(Some(&res), seq.get("ABCD"));
        assert_eq!(Some(&res_2), seq.get("ABCE"));
    }
}
