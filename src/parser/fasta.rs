//! A module for parsing fasta files.

use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

use indexmap::{IndexMap, IndexSet};

use crate::helper::alphabet;
use crate::helper::sequence::SeqCheck;
use crate::helper::types::{DataType, Header, SeqMatrix};

pub fn parse_only_id(input: &Path) -> IndexSet<String> {
    let file = File::open(input).expect("Failed opening a fasta file.");
    let buff = BufReader::new(file);
    let mut ids = IndexSet::new();
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

pub struct Fasta<'a> {
    input: &'a Path,
    datatype: &'a DataType,
    pub matrix: SeqMatrix,
    pub header: Header,
}

impl<'a> Fasta<'a> {
    pub fn new(input: &'a Path, datatype: &'a DataType) -> Self {
        Self {
            input,
            datatype,
            matrix: IndexMap::new(),
            header: Header::new(),
        }
    }

    pub fn parse(&mut self) {
        let file = File::open(self.input).expect("Failed opening a fasta file");
        let buff = BufReader::new(file);
        self.parse_matrix(buff);
        let mut seq_info = SeqCheck::new();
        assert!(
            !self.matrix.is_empty(),
            "{} is empty. \
        Make sure the file format is fasta or it is not an empty file!",
            self.input.display()
        );
        seq_info.check(&self.matrix);
        self.header.aligned = seq_info.is_alignment;
        self.header.nchar = seq_info.longest;
        self.header.ntax = self.matrix.len();
        self.match_header_datatype();
    }

    fn parse_matrix<R: Read>(&mut self, buff: R) {
        let fasta = FastaReader::new(buff);
        fasta
            .into_iter()
            .for_each(|fas| match self.matrix.get(&fas.id) {
                Some(original_seq) => panic!(
                    "\nFound duplicate IDs for file {}. \
                    First duplicate ID found: {}. \
                Both sequences are the same: {}.\n
                ",
                    self.input.display(),
                    fas.id,
                    *original_seq == fas.seq
                ),
                None => {
                    alphabet::check_valid_seq(self.input, self.datatype, &fas.id, &fas.seq);
                    self.matrix.insert(fas.id, fas.seq);
                }
            });
    }

    #[inline]
    fn match_header_datatype(&mut self) {
        if let DataType::Aa = self.datatype {
            self.header.datatype = String::from("protein")
        };
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
}

impl<R: Read> FastaReader<R> {
    fn new(file: R) -> Self {
        Self {
            reader: BufReader::new(file),
            id: String::new(),
            seq: String::new(),
        }
    }

    fn next_seq(&mut self) -> Option<Records> {
        while let Some(Ok(line)) = self.reader.by_ref().lines().next() {
            if let Some(id) = line.strip_prefix('>') {
                if self.id.is_empty() {
                    self.id = String::from(id);
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

        // Return the last found record.
        if !self.id.is_empty() {
            let recs = self.get_recs(&self.id, &self.seq);
            self.id.clear();
            self.seq.clear();
            Some(recs)
        } else {
            None
        }
    }

    #[inline]
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

    const DNA: DataType = DataType::Dna;
    const AA: DataType = DataType::Aa;

    #[test]
    fn read_fasta_simple_test() {
        let path = Path::new("tests/files/simple.fas");
        let mut fasta = Fasta::new(path, &DNA);
        fasta.parse();

        assert_eq!(2, fasta.matrix.len());
    }

    #[test]
    fn check_is_alignment_test() {
        let path = Path::new("tests/files/simple.fas");
        let mut fasta = Fasta::new(path, &DNA);
        fasta.parse();

        assert_eq!(true, fasta.header.aligned);
    }

    #[test]
    fn check_isnot_alignment_test() {
        let path = Path::new("tests/files/unaligned.fas");
        let mut fasta = Fasta::new(path, &DNA);
        fasta.parse();

        assert_eq!(false, fasta.header.aligned);
    }

    #[test]
    fn interleaved_fas_test() {
        let path = Path::new("tests/files/interleave.fas");
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

    #[test]
    fn simple_aa_fas_test() {
        let sample = Path::new("tests/files/simple_aa.fas");
        let mut fas = Fasta::new(sample, &AA);
        fas.parse();
        let key = String::from("ABCE");
        let res = String::from("MAYPMQLGFQDATSPI");
        assert_eq!(Some(&res), fas.matrix.get(&key));
        assert_eq!(String::from("protein"), fas.header.datatype);
    }
}
