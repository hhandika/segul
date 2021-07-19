use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Lines, Result};
use std::path::Path;

use indexmap::{IndexMap, IndexSet};
use nom::{character::complete, sequence, IResult};

use crate::helper::common::{self, DataType, Header, SeqCheck};

pub struct Phylip<'a> {
    input: &'a Path,
    datatype: &'a DataType,
    pub matrix: IndexMap<String, String>,
    pub header: Header,
}

impl<'a> Phylip<'a> {
    pub fn new(input: &'a Path, datatype: &'a DataType) -> Self {
        Self {
            input,
            datatype,
            matrix: IndexMap::new(),
            header: Header::new(),
        }
    }

    pub fn parse(&mut self) -> Result<()> {
        self.parse_matrix()?;
        let mut seq_info = SeqCheck::new();
        seq_info.get_sequence_info(&self.matrix);
        self.header.aligned = seq_info.is_alignment;
        self.match_header_datatype();
        self.check_ntax_matches();
        self.check_nchar_matches(seq_info.longest);
        Ok(())
    }

    pub fn parse_only_id(&mut self) -> IndexSet<String> {
        let buff = self.get_header().expect("CANNOT READ THE FILE");
        let records = Reader::new(buff, self.header.ntax);
        let mut ids = IndexSet::new();
        records
            .into_iter()
            .take_while(|rec| !rec.interleave)
            .for_each(|rec| {
                if let Some(id) = rec.id {
                    ids.insert(id);
                }
            });

        assert!(
            ids.len() == self.header.ntax,
            "Failed parsing {}. \
        The number of taxa does not match the information in the header.\
            In the header: {} \
            and taxa found: {}. \
            Please check if the header matches the sample count in the phylip file. \
            Or no duplicate IDs are present.
        ",
            self.input.display(),
            self.header.ntax,
            ids.len()
        );

        ids
    }

    fn parse_matrix(&mut self) -> Result<()> {
        let buff = self.get_header()?;
        let records = Reader::new(buff, self.header.ntax);
        let mut ids: IndexMap<usize, String> = IndexMap::new();
        records.into_iter().for_each(|rec| match rec.id {
            Some(id) => {
                common::check_valid_seq(&self.input, &self.datatype, &id, &rec.seq);
                ids.insert(rec.pos, id.clone());
                self.insert_matrix(id, rec.seq);
            }
            None => {
                if let Some(id) = ids.get(&rec.pos) {
                    common::check_valid_seq(&self.input, &self.datatype, &id, &rec.seq);
                    if let Some(value) = self.matrix.get_mut(id) {
                        value.push_str(&rec.seq);
                    }
                }
            }
        });

        Ok(())
    }

    // We return the buffer after reading the header.
    // So we can use the buffer to read the rest of the file.
    fn get_header(&mut self) -> Result<BufReader<File>> {
        let file = File::open(self.input)?;
        let mut buff = BufReader::new(file);
        let mut header_line = String::new();
        buff.read_line(&mut header_line)?;
        self.parse_header(&header_line.trim());

        Ok(buff)
    }

    fn insert_matrix(&mut self, id: String, dna: String) {
        match self.matrix.get(&id) {
            Some(_) => panic!(
                "Found uplicate samples for file {}. First duplicate found: {}",
                self.input.display(),
                id
            ),
            None => {
                self.matrix.insert(id, dna);
            }
        }
    }

    fn parse_header(&mut self, header_line: &str) {
        let header: IResult<&str, (&str, &str)> =
            sequence::separated_pair(complete::digit0, complete::space0, complete::digit0)(
                header_line,
            );

        match header {
            Ok((_, (tax, chars))) => self.parse_num(tax, chars),
            Err(_) => eprintln!("Unknown header! Failed to parse. Ignoring it..."),
        };
    }

    fn parse_num(&mut self, tax: &str, seq: &str) {
        self.header.ntax = tax
            .parse::<usize>()
            .expect("Header taxa number is not a number.");
        self.header.nchar = seq
            .parse::<usize>()
            .expect("Header char length is not a number.");
    }

    fn match_header_datatype(&mut self) {
        if let DataType::Aa = self.datatype {
            self.header.datatype = String::from("protein")
        };
    }

    fn check_ntax_matches(&self) {
        if self.matrix.len() != self.header.ntax {
            panic!(
                "Error reading phylip file: {}. \
            The number of taxa does not match the information in the header.\
            In the header: {} \
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
                "Error reading phylip file {}, \
            The NCHAR value in the header does not match the sequence length. \
            The value in the header {}. \
            The sequence length {}.",
                self.input.display(),
                self.header.nchar,
                longest
            );
        }
    }
}

struct Records {
    id: Option<String>,
    seq: String,
    pos: usize,
    interleave: bool,
}

impl Records {
    fn new() -> Self {
        Self {
            id: None,
            seq: String::new(),
            pos: 0,
            interleave: false,
        }
    }
}

struct Reader<R> {
    reader: Lines<BufReader<R>>,
    interleave: bool,
    pos: usize,
    ntax: usize,
}

impl<R: Read> Reader<R> {
    fn new(file: R, ntax: usize) -> Self {
        Self {
            reader: BufReader::new(file).lines(),
            interleave: false,
            pos: 1,
            ntax,
        }
    }

    fn next_seq(&mut self) -> Option<Records> {
        while let Some(Ok(lines)) = self.reader.by_ref().next() {
            let line = lines.trim();
            if !line.is_empty() {
                let mut records = Records::new();
                if !self.interleave {
                    let seq: Vec<&str> = line.split_whitespace().collect();
                    if seq.len() == 2 {
                        records.id = Some(seq[0].to_string());
                        records.seq = seq[1].to_string();
                    }
                } else {
                    records.id = None;
                    records.seq = line.to_string();
                }

                records.interleave = self.interleave;
                records.pos = self.pos;
                self.pos += 1;

                if self.pos == self.ntax + 1 {
                    self.pos = 1;
                    self.interleave = true;
                }

                return Some(records);
            }
            continue;
        }
        None
    }
}

impl<R: Read> Iterator for Reader<R> {
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
    fn read_phylip_simple_test() {
        let path = Path::new("test_files/simple.phy");
        let mut phylip = Phylip::new(path, &DNA);
        phylip.parse().unwrap();

        assert_eq!(2, phylip.header.ntax);
        assert_eq!(4, phylip.header.nchar);
        assert_eq!(2, phylip.matrix.len());
    }

    #[test]
    #[should_panic]
    fn read_phylip_invalid_test() {
        let path = Path::new("test_files/invalid.phy");
        let mut phylip = Phylip::new(path, &DNA);
        phylip.parse().unwrap();
    }

    #[test]
    fn read_phylip_whitespace_test() {
        let path = Path::new("test_files/whitespaces.phy");
        let mut phylip = Phylip::new(path, &DNA);
        phylip.parse().unwrap();
        assert_eq!(2, phylip.header.ntax);
        assert_eq!(4, phylip.header.nchar);
        assert_eq!(2, phylip.matrix.len());
    }

    #[test]
    fn parse_phylip_header_test() {
        let header = "2 24";
        let mut phy = Phylip::new(Path::new("."), &DNA);
        phy.parse_header(header);

        assert_eq!(2, phy.header.ntax);
        assert_eq!(24, phy.header.nchar);
    }

    #[test]
    fn parse_matrix_interleave_phylip_test() {
        let path = Path::new("test_files/interleave.phy");
        let mut phy = Phylip::new(path, &DNA);
        phy.parse_matrix().unwrap();
        let res = phy.matrix.get("ABCD");
        assert_eq!(Some(&String::from("agccatggaa")), res);
    }

    #[test]
    fn read_interleave_phylip_test() {
        let path = Path::new("test_files/interleave.phy");
        let mut phy = Phylip::new(path, &DNA);
        phy.parse().unwrap();
        let res = phy.matrix.get("ABCD");
        assert_eq!(Some(&String::from("agccatggaa")), res);
    }

    #[test]
    fn reade_int_phylip_whitespaces_test() {
        let path = Path::new("test_files/interleave_whitespaces.phy");
        let mut phy = Phylip::new(path, &DNA);
        phy.parse().unwrap();
        let res = phy.matrix.get("ABCD");
        assert_eq!(Some(&String::from("agccatggcc")), res);
    }

    #[test]
    #[should_panic]
    fn read_interleave_phylip_panic_test() {
        // The header does not matches the character length.
        let path = Path::new("test_files/invalid_interleave.phy");
        let mut phy = Phylip::new(path, &DNA);
        phy.parse().unwrap();
    }

    #[test]
    fn read_interleave_phylip_id_test() {
        let path = Path::new("test_files/interleave.phy");
        let mut phy = Phylip::new(path, &DNA);
        let res = phy.parse_only_id();
        assert_eq!(2, res.len());
    }

    #[test]
    fn phylip_simple_aa_test() {
        let sample = Path::new("test_files/simple_aa.phy");
        let mut phy = Phylip::new(sample, &AA);
        phy.parse().unwrap();
        let key = String::from("ABCE");
        let res = String::from("MAYPMQLGFQDATSPI");
        assert_eq!(Some(&res), phy.matrix.get(&key));
        assert_eq!(String::from("protein"), phy.header.datatype);
    }
}
