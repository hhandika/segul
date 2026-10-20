//! Phylip parser, supports interleaved and sequential formats.
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Result};
use std::path::Path;

use colored::Colorize;
use indexmap::{IndexMap, IndexSet};
use nom::{character::complete, sequence, IResult, Parser};

use crate::helper::alphabet;
use crate::helper::sequence::SeqCheck;
use crate::helper::types::{DataType, Header, SeqMatrix};
use crate::parser;

/// Parse a phylip file sequence IDs.
/// return an instance of `IndexSet<String>` containing sequence IDs
///  and header information.
pub struct Phylip<'a> {
    input: &'a Path,
    datatype: &'a DataType,
    pub matrix: SeqMatrix,
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

    pub fn parse(&mut self) {
        self.parse_matrix().expect("Failed parsing phylip file");
        let mut seq_info = SeqCheck::new();
        seq_info.check(&self.matrix);
        self.header.aligned = seq_info.is_alignment;
        self.match_header_datatype();
        self.check_ntax_matches();
        self.check_nchar_matches(seq_info.longest);
    }

    pub fn parse_only_id(&mut self) -> IndexSet<String> {
        let mut reader = self.get_header().expect("Cannot read phylip file");
        let mut line = String::new();
        let mut ids = IndexSet::new();
        let mut rows = 0;
        // IDs are only in the first block, stop reading after it.
        while rows < self.header.ntax && read_next_line(&mut reader, &mut line) {
            let text = line.trim();
            if !text.is_empty() {
                let (id, _) = split_id_seq(text);
                ids.insert(id.to_string());
                rows += 1;
            }
        }

        parser::warn_duplicate_ids!(self, ids);

        ids
    }

    fn parse_matrix(&mut self) -> Result<()> {
        let mut reader = self.get_header()?;
        let ntax = self.header.ntax;
        self.matrix.reserve(ntax);
        let mut line = String::new();
        // Row position within the current block.
        // Rows keep the same order across blocks,
        // so it matches the matrix index.
        let mut pos = 0;
        let mut interleave = false;
        let mut reserved = false;
        while read_next_line(&mut reader, &mut line) {
            let text = line.trim();
            if text.is_empty() {
                continue;
            }
            if !interleave {
                let (id, seq) = split_id_seq(text);
                alphabet::check_valid_seq(self.input, self.datatype, id, seq);
                parser::insert_matrix!(self, id, seq);
            } else {
                if !reserved {
                    self.reserve_full_length();
                    reserved = true;
                }
                if let Some((id, seq)) = self.matrix.get_index_mut(pos) {
                    alphabet::check_valid_seq(self.input, self.datatype, id, text);
                    seq.push_str(text);
                }
            }

            pos += 1;
            if pos == ntax {
                pos = 0;
                interleave = true;
            }
        }

        Ok(())
    }

    // Reserve the full sequence length once the file is known
    // to be interleaved, avoiding repeated reallocations.
    fn reserve_full_length(&mut self) {
        let nchar = self.header.nchar.min(parser::MAX_RESERVE);
        self.matrix
            .values_mut()
            .for_each(|seq| seq.reserve_exact(nchar.saturating_sub(seq.len())));
    }

    // We return the buffer after reading the header.
    // So we can use the buffer to read the rest of the file.
    fn get_header(&mut self) -> Result<BufReader<File>> {
        let file = File::open(self.input)?;
        let mut buff = BufReader::new(file);
        let mut header_line = String::new();
        buff.read_line(&mut header_line)?;
        self.parse_header(header_line.trim());

        Ok(buff)
    }

    fn parse_header(&mut self, header_line: &str) {
        let header: IResult<&str, (&str, &str)> =
            sequence::separated_pair(complete::digit0, complete::space0, complete::digit0)
                .parse(header_line);

        match header {
            Ok((_, (tax, chars))) => self.parse_num(tax, chars),
            Err(_) => eprintln!("Unknown header! Failed to parse. Ignoring it..."),
        };
    }

    fn parse_num(&mut self, tax: &str, seq: &str) {
        self.header.ntax = tax.parse::<usize>().unwrap_or_else(|_| {
            panic!(
                "Failed parsing taxa counts in file: {}. \
            It is not a number. \
            Make sure it is a relaxed phylip file.",
                self.input.display()
            )
        });
        self.header.nchar = seq.parse::<usize>().unwrap_or_else(|_| {
            panic!(
                "Failed parsing char length in file: {}. \
            It is not a number. \
            Make sure it is a relaxed phylip file.",
                self.input.display()
            )
        });
    }

    #[inline]
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

    #[inline]
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

/// Read the next line into a reused buffer.
/// Returns false at the end of the file.
fn read_next_line<R: BufRead>(reader: &mut R, line: &mut String) -> bool {
    line.clear();
    matches!(reader.read_line(line), Ok(n) if n > 0)
}

/// Split a sequential (or first interleaved block) line
/// into its ID and sequence.
fn split_id_seq(line: &str) -> (&str, &str) {
    parser::split_id_seq(line)
        .unwrap_or_else(|| panic!("An invalid sequence found. Check if IDs contain spaces."))
}

#[cfg(test)]
mod test {
    use super::*;

    const DNA: DataType = DataType::Dna;
    const AA: DataType = DataType::Aa;

    #[test]
    fn test_read_phylip_simple() {
        let path = Path::new("tests/files/simple.phy");
        let mut phylip = Phylip::new(path, &DNA);
        phylip.parse();

        assert_eq!(2, phylip.header.ntax);
        assert_eq!(4, phylip.header.nchar);
        assert_eq!(2, phylip.matrix.len());
    }

    #[test]
    #[should_panic]
    fn test_read_phylip_invalid() {
        let path = Path::new("tests/files/invalid.phy");
        let mut phylip = Phylip::new(path, &DNA);
        phylip.parse();
    }

    #[test]
    fn test_read_phylip_whitespace() {
        let path = Path::new("tests/files/whitespaces.phy");
        let mut phylip = Phylip::new(path, &DNA);
        phylip.parse();
        assert_eq!(2, phylip.header.ntax);
        assert_eq!(4, phylip.header.nchar);
        assert_eq!(2, phylip.matrix.len());
    }

    #[test]
    fn test_parse_phylip_header() {
        let header = "2 24";
        let mut phy = Phylip::new(Path::new("."), &DNA);
        phy.parse_header(header);

        assert_eq!(2, phy.header.ntax);
        assert_eq!(24, phy.header.nchar);
    }

    #[test]
    fn test_parse_matrix_interleave_phylip() {
        let path = Path::new("tests/files/interleave.phy");
        let mut phy = Phylip::new(path, &DNA);
        phy.parse_matrix().unwrap();
        let res = phy.matrix.get("ABCD");
        assert_eq!(Some(&String::from("agccatggaa")), res);
    }

    #[test]
    fn test_read_interleave_phylip() {
        let path = Path::new("tests/files/interleave.phy");
        let mut phy = Phylip::new(path, &DNA);
        phy.parse();
        let res = phy.matrix.get("ABCD");
        assert_eq!(Some(&String::from("agccatggaa")), res);
    }

    #[test]
    fn test_read_int_phylip_whitespaces() {
        let path = Path::new("tests/files/interleave_whitespaces.phy");
        let mut phy = Phylip::new(path, &DNA);
        phy.parse();
        let res = phy.matrix.get("ABCD");
        assert_eq!(Some(&String::from("agccatggcc")), res);
    }

    #[test]
    #[should_panic]
    fn test_read_interleave_phylip_panic() {
        // The header does not matches the character length.
        let path = Path::new("tests/files/invalid_interleave.phy");
        let mut phy = Phylip::new(path, &DNA);
        phy.parse();
    }

    #[test]
    fn test_read_interleave_phylip_id() {
        let path = Path::new("tests/files/interleave.phy");
        let mut phy = Phylip::new(path, &DNA);
        let res = phy.parse_only_id();
        assert_eq!(2, res.len());
    }

    #[test]
    fn test_read_sequential_phylip_id() {
        let path = Path::new("tests/files/simple.phy");
        let mut phy = Phylip::new(path, &DNA);
        let res = phy.parse_only_id();
        assert_eq!(2, res.len());
        assert!(phy.matrix.is_empty());
    }

    #[test]
    fn test_interleave_phylip_capacity() {
        let path = Path::new("tests/files/interleave.phy");
        let mut phy = Phylip::new(path, &DNA);
        phy.parse();
        assert!(phy.matrix.values().all(|seq| seq.capacity() == 10));
    }

    #[test]
    fn test_phylip_simple_aa() {
        let sample = Path::new("tests/files/simple_aa.phy");
        let mut phy = Phylip::new(sample, &AA);
        phy.parse();
        let key = String::from("ABCE");
        let res = String::from("MAYPMQLGFQDATSPI");
        assert_eq!(Some(&res), phy.matrix.get(&key));
        assert_eq!(String::from("protein"), phy.header.datatype);
    }
}
