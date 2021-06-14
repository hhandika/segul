use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, Result};
use std::path::Path;

use indexmap::IndexMap;
use nom::{character::complete, sequence, IResult};

use crate::common::{Header, SeqCheck, SeqFormat, SeqPartition};
use crate::writer::SeqWriter;

pub fn convert_phylip(path: &str, filetype: SeqFormat) {
    let input = Path::new(path);
    let mut phylip = Phylip::new(input);
    phylip.read().expect("CANNOT READ PHYLIP FILES");
    let header = phylip.get_header();
    let mut convert = SeqWriter::new(
        Path::new(path),
        &phylip.matrix,
        header,
        None,
        SeqPartition::None,
    );

    match filetype {
        SeqFormat::Nexus => convert.write_sequence(&filetype),
        SeqFormat::Fasta => convert.write_fasta(),
        _ => (),
    }
}

pub struct Phylip<'a> {
    input: &'a Path,
    pub matrix: IndexMap<String, String>,
    pub ntax: usize,
    pub nchar: usize,
    pub is_alignment: bool,
}

impl SeqCheck for Phylip<'_> {}

impl<'a> Phylip<'a> {
    pub fn new(input: &'a Path) -> Self {
        Self {
            input,
            matrix: IndexMap::new(),
            ntax: 0,
            nchar: 0,
            is_alignment: false,
        }
    }

    pub fn read(&mut self) -> Result<()> {
        let file = File::open(self.input).expect("CANNOT OPEN THE INPUT FILE.");
        let mut buff = BufReader::new(file);

        let mut header_line = String::new();
        buff.read_line(&mut header_line)?;
        self.parse_header(&header_line.trim());

        buff.lines().filter_map(|ok| ok.ok()).for_each(|line| {
            let seq: Vec<&str> = line.split_whitespace().collect();
            match seq.len() {
                2 => self.matrix.insert(seq[0].to_string(), seq[1].to_string()),
                _ => {
                    panic!(
                        "UNSUPPORTED PHYLIP. \
                    THE PROGRAM ONLY WORK WITH NON-INTERLEAVED PHYLIP."
                    )
                }
            };
        });
        self.check_is_alignment(&self.matrix);
        Ok(())
    }

    fn get_header(&self) -> Header {
        let mut header = Header::new();
        header.ntax = Some(self.ntax);
        header.nchar = Some(self.nchar);
        header
    }

    fn parse_header(&mut self, header_line: &str) {
        let header: IResult<&str, (&str, &str)> =
            sequence::separated_pair(complete::digit0, complete::space0, complete::digit0)(
                header_line,
            );

        match header {
            Ok((_, (tax, chars))) => self.parse_num(tax, chars),
            Err(_) => eprintln!("UNKNOWN HEADER! FAILED TO PARSE"),
        };
    }

    fn parse_num(&mut self, tax: &str, seq: &str) {
        self.ntax = tax
            .parse::<usize>()
            .expect("HEADER TAXA NUMBER IS NOT A NUMBER");
        self.nchar = seq
            .parse::<usize>()
            .expect("HEADER CHARS LENGTH IS NOT A NUMBER");
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn read_phylip_simple_test() {
        let path = Path::new("test_files/simple.phy");
        let mut phylip = Phylip::new(path);
        phylip.read().unwrap();

        assert_eq!(2, phylip.ntax);
        assert_eq!(4, phylip.nchar);
        assert_eq!(2, phylip.matrix.len());
    }
}
