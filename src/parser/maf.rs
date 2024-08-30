//! Parse Multiple Alignment Format (MAF) files.
//!
//!
//! MAF format is developed by UCSC. It stores multiple sequence alignment in a
//! line-base human readable format. It uses by several aligners, such as LASTZ
//! and Progressive Cactus.
//!
//! # Example
//!
//! ```text
//! ##maf version=1 scoring=roast.v2.0
//! a score=0.000000
//! s scaffold_1 0 10 + 10 ATCGATCGAT
//! s scaffold_2 0 10 + 10 ATCGATCGAT
//!
//! a score=0.000000
//! s scaffold_1 10 10 + 10 ATCGATCGAT
//! s scaffold_2 10 10 + 10 ATCGATCGAT
//! ```
//!
//! The line with `##maf` is the header line. The line with `a` is the alignment
//! line. The line with `s` is the sequence line. It some cases, there are `q`
//! lines, which are the quality lines.
use std::{
    error::Error,
    fmt::Display,
    io::{prelude::*, BufReader},
    str::FromStr,
};

use nom::{bytes::complete, character, combinator, sequence, IResult};

const END_OF_LINE: u8 = b'\n';
const EOF: usize = 0;

/// We can think of a paragraph as a block of text that is separated by a blank
/// line. In MAF format, there are two types of paragraphs: header and alignment.
/// Source: http://genome.ucsc.edu/FAQ/FAQformat.html#format5
pub enum MafParagraph {
    Track(TrackLine),
    Header(MafHeader),
    Comments(String),
    Alignment(MafAlignment),
    Empty,
    Unknown,
}

pub struct MafHeader {
    pub version: String,
    pub scoring: String,
    pub program: Option<String>,
}

impl MafHeader {
    pub fn new() -> Self {
        MafHeader {
            version: String::new(),
            scoring: String::new(),
            program: None,
        }
    }
    pub fn from_str(&mut self, line: &str) -> Result<(), Box<dyn Error>> {
        let mut parts = line.split_whitespace();

        let _ = parts.next(); // Skip the first character

        self.version = self.parse_version(parts.next().unwrap())?;
        self.scoring = self.parse_scoring(parts.next().unwrap())?;
        self.program = parts.next().map(|s| s.to_string());

        Ok(())
    }

    fn parse_version(&self, value: &str) -> Result<String, Box<dyn Error>> {
        let tag: IResult<&str, &str> = sequence::preceded(
            complete::tag("version="),
            character::complete::alphanumeric1,
        )(value);
        match tag {
            Ok((_, value)) => Ok(value.to_string()),
            Err(_) => Err("Error parsing version".into()),
        }
    }

    fn parse_scoring(&self, value: &str) -> Result<String, Box<dyn Error>> {
        let tag: IResult<&str, &str> = sequence::preceded(
            complete::tag("scoring="),
            complete::take_while(|c: char| c.is_alphanumeric() || c == '.'),
        )(value);
        match tag {
            Ok((_, value)) => Ok(value.to_string()),
            Err(_) => Err("Error parsing scoring".into()),
        }
    }
}

pub struct TrackLine {
    pub name: String,
    pub description: Option<String>,
    pub frames: Option<String>,
    pub maf_dot: bool,
    pub visibility: Option<String>,
    pub species_order: Option<String>,
}

impl TrackLine {
    pub fn new() -> Self {
        TrackLine {
            name: String::new(),
            description: None,
            frames: None,
            maf_dot: false,
            visibility: None,
            species_order: None,
        }
    }

    pub fn from_str(&mut self, line: &str) -> Result<(), Box<dyn Error>> {
        let track = self.parse_track(line);

        match track {
            Ok(value) => {
                self.parse_attributes(&value);
                Ok(())
            }
            Err(_) => Err("Error parsing track".into()),
        }
    }

    fn parse_track(&mut self, input: &str) -> Result<String, Box<dyn Error>> {
        let track: IResult<&str, &str> =
            sequence::preceded(complete::tag("track "), complete::take_until("\n"))(input);

        match track {
            Ok((_, value)) => Ok(value.to_string()),
            Err(_) => Err("Error parsing track".into()),
        }
    }

    fn parse_attributes(&mut self, input: &str) {
        let tag: IResult<
            &str,
            (
                &str,
                Option<&str>,
                // Option<&str>,
                // Option<&str>,
                // Option<&str>,
                // Option<&str>,
            ),
        > = sequence::tuple((
            sequence::preceded(complete::tag("name="), character::complete::alphanumeric1),
            combinator::opt(sequence::preceded(
                complete::tag(" description=\""),
                complete::take_until("\""),
            )),
            // combinator::opt(sequence::preceded(
            //     complete::tag(" frames="),
            //     character::complete::alphanumeric1,
            // )),
            // combinator::opt(sequence::preceded(
            //     complete::tag(" mafDot="),
            //     character::complete::alphanumeric1,
            // )),
            // combinator::opt(sequence::preceded(
            //     complete::tag(" visibility="),
            //     character::complete::alphanumeric1,
            // )),
            // combinator::opt(sequence::preceded(
            //     complete::tag(" speciesOrder="),
            //     sequence::delimited(
            //         character::complete::char('"'),
            //         complete::take_until("\""),
            //         character::complete::char('"'),
            //     ),
        ))(input);

        match tag {
            Ok((_, (name, description))) => {
                self.name = name.to_string();
                self.description = description.map(|s| s.to_string());

                // self.frames = frames.map(|s| s.to_string());
                // self.maf_dot = maf_dot.map(|s| s == "on").unwrap_or_default();
                // self.species_order = species_order.map(|s| s.to_string());
                // self.parse_visibility(visibility);
            }
            Err(_) => {}
        }

        self.parse_visibility(input);
    }

    fn parse_visibility(&mut self, value: &str) {
        let tag: IResult<&str, Option<&str>> = combinator::opt(sequence::preceded(
            complete::tag("visibility="),
            character::complete::alphanumeric1,
        ))(value);

        match tag {
            Ok((_, value)) => {
                self.visibility = value.map(|s| s.to_string());
            }
            Err(_) => {}
        }
    }
}

#[derive(Debug, PartialEq)]
pub enum MafVisibility {
    Dense,
    Pack,
    Full,
}

impl Display for MafVisibility {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MafVisibility::Dense => write!(f, "dense"),
            MafVisibility::Pack => write!(f, "pack"),
            MafVisibility::Full => write!(f, "full"),
        }
    }
}

impl FromStr for MafVisibility {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "dense" => Ok(MafVisibility::Dense),
            "pack" => Ok(MafVisibility::Pack),
            "full" => Ok(MafVisibility::Full),
            _ => Err(format!("Unknown visibility: {}", s)),
        }
    }
}

/// The `a` line in MAF format.
/// It contains:
/// - `score`: the score of the alignment. At the the first line
/// - `sequences`: the `s` lines in the alignment.
/// - `information`: the `i` lines in the alignment. It is optional.
/// - `quality`: the `q` lines in the alignment. It is optional.
///
pub struct MafAlignment {
    pub score: f64,
    pub sequences: Vec<MafSequence>,
    pub quality: Option<Quality>,
    pub information: Option<MafInformation>,
    pub empty: Option<MafEmptyLine>,
}

impl MafAlignment {
    pub fn new() -> Self {
        MafAlignment {
            score: 0.0,
            sequences: Vec::new(),
            quality: None,
            information: None,
            empty: None,
        }
    }

    pub fn add_sequence(&mut self, sequence: MafSequence) {
        self.sequences.push(sequence);
    }

    pub fn parse_scores(&self, value: &str) -> Result<f64, Box<dyn Error>> {
        let tag: IResult<&str, &str> = sequence::preceded(
            complete::tag("a score="),
            complete::take_while(|c: char| c.is_ascii_digit() || c == '.'),
        )(value);
        match tag {
            Ok((_, value)) => Ok(value.parse()?),
            Err(_) => Err("Error parsing score".into()),
        }
    }
}

pub struct MafSequence {
    pub source: String,
    pub start: u64,
    pub size: u64,
    pub strand: char,
    pub src_size: u64,
    pub text: Vec<u8>,
}

impl MafSequence {
    pub fn from_str(line: &[u8]) -> Result<Self, Box<dyn Error>> {
        let mut parts = line
            .split(|b| b.is_ascii_whitespace())
            .filter(|b| !b.is_empty());

        let _ = parts.next(); // Skip the first character

        let source = parts.next().unwrap_or_default();
        let start = parts.next().unwrap_or_default();
        let size = parts.next().unwrap_or_default();
        let strand = parts.next().unwrap_or_default();
        let src_size = parts.next().unwrap_or_default();
        let text = parts.next().unwrap_or_default();

        Ok(MafSequence {
            source: String::from_utf8_lossy(source).to_string(),
            start: String::from_utf8_lossy(start).parse().unwrap_or_default(),
            size: String::from_utf8_lossy(size).parse().unwrap_or_default(),
            strand: String::from_utf8_lossy(strand)
                .chars()
                .next()
                .unwrap_or_default(),
            src_size: String::from_utf8_lossy(src_size)
                .parse()
                .unwrap_or_default(),
            text: text.to_vec(),
        })
    }
}

/// The `i` line in MAF format.
/// It describes the information before and after the block of sequences.
/// It contains:
/// - `src`: the source of the information.
/// - `leftStatus`: a character that specifies the relationship between the
///     sequence in the current block and the sequence in the previous block.
/// - `leftCount`: the number of bases in the aligning species between the start and
///     the end of the previous block.
/// - `rightStatus`: a character that specifies the relationship between the
///     sequence in the current block and the sequence in the subsequent block.
/// - `leftValues`: the values of the information before the block of sequences.
/// - `rightCount`: the number of bases in the aligning species between the start of the
///    current block and the start of the next block.
pub struct MafInformation {
    pub source: String,
    pub left_status: char,
    pub left_count: usize,
    pub right_status: char,
    pub right_count: usize,
}

impl MafInformation {
    pub fn from_str(line: &[u8]) -> Result<Self, Box<dyn Error>> {
        let mut parts = line.split(|b| b.is_ascii_whitespace());

        let _ = parts.next(); // Skip the first character

        let source = parts.next().unwrap_or_default();
        let left_status = parts.next().unwrap_or_default();
        let left_count = parts.next().unwrap_or_default();
        let right_status = parts.next().unwrap_or_default();
        let right_count = parts.next().unwrap_or_default();

        Ok(MafInformation {
            source: String::from_utf8_lossy(source).to_string(),
            left_status: String::from_utf8_lossy(left_status)
                .chars()
                .next()
                .unwrap_or_default(),
            left_count: String::from_utf8_lossy(left_count)
                .parse()
                .unwrap_or_default(),
            right_status: String::from_utf8_lossy(right_status)
                .chars()
                .next()
                .unwrap_or_default(),
            right_count: String::from_utf8_lossy(right_count)
                .parse()
                .unwrap_or_default(),
        })
    }
}

/// This is the `e` line in MAF format.
pub struct MafEmptyLine {
    pub source: String,
    pub size: u64,
    pub strand: char,
    pub src_size: u64,
}

impl MafEmptyLine {
    pub fn from_str(line: &[u8]) -> Result<Self, Box<dyn Error>> {
        let mut parts = line.split(|b| b.is_ascii_whitespace());

        let _ = parts.next(); // Skip the first character

        let source = parts.next().unwrap_or_default();
        let size = parts.next().unwrap_or_default();
        let strand = parts.next().unwrap_or_default();
        let src_size = parts.next().unwrap_or_default();

        Ok(MafEmptyLine {
            source: String::from_utf8_lossy(source).to_string(),
            size: String::from_utf8_lossy(size).parse().unwrap_or_default(),
            strand: String::from_utf8_lossy(strand)
                .chars()
                .next()
                .unwrap_or_default(),
            src_size: String::from_utf8_lossy(src_size)
                .parse()
                .unwrap_or_default(),
        })
    }
}

pub struct Quality {
    /// The name of the source sequence.
    pub source: String,
    /// We store the quality values as characters.
    /// It can be 0-9 or F. The F means finished.
    pub values: Vec<char>,
}

impl Quality {
    pub fn from_str(line: &[u8]) -> Result<Self, Box<dyn Error>> {
        let mut parts = line.split(|b| b.is_ascii_whitespace());

        let _ = parts.next(); // Skip the first character

        let source = parts.next().unwrap_or_default();
        let values = parts.next().unwrap_or_default();

        Ok(Quality {
            source: String::from_utf8_lossy(source).to_string(),
            values: String::from_utf8_lossy(values).chars().collect(),
        })
    }
}

pub struct MafReader<R> {
    pub reader: BufReader<R>,
    pub buf: Vec<u8>,
}

impl<R: Read> MafReader<R> {
    pub fn new(file: R) -> Self {
        MafReader {
            reader: BufReader::new(file),
            buf: Vec::new(),
        }
    }

    pub fn next_paragraph(&mut self) -> Option<MafParagraph> {
        let bytes = self
            .reader
            .read_until(END_OF_LINE, &mut self.buf)
            .expect("Error reading file");
        if bytes == EOF {
            return None;
        }

        // Check the first character of the line
        let paragraph = match self.buf[0] {
            b'#' => self.parse_header(),
            b't' => self.parse_track_line(),
            b'a' => self.parse_alignments(),
            b'\n' => Some(MafParagraph::Empty),
            _ => Some(MafParagraph::Unknown),
        };

        self.buf.clear();

        paragraph
    }

    fn parse_header(&mut self) -> Option<MafParagraph> {
        let line = String::from_utf8_lossy(&self.buf);
        if line.starts_with("##maf") {
            let mut header = MafHeader::new();
            header.from_str(&line).unwrap();
            self.buf.clear();
            Some(MafParagraph::Header(header))
        } else {
            self.parse_comments()
        }
    }

    fn parse_comments(&mut self) -> Option<MafParagraph> {
        let line = String::from_utf8_lossy(&self.buf);
        Some(MafParagraph::Comments(line.to_string()))
    }

    fn parse_track_line(&mut self) -> Option<MafParagraph> {
        let line = String::from_utf8_lossy(&self.buf);
        let mut track = TrackLine::new();
        track.from_str(&line).unwrap();
        Some(MafParagraph::Track(track))
    }

    fn parse_alignments(&mut self) -> Option<MafParagraph> {
        let mut alignment = MafAlignment::new();
        alignment
            .parse_scores(&String::from_utf8_lossy(&self.buf))
            .unwrap();
        self.buf.clear();
        loop {
            let bytes = self
                .reader
                .read_until(END_OF_LINE, &mut self.buf)
                .expect("Error reading file");
            if bytes == EOF {
                break;
            }

            match self.buf[0] {
                b's' => {
                    let sequence =
                        MafSequence::from_str(&self.buf).expect("Error parsing sequence");
                    alignment.sequences.push(sequence);
                }
                b'q' => {
                    let quality = Quality::from_str(&self.buf).unwrap();
                    alignment.quality = Some(quality);
                }
                b'i' => {
                    let information = MafInformation::from_str(&self.buf).unwrap();
                    alignment.information = Some(information);
                }
                b'e' => {
                    let empty = MafEmptyLine::from_str(&self.buf).unwrap();
                    alignment.empty = Some(empty);
                }
                b'\n' | b' ' => break,
                _ => {}
            }
            self.buf.clear();
        }

        Some(MafParagraph::Alignment(alignment))
    }
}

impl<R: Read> Iterator for MafReader<R> {
    type Item = MafParagraph;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_paragraph()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_header() {
        let line = "##maf version=1 scoring=roast.v2.0\n";
        let mut header = MafHeader::new();
        header.from_str(line).unwrap();
        assert_eq!(header.version, "1");
        assert_eq!(header.scoring, "roast.v2.0");
    }

    #[test]
    fn test_parse_track_line() {
        let line = "track name=chr1 description=\"Human chromosome\" visibility=pack\n";
        let mut track = TrackLine::new();
        track.from_str(line).unwrap();
        let track_line = track.parse_track(line);
        assert_eq!(
            track_line.unwrap(),
            "name=chr1 description=\"Human chromosome\" visibility=pack"
        );
        // assert_eq!(track.name, "chr1");
        assert_eq!(track.description.unwrap(), "Human chromosome");
        // assert_eq!(track.visibility.unwrap(), MafVisibility::Pack);
    }

    #[test]
    fn test_parse_maf_file() {
        let file = std::fs::File::open("tests/files/maf/simple.maf").unwrap();
        let reader = MafReader::new(file);
        let mut alignments = Vec::new();
        for paragraph in reader {
            match paragraph {
                MafParagraph::Track(track) => {
                    assert_eq!(track.name, "euArc");
                    assert_eq!(track.description, Some(String::from("Primate chromosomes")));
                    // assert_eq!(track.visibility, Some(String::from("pack")));
                }
                MafParagraph::Header(header) => {
                    assert_eq!(header.version, "1");
                    assert_eq!(header.scoring, "tba.v8");
                }
                MafParagraph::Alignment(alignment) => {
                    alignments.push(alignment);
                }
                _ => {}
            }
        }

        assert_eq!(alignments.len(), 3);
    }

    #[test]
    fn parse_visibility() {
        let line = "track name=chr1 description=\"Human chromosome\" visibility=pack\n";
        let mut track = TrackLine::new();
        track.parse_visibility(line);
        assert_eq!(track.visibility, Some(String::from("pack")));
    }
}
