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

use nom::{bytes::complete, character, sequence, IResult};

const END_OF_LINE: u8 = b'\n';
const EOF: usize = 0;

/// We can think of a paragraph as a block of text that is separated by a blank
/// line. In MAF format, there are two types of paragraphs: header and alignment.
/// Source: http://genome.ucsc.edu/FAQ/FAQformat.html#format5
pub enum MafParagraf {
    Track(TrackLine),
    Header(MafHeader),
    Comments(String),
    Alignment(MafAlignment),
    Empty,
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
    pub maf_dot: Option<String>,
    pub visiblity: Option<String>,
    pub species_order: Option<String>,
}

impl TrackLine {
    pub fn new() -> Self {
        TrackLine {
            name: String::new(),
            description: None,
            frames: None,
            maf_dot: None,
            visiblity: None,
            species_order: None,
        }
    }

    pub fn from_str(&mut self, line: &str) -> Result<(), Box<dyn Error>> {
        let mut parts = line.split_whitespace();

        let _ = parts.next(); // Skip the first character

        self.name = parts.next().unwrap_or_default().to_string();
        self.description = parts.next().map(|s| s.to_string());
        self.frames = parts.next();
        self.maf_dot = parts.next();
        self.visiblity = parts.next();
        self.species_order = parts.next();

        Ok(())
    }
}

pub enum MafVisiblity {
    Dense,
    Pack,
    Full,
}

impl Display for MafVisiblity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MafVisiblity::Dense => write!(f, "dense"),
            MafVisiblity::Pack => write!(f, "pack"),
            MafVisiblity::Full => write!(f, "full"),
        }
    }
}

impl FromStr for MafVisiblity {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "dense" => Ok(MafVisiblity::Dense),
            "pack" => Ok(MafVisiblity::Pack),
            "full" => Ok(MafVisiblity::Full),
            _ => Err(format!("Unknown visiblity: {}", s)),
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
        let mut parts = line.split(|b| b.is_ascii_whitespace());

        let _ = parts.next(); // Skip the first character

        let source = parts.next().unwrap_or_default();
        let start = parts.next().unwrap_or_default();
        let size = parts.next().unwrap_or_default();
        let strand = parts.next().unwrap_or_default();
        let src_size = parts.next().unwrap_or_default();

        let text = parts.next().unwrap();

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
/// - `righCount`: the number of bases in the aligning species between the start of the
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

    pub fn next_paragraph(&mut self) -> Option<MafParagraf> {
        let bytes = self
            .reader
            .read_until(END_OF_LINE, &mut self.buf)
            .expect("Error reading file");
        if bytes == EOF {
            return None;
        }

        // Check the first character of the line
        match self.buf[0] {
            b'#' => {
                if self.buf[1] == b'#' {
                    self.parse_header()
                } else {
                    self.parse_track_line()
                }
            }
            b'a' => self.parse_alignments(),
            b'\n' => Some(MafParagraf::Empty),
            _ => None,
        }
    }

    fn parse_header(&mut self) -> Option<MafParagraf> {
        let line = String::from_utf8_lossy(&self.buf);
        let mut header = MafHeader::new();
        header.from_str(&line).unwrap();
        Some(MafParagraf::Header(header))
    }

    fn parse_track_line(&mut self) -> Option<MafParagraf> {
        let line = String::from_utf8_lossy(&self.buf);
        let mut track = TrackLine::new();
        track.from_str(&line).unwrap();
        Some(MafParagraf::Track(track))
    }

    fn parse_alignments(&mut self) -> Option<MafParagraf> {
        let mut alignment = MafAlignment::new();

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
                b'\n' => {
                    break;
                }
                _ => {}
            }
        }

        Some(MafParagraf::Alignment(alignment))
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
        let line = "track name=chr1 description=\"Human chromosome
        # 1\" visibility=pack\n";
        let mut track = TrackLine::new();
        track.from_str(line).unwrap();
        assert_eq!(track.name, "chr1");
        assert_eq!(track.description.unwrap(), "Human chromosome 1");
        assert_eq!(track.visiblity.unwrap(), "pack");
    }
}
