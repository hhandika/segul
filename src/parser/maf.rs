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

use nom::{
    bytes::complete,
    character,
    combinator::{self},
    sequence, IResult,
};

#[cfg(target_os = "windows")]
const CAR_RETURN: u8 = b'\r';

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
    pub scoring: Option<String>,
    pub program: Option<String>,
}

impl MafHeader {
    pub fn new() -> Self {
        MafHeader {
            version: String::new(),
            scoring: None,
            program: None,
        }
    }

    pub fn from_str(&mut self, line: &str) -> Result<(), Box<dyn Error>> {
        let mut parts = line.split_whitespace();

        let _ = parts.next(); // Skip the first character

        self.version = self.parse_version(parts.next().unwrap())?;
        self.scoring = parts.next().map(|s| self.parse_scoring(s)).transpose()?;
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
        let (_, tract) = self.parse_track(line).unwrap();
        let (_, name) = self.parse_name(tract).unwrap();

        self.name = name.to_string();

        if tract.contains("description") {
            let input = self.capture_tags(tract, "description=");
            if let Some(input) = input {
                let (_, description) = self.parse_description(&input).unwrap_or_default();
                self.description = description.map(|s| s.to_string());
            }
        }

        if tract.contains("frames") {
            let input = self.capture_tags(tract, "frames=");
            if let Some(input) = input {
                let (_, frames) = self.parse_frames(input.as_str()).unwrap_or_default();
                self.frames = frames.map(|s| s.to_string());
            }
        }

        if tract.contains("mafDot") {
            let input = self.capture_tags(tract, "mafDot=");
            if let Some(input) = input {
                let (_, maf_dot) = self.parse_maf_dot(input.as_str()).unwrap_or_default();
                self.maf_dot = maf_dot;
            }
        }

        if tract.contains("visibility") {
            let input = self.capture_tags(tract, "visibility=");
            if let Some(input) = input {
                let (_, visibility) = self.parse_visibility(input.as_str()).unwrap_or_default();
                self.visibility = visibility.map(|s| s.to_string());
            }
        }

        if tract.contains("speciesOrder") {
            let input = self.capture_tags(tract, "speciesOrder=");
            if let Some(input) = input {
                let (_, species_order) =
                    self.parse_species_order(input.as_str()).unwrap_or_default();
                self.species_order = species_order.map(|s| s.to_string());
            }
        }

        Ok(())
    }

    fn parse_track<'a>(&mut self, input: &'a str) -> IResult<&'a str, &'a str> {
        let track: IResult<&str, &str> =
            sequence::preceded(complete::tag("track "), complete::take_until("\n"))(input);
        track
    }

    fn parse_name<'a>(&mut self, input: &'a str) -> IResult<&'a str, &'a str> {
        let name: IResult<&str, &str> =
            sequence::preceded(complete::tag("name="), character::complete::alphanumeric1)(input);
        name
    }

    fn parse_description<'a>(&mut self, input: &'a str) -> IResult<&'a str, Option<&'a str>> {
        let description: IResult<&str, Option<&str>> = combinator::opt(sequence::preceded(
            complete::tag("description="),
            sequence::delimited(
                character::complete::char('"'),
                complete::take_while(|c: char| c != '"'),
                character::complete::char('"'),
            ),
        ))(input);

        description
    }

    fn parse_frames<'a>(&mut self, input: &'a str) -> IResult<&'a str, Option<&'a str>> {
        let (input, _) = character::complete::space0(input)?;
        let frames: IResult<&str, Option<&str>> = combinator::opt(sequence::preceded(
            complete::tag("frames="),
            character::complete::alphanumeric1,
        ))(input);
        frames
    }

    fn parse_maf_dot<'a>(&mut self, input: &'a str) -> IResult<&'a str, bool> {
        let (input, _) = character::complete::space0(input)?;
        let maf_dot: IResult<&str, &str> =
            sequence::preceded(complete::tag("mafDot="), character::complete::alphanumeric1)(input);
        match maf_dot {
            Ok((_, value)) => Ok((input, value == "on")),
            Err(_) => Ok((input, false)),
        }
    }

    fn parse_visibility<'a>(&mut self, input: &'a str) -> IResult<&'a str, Option<&'a str>> {
        let (input, _) = character::complete::space0(input)?;
        combinator::opt(sequence::preceded(
            complete::tag("visibility="),
            character::complete::alphanumeric1,
        ))(input)
    }

    fn parse_species_order<'a>(&mut self, input: &'a str) -> IResult<&'a str, Option<&'a str>> {
        let (input, _) = character::complete::space0(input)?;
        let species_order: IResult<&str, Option<&str>> = combinator::opt(sequence::preceded(
            complete::tag("speciesOrder="),
            sequence::delimited(
                character::complete::char('"'),
                complete::take_until("\""),
                character::complete::char('"'),
            ),
        ))(input);
        species_order
    }

    // Split the line at tags and the second part is the value.
    fn capture_tags(&self, input: &str, tag: &str) -> Option<String> {
        // use find pos to get the position of the first space
        let pos = input.find(tag);
        match pos {
            Some(pos) => Some(input[pos..].to_string()),
            None => None,
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
    pub score: Option<f64>,
    pub sequences: Vec<MafSequence>,
    pub quality: Option<Quality>,
    pub information: Option<MafInformation>,
    pub empty: Option<MafEmptyLine>,
}

impl MafAlignment {
    pub fn new() -> Self {
        MafAlignment {
            score: None,
            sequences: Vec::new(),
            quality: None,
            information: None,
            empty: None,
        }
    }

    pub fn add_sequence(&mut self, sequence: MafSequence) {
        self.sequences.push(sequence);
    }

    pub fn parse_scores(&self, value: &str) -> Option<f64> {
        let tag: IResult<&str, &str> = sequence::preceded(
            complete::tag("a score="),
            complete::take_while(|c: char| c.is_ascii_digit() || c == '.'),
        )(value);
        match tag {
            Ok((_, value)) => Some(
                value
                    .parse()
                    .expect("Error parsing score. It must be a number"),
            ),
            Err(_) => None,
        }
    }
}

pub struct MafSequence {
    pub source: String,
    pub start: usize,
    pub size: usize,
    pub strand: char,
    pub src_size: usize,
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

#[derive(Debug)]
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

    /// Read the next paragraph from the file.
    /// In MAF terms, a paragraph is a block of text that is separated by a blank line.
    /// It can be a header, a track line,or an alignment.
    pub fn next_paragraph(&mut self) -> Option<MafParagraph> {
        let bytes = self
            .reader
            .read_until(END_OF_LINE, &mut self.buf)
            .expect("Error reading file");
        if bytes == EOF {
            return None;
        }

        #[cfg(target_os = "windows")]
        self.check_carriage_return();

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

    // Check carriage return for windows
    #[cfg(target_os = "windows")]
    fn check_carriage_return(&mut self) {
        if self.buf.contains(&CAR_RETURN) {
            // Remove the carriage return
            // and leave only the line feed
            self.buf.retain(|&c| c != CAR_RETURN);
        }
    }

    fn parse_header(&mut self) -> Option<MafParagraph> {
        // We convert to a string since this line is short
        // and we can parse it easily.
        let line = String::from_utf8_lossy(&self.buf);
        if line.starts_with("##maf") {
            let mut header = MafHeader::new();
            header.from_str(&line).unwrap();
            self.buf.clear();
            Some(MafParagraph::Header(header))
        } else {
            // It must have been a comment if it does not start with ##maf
            // Typically it starts with a single #
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
        alignment.parse_scores(&String::from_utf8_lossy(&self.buf));
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
                END_OF_LINE | b' ' => break,
                #[cfg(target_os = "windows")]
                CAR_RETURN => break,
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
        assert_eq!(header.scoring, Some(String::from("roast.v2.0")));
    }

    #[test]
    fn test_parse_track_line() {
        let line = "track name=chr1 description=\"Human chromosome\" visibility=pack\n";
        let mut track = TrackLine::new();
        track.from_str(line).unwrap();
        track.from_str(line).unwrap();
        assert_eq!(track.name, "chr1");
        assert_eq!(track.description, Some(String::from("Human chromosome")));
        assert_eq!(track.visibility, Some(String::from("pack")));
    }

    #[test]
    fn parse_name() {
        let line = "name=chr1";
        let mut track = TrackLine::new();
        let name = track.parse_name(line).unwrap();
        let (_, res) = name;
        assert_eq!(res, "chr1");
    }

    #[test]
    fn parse_description() {
        let line = " description=\"Human chromosome\"";
        let line2 = " description=\"Human chromosome\" visibility=pack";
        let mut track = TrackLine::new();
        let tags = track.capture_tags(line, "description=").unwrap();
        let description = track.parse_description(&tags).unwrap();
        let (_, res) = description;
        assert_eq!(res, Some("Human chromosome"));
        let tags2 = track.capture_tags(line2, "description=").unwrap();
        let description2 = track.parse_description(&tags2).unwrap();
        let (_, res2) = description2;
        assert_eq!(res2, Some("Human chromosome"));
    }

    #[test]
    fn parse_frames() {
        let line = " frames=multiz28wayFrames";
        let mut track = TrackLine::new();
        let frames = track.parse_frames(line).unwrap();
        let (_, res) = frames;
        assert_eq!(res, Some("multiz28wayFrames"));
    }

    #[test]
    fn parse_maf_dot() {
        let line = " mafDot=on";
        let mut track = TrackLine::new();
        let tags = track.capture_tags(line, "mafDot=").unwrap();
        let maf_dot = track.parse_maf_dot(&tags).unwrap();
        let (_, res) = maf_dot;
        assert_eq!(res, true);
    }

    #[test]
    fn parse_visibility() {
        let line = " visibility=pack";
        let mut track = TrackLine::new();
        let tags = track.capture_tags(line, "visibility=").unwrap();
        let visibility = track.parse_visibility(&tags).unwrap();
        let (_, res) = visibility;
        assert_eq!(res, Some("pack"));

        let complete_line = " description=\"Human chromosome\" visibility=pack\n";
        let tags2 = track.capture_tags(complete_line, "visibility=").unwrap();
        let visibility2 = track.parse_visibility(&tags2).unwrap();
        let (_, res2) = visibility2;
        assert_eq!(res2, Some("pack"));
    }

    #[test]
    fn parse_species_order() {
        let line = " speciesOrder=\"hg18 panTro2\"";
        let mut track = TrackLine::new();
        let species_order = track.parse_species_order(line).unwrap();
        let (_, res) = species_order;
        assert_eq!(res, Some("hg18 panTro2"));
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
                    assert_eq!(track.visibility, Some(String::from("pack")));
                }
                MafParagraph::Header(header) => {
                    assert_eq!(header.version, "1");
                    assert_eq!(header.scoring, Some(String::from("tba.v8")));
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
    fn test_capture_tags() {
        let line = "track name=chr1 description=\"Human chromosome\" visibility=pack";
        let track = TrackLine::new();
        let tags = track.capture_tags(line, "description=");
        assert_eq!(
            tags,
            Some("description=\"Human chromosome\" visibility=pack".to_string())
        );
    }
}
