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
    io::{BufReader, prelude::*},
    str::FromStr,
};

use anyhow::{Context, anyhow, bail};
use nom::{
    IResult, Parser,
    bytes::complete,
    character,
    combinator::{self},
    sequence,
};

use crate::helper::types::DnaStrand;

use super::{END_OF_LINE, EOF};

/// We can think of a paragraph as a block of text that is separated by a blank
/// line. In MAF format, there are two types of paragraphs: header and alignment.
/// Source: <http://genome.ucsc.edu/FAQ/FAQformat.html#format5>
pub enum MafParagraph {
    Track(TrackLine),
    Header(MafHeader),
    Comments(String),
    Alignment(MafAlignment),
    Empty,
    Unknown,
}

#[derive(Debug, Default)]
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

        let version = parts.next().ok_or("Missing MAF version")?;
        self.version = self.parse_version(version)?;
        self.scoring = parts.next().map(|s| self.parse_scoring(s)).transpose()?;
        self.program = parts.next().map(|s| s.to_string());

        Ok(())
    }

    fn parse_version(&self, value: &str) -> Result<String, Box<dyn Error>> {
        let tag: IResult<&str, &str> = sequence::preceded(
            complete::tag("version="),
            character::complete::alphanumeric1,
        ).parse(value);
        match tag {
            Ok((_, value)) => Ok(value.to_string()),
            Err(_) => Err("Error parsing version".into()),
        }
    }

    fn parse_scoring(&self, value: &str) -> Result<String, Box<dyn Error>> {
        let tag: IResult<&str, &str> = sequence::preceded(
            complete::tag("scoring="),
            complete::take_while(|c: char| c.is_alphanumeric() || c == '.'),
        ).parse(value);
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

impl Default for TrackLine {
    fn default() -> Self {
        Self::new()
    }
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
        let line = line.trim_end();
        let (_, tract) = self
            .parse_track(line)
            .map_err(|_| "Track line must start with 'track '")?;
        let (_, name) = self
            .parse_name(tract)
            .map_err(|_| "Track line must have a name")?;

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
            sequence::preceded(complete::tag("track "), combinator::rest).parse(input);
        track
    }

    fn parse_name<'a>(&mut self, input: &'a str) -> IResult<&'a str, &'a str> {
        let name: IResult<&str, &str> =
            sequence::preceded(complete::tag("name="), character::complete::alphanumeric1).parse(input);
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
        )).parse(input);

        description
    }

    fn parse_frames<'a>(&mut self, input: &'a str) -> IResult<&'a str, Option<&'a str>> {
        let (input, _) = character::complete::space0(input)?;
        let frames: IResult<&str, Option<&str>> = combinator::opt(sequence::preceded(
            complete::tag("frames="),
            character::complete::alphanumeric1,
        )).parse(input);
        frames
    }

    fn parse_maf_dot<'a>(&mut self, input: &'a str) -> IResult<&'a str, bool> {
        let (input, _) = character::complete::space0(input)?;
        let maf_dot: IResult<&str, &str> =
            sequence::preceded(complete::tag("mafDot="), character::complete::alphanumeric1).parse(input);
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
        )).parse(input)
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
        )).parse(input);
        species_order
    }

    // Split the line at tags and the second part is the value.
    fn capture_tags(&self, input: &str, tag: &str) -> Option<String> {
        // use find pos to get the position of the first space
        let pos = input.find(tag);
        pos.map(|pos| input[pos..].to_string())
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
#[derive(Debug, Clone)]
pub struct MafAlignment {
    pub score: Option<f64>,
    pub sequences: Vec<MafSequence>,
    pub quality: Option<Quality>,
    pub information: Option<MafInformation>,
    pub empty: Option<MafEmptyLine>,
}

impl Default for MafAlignment {
    fn default() -> Self {
        Self::new()
    }
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

    /// Parse the `key=value` pairs on the `a` line and keep the score.
    /// Negative and exponential scores are supported.
    /// The score is `None` when missing or invalid.
    pub fn parse_scores(&mut self, value: &str) {
        self.score = value
            .split_whitespace()
            .filter_map(|field| field.split_once('='))
            .find(|(key, _)| *key == "score")
            .and_then(|(_, score)| score.parse::<f64>().ok());
    }
}

#[derive(Debug, Clone)]
pub struct MafSequence {
    pub source: String,
    pub start: usize,
    pub size: usize,
    pub strand: DnaStrand,
    pub src_size: usize,
    pub text: Vec<u8>,
}

impl MafSequence {
    /// Parse an `s` line: `s src start size strand srcSize text`.
    /// Every field is required and the strand must be `+` or `-`.
    pub fn from_buf(line: &[u8]) -> anyhow::Result<Self> {
        let fields: Vec<&[u8]> = line
            .split(|b| b.is_ascii_whitespace())
            .filter(|b| !b.is_empty())
            .collect();
        if fields.len() != 7 {
            bail!(
                "expected 7 fields in a sequence line (s src start size strand srcSize text), \
                found {}",
                fields.len()
            );
        }
        let strand = match fields[4] {
            b"+" => DnaStrand::Forward,
            b"-" => DnaStrand::Reverse,
            other => bail!("invalid strand '{}'", String::from_utf8_lossy(other)),
        };
        Ok(MafSequence {
            source: String::from_utf8_lossy(fields[1]).to_string(),
            start: Self::parse_number(fields[2], "start")?,
            size: Self::parse_number(fields[3], "size")?,
            strand,
            src_size: Self::parse_number(fields[5], "srcSize")?,
            text: fields[6].to_vec(),
        })
    }

    fn parse_number(field: &[u8], name: &str) -> anyhow::Result<usize> {
        let text = String::from_utf8_lossy(field);
        text.parse()
            .with_context(|| format!("invalid {} '{}'", name, text))
    }
}

/// The `i` line in MAF format.
/// It describes the information before and after the block of sequences.
/// It contains:
/// - `src`: the source of the information.
/// - `leftStatus`: a character that specifies the relationship between the
///   sequence in the current block and the sequence in the previous block.
/// - `leftCount`: the number of bases in the aligning species between the start and
///   the end of the previous block.
/// - `rightStatus`: a character that specifies the relationship between the
///   sequence in the current block and the sequence in the subsequent block.
/// - `leftValues`: the values of the information before the block of sequences.
/// - `rightCount`: the number of bases in the aligning species between the start of the
///   current block and the start of the next block.
#[derive(Debug, Clone)]
pub struct MafInformation {
    pub source: String,
    pub left_status: char,
    pub left_count: usize,
    pub right_status: char,
    pub right_count: usize,
}

impl MafInformation {
    pub fn from_buf(line: &[u8]) -> Result<Self, Box<dyn Error>> {
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
#[derive(Debug, Clone)]
pub struct MafEmptyLine {
    pub source: String,
    pub size: u64,
    pub strand: char,
    pub src_size: u64,
}

impl MafEmptyLine {
    pub fn from_buf(line: &[u8]) -> Result<Self, Box<dyn Error>> {
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

#[derive(Debug, Clone)]
pub struct Quality {
    /// The name of the source sequence.
    pub source: String,
    /// We store the quality values as characters.
    /// It can be 0-9 or F. The F means finished.
    pub values: Vec<char>,
}

impl Quality {
    pub fn from_buf(line: &[u8]) -> Result<Self, Box<dyn Error>> {
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
    line_number: usize,
    // An `a` line read while finishing the previous block.
    pending: Option<Vec<u8>>,
}

impl<R: Read> MafReader<R> {
    pub fn new(file: R) -> Self {
        MafReader {
            reader: BufReader::new(file),
            line_number: 0,
            pending: None,
        }
    }

    /// Read the next paragraph from the file.
    /// In MAF terms, a paragraph is a block of text that is separated by a blank line.
    /// It can be a header, a track line, or an alignment.
    /// Returns `None` at the end of the file.
    pub fn next_paragraph(&mut self) -> Option<anyhow::Result<MafParagraph>> {
        let line = match self.read_line() {
            Ok(Some(line)) => line,
            Ok(None) => return None,
            Err(e) => return Some(Err(e)),
        };
        if line.trim_ascii().is_empty() {
            return Some(Ok(MafParagraph::Empty));
        }
        let paragraph = match line[0] {
            b'#' => Ok(self.parse_header(&line)),
            b't' => Ok(self.parse_track_line(&line)),
            b'a' => self.parse_alignments(&line),
            _ => Ok(MafParagraph::Unknown),
        };
        Some(paragraph)
    }

    /// Read one line without its line ending (`\n` or `\r\n`).
    fn read_line(&mut self) -> anyhow::Result<Option<Vec<u8>>> {
        if let Some(line) = self.pending.take() {
            return Ok(Some(line));
        }
        let mut buf = Vec::new();
        let bytes = self
            .reader
            .read_until(END_OF_LINE, &mut buf)
            .with_context(|| format!("Failed reading MAF line {}", self.line_number + 1))?;
        if bytes == EOF {
            return Ok(None);
        }
        self.line_number += 1;
        while matches!(buf.last(), Some(b'\n' | b'\r')) {
            buf.pop();
        }
        Ok(Some(buf))
    }

    fn parse_header(&self, line: &[u8]) -> MafParagraph {
        // We convert to a string since this line is short
        // and we can parse it easily.
        let line = String::from_utf8_lossy(line);
        if line.starts_with("##maf") {
            let mut header = MafHeader::new();
            match header.from_str(&line) {
                Ok(()) => MafParagraph::Header(header),
                Err(e) => {
                    log::warn!("Skipping MAF header at line {}: {}", self.line_number, e);
                    MafParagraph::Unknown
                }
            }
        } else {
            // It must have been a comment if it does not start with ##maf
            // Typically it starts with a single #
            MafParagraph::Comments(line.to_string())
        }
    }

    fn parse_track_line(&self, line: &[u8]) -> MafParagraph {
        let line = String::from_utf8_lossy(line);
        let mut track = TrackLine::new();
        match track.from_str(&line) {
            Ok(()) => MafParagraph::Track(track),
            Err(e) => {
                log::warn!("Skipping MAF track line {}: {}", self.line_number, e);
                MafParagraph::Unknown
            }
        }
    }

    fn parse_alignments(&mut self, first_line: &[u8]) -> anyhow::Result<MafParagraph> {
        let mut alignment = MafAlignment::new();
        alignment.parse_scores(&String::from_utf8_lossy(first_line));
        while let Some(line) = self.read_line()? {
            if line.trim_ascii().is_empty() {
                break;
            }
            match line[0] {
                b's' => {
                    let sequence = MafSequence::from_buf(&line).with_context(|| {
                        format!("Invalid MAF sequence line {}", self.line_number)
                    })?;
                    alignment.sequences.push(sequence);
                }
                b'q' => alignment.quality = Quality::from_buf(&line).ok(),
                b'i' => alignment.information = MafInformation::from_buf(&line).ok(),
                b'e' => alignment.empty = MafEmptyLine::from_buf(&line).ok(),
                // A new block without a separating blank line.
                b'a' => {
                    self.pending = Some(line);
                    break;
                }
                _ => {}
            }
        }
        if alignment.sequences.is_empty() {
            return Err(anyhow!(
                "MAF alignment block ending at line {} has no sequence lines",
                self.line_number
            ));
        }
        Ok(MafParagraph::Alignment(alignment))
    }
}

impl<R: Read> Iterator for MafReader<R> {
    type Item = anyhow::Result<MafParagraph>;

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
    fn test_score() {
        let line = "a score=0.000000\n";
        let mut alignment = MafAlignment::new();
        alignment.parse_scores(line);
        assert_eq!(alignment.score, Some(0.0));
    }

    #[test]
    fn test_score_int() {
        let line = "a score=1234\n";
        let mut alignment = MafAlignment::new();
        alignment.parse_scores(line);
        assert_eq!(alignment.score, Some(1234.0));
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
        assert!(res);
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
            match paragraph.unwrap() {
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
        assert_eq!(alignments[0].score, Some(23262.0));
    }

    fn read_blocks(text: &str) -> anyhow::Result<Vec<MafAlignment>> {
        let mut blocks = Vec::new();
        for paragraph in MafReader::new(text.as_bytes()) {
            if let MafParagraph::Alignment(aln) = paragraph? {
                blocks.push(aln);
            }
        }
        Ok(blocks)
    }

    #[test]
    fn test_crlf_blocks_are_not_fused() {
        let text = std::fs::read_to_string("tests/files/maf/simple.maf")
            .unwrap()
            .replace('\n', "\r\n");
        let blocks = read_blocks(&text).unwrap();
        assert_eq!(blocks.len(), 3);
        assert_eq!(blocks[1].sequences.len(), 5);
        assert_eq!(blocks[0].sequences[0].text.last(), Some(&b'G'));
    }

    #[test]
    fn test_whitespace_only_separator() {
        let text = "a score=1\ns a.c 0 2 + 10 AC\n   \na score=2\ns a.c 2 2 + 10 GT\n";
        assert_eq!(read_blocks(text).unwrap().len(), 2);
    }

    #[test]
    fn test_blocks_without_blank_line() {
        let text = "a score=1\ns a.c 0 2 + 10 AC\na score=2\ns a.c 2 2 + 10 GT";
        let blocks = read_blocks(text).unwrap();
        assert_eq!(blocks.len(), 2);
        assert_eq!(blocks[1].score, Some(2.0));
        assert_eq!(blocks[1].sequences[0].text, b"GT".to_vec());
    }

    #[test]
    fn test_negative_and_exponent_scores() {
        let mut alignment = MafAlignment::new();
        alignment.parse_scores("a score=-12.5 pass=2");
        assert_eq!(alignment.score, Some(-12.5));
        alignment.parse_scores("a pass=2 score=1e3");
        assert_eq!(alignment.score, Some(1000.0));
        alignment.parse_scores("a");
        assert_eq!(alignment.score, None);
    }

    #[test]
    fn test_invalid_sequence_lines() {
        let missing = read_blocks("a score=1\ns a.c 0 2 + AC\n").unwrap_err();
        assert!(format!("{:#}", missing).contains("line 2"));
        let strand = read_blocks("a score=1\ns a.c 0 2 x 10 AC\n").unwrap_err();
        assert!(format!("{:#}", strand).contains("invalid strand"));
        assert!(read_blocks("a score=1\ns a.c zero 2 + 10 AC\n").is_err());
    }

    #[test]
    fn test_track_line_without_newline() {
        let mut track = TrackLine::new();
        track.from_str("track name=chr1").unwrap();
        assert_eq!(track.name, "chr1");
        let blocks = read_blocks("track name=x\n\na score=1\ns a.c 0 2 + 10 AC").unwrap();
        assert_eq!(blocks.len(), 1);
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
