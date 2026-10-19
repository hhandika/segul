//! Simple BED parser.
//!
//! Supports BED3 to BED12 records separated by tabs or whitespace.
//! Blank lines, comments (`#`), and `track` or `browser` lines are skipped.
//! Gzip compressed files (`.bed.gz`) are read directly.

use std::{io::BufRead, path::Path};

use anyhow::{Context, bail};
use serde::Deserialize;

use crate::helper::{files, types::DnaStrand};

#[derive(Debug, Deserialize, Default, Clone, PartialEq)]
pub struct BedRecord {
    pub chrom: String,
    pub chrom_start: usize,
    pub chrom_end: usize,
    #[serde(default)]
    pub name: Option<String>,
    pub score: Option<u16>,
    pub strand: Option<DnaStrand>,
    pub thick_start: Option<usize>,
    pub thick_end: Option<usize>,
    pub item_rgb: Option<String>,
    pub block_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub block_sizes: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub block_starts: Option<String>,
}

impl BedRecord {
    pub fn new(chrom: String, chrom_start: usize, chrom_end: usize, name: Option<String>) -> Self {
        Self {
            chrom,
            chrom_start,
            chrom_end,
            name,
            ..Self::default()
        }
    }

    /// Parse a BED record from a single line.
    ///
    /// The line may end with `\n` or `\r\n`. Fields are split on tabs when the line has
    /// any, so names may contain spaces; otherwise they are split on whitespace.
    /// The first three columns are required. Optional columns are read leniently:
    /// a missing or invalid value, such as `.`, becomes `None`.
    pub fn from_vec_bytes(line: &[u8]) -> anyhow::Result<Self> {
        let line = line.trim_ascii_end();
        let fields: Vec<&[u8]> = if line.contains(&b'\t') {
            line.split(|&b| b == b'\t')
                .map(|f| f.trim_ascii())
                .collect()
        } else {
            line.split(|b| b.is_ascii_whitespace())
                .filter(|f| !f.is_empty())
                .collect()
        };
        if fields.len() < 3 {
            bail!(
                "expected at least 3 columns (chrom, start, end), found {}",
                fields.len()
            );
        }
        let chrom = Self::to_text(fields[0]);
        if chrom.is_empty() {
            bail!("missing chromosome name");
        }
        let chrom_start = Self::parse_position(fields[1], "start")?;
        let chrom_end = Self::parse_position(fields[2], "end")?;
        if chrom_start >= chrom_end {
            bail!(
                "start ({}) must be smaller than end ({})",
                chrom_start,
                chrom_end
            );
        }
        let optional = |i: usize| fields.get(i).filter(|f| !f.is_empty() && **f != b".");
        Ok(Self {
            chrom,
            chrom_start,
            chrom_end,
            name: optional(3).map(|f| Self::to_text(f)),
            score: optional(4).and_then(|f| Self::to_text(f).parse().ok()),
            strand: fields
                .get(5)
                .map(|f| Self::to_text(f).parse().unwrap_or(DnaStrand::Missing)),
            thick_start: optional(6).and_then(|f| Self::to_text(f).parse().ok()),
            thick_end: optional(7).and_then(|f| Self::to_text(f).parse().ok()),
            item_rgb: optional(8).map(|f| Self::to_text(f)),
            block_count: optional(9).and_then(|f| Self::to_text(f).parse().ok()),
            block_sizes: optional(10).map(|f| Self::to_text(f)),
            block_starts: optional(11).map(|f| Self::to_text(f)),
        })
    }

    fn to_text(field: &[u8]) -> String {
        String::from_utf8_lossy(field).to_string()
    }

    fn parse_position(field: &[u8], column: &str) -> anyhow::Result<usize> {
        let text = String::from_utf8_lossy(field);
        text.parse::<usize>()
            .with_context(|| format!("invalid {} position '{}'", column, text))
    }
}

pub struct BedParser<'a> {
    pub file: &'a Path,
    /// Kept for compatibility. Headers are always detected and skipped.
    pub has_header: bool,
}

impl<'a> BedParser<'a> {
    pub fn new(file: &'a Path) -> Self {
        Self {
            file,
            has_header: false,
        }
    }

    /// Parse every record in the file. Gzip files are detected by extension.
    pub fn parse(&self) -> anyhow::Result<Vec<BedRecord>> {
        let reader = files::open_reader(self.file)?;
        Self::parse_reader(reader)
            .with_context(|| format!("Failed parsing BED file {}", self.file.display()))
    }

    /// Parse every record from a buffered reader.
    pub fn parse_reader<R: BufRead>(mut reader: R) -> anyhow::Result<Vec<BedRecord>> {
        let mut records = Vec::new();
        let mut line = Vec::new();
        let mut line_number = 0;
        loop {
            line.clear();
            if reader.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            line_number += 1;
            let content = line.trim_ascii();
            if Self::is_skipped(content) {
                continue;
            }
            let record = BedRecord::from_vec_bytes(content)
                .with_context(|| format!("Invalid BED record at line {}", line_number))?;
            records.push(record);
        }
        Ok(records)
    }

    fn is_skipped(line: &[u8]) -> bool {
        line.is_empty()
            || line.starts_with(b"#")
            || line.starts_with(b"track")
            || line.starts_with(b"browser")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::path::PathBuf;

    use flate2::{Compression, write::GzEncoder};

    fn parse_text(text: &str) -> anyhow::Result<Vec<BedRecord>> {
        BedParser::parse_reader(text.as_bytes())
    }

    #[test]
    fn test_bed_parser() {
        let bed_file = PathBuf::from("tests/files/bed/common.bed");
        let records = BedParser::new(&bed_file).parse().unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[0].chrom_start, 1000);
        assert_eq!(records[0].chrom_end, 5000);
        assert_eq!(records[0].name, Some("gene1".to_string()));
        assert_eq!(records[0].score, Some(0));
        assert_eq!(records[0].strand, Some(DnaStrand::Forward));
    }

    #[test]
    fn test_bed_with_header() {
        let bed_file = PathBuf::from("tests/files/bed/with_header.bed");
        let records = BedParser::new(&bed_file).parse().unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[0].name, Some("gene1".to_string()));
        assert_eq!(records[1].strand, Some(DnaStrand::Reverse));
    }

    #[test]
    fn test_four_columns_with_trailing_newline() {
        let records = parse_text("chr1\t100\t200\tlocus_a\nchr1\t300\t400\tlocus_b\n").unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[1].name, Some("locus_b".to_string()));
        assert_eq!(records[1].score, None);
        assert_eq!(records[1].strand, None);
    }

    #[test]
    fn test_crlf_and_blank_lines() {
        let records = parse_text("chr1 100 200 a\r\n\r\n   \r\nchr2 5 10 b\r\n").unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, Some("a".to_string()));
        assert_eq!(records[1].chrom, "chr2");
    }

    #[test]
    fn test_bed3() {
        let records = parse_text("chr1\t0\t10").unwrap();
        assert_eq!(records[0].name, None);
        assert_eq!(records[0].chrom_end, 10);
    }

    #[test]
    fn test_tab_delimited_name_with_spaces() {
        let records = parse_text("chr1\t0\t10\tgene one\t0\t-\n").unwrap();
        assert_eq!(records[0].name, Some("gene one".to_string()));
        assert_eq!(records[0].strand, Some(DnaStrand::Reverse));
    }

    #[test]
    fn test_lenient_optional_columns() {
        let records = parse_text("chr1\t0\t10\tg\t1.5\t.\n").unwrap();
        assert_eq!(records[0].score, None);
        assert_eq!(records[0].strand, Some(DnaStrand::Missing));
    }

    #[test]
    fn test_short_line_reports_line_number() {
        let err = parse_text("chr1\t0\t10\n\nchr1\t5\n").unwrap_err();
        let message = format!("{:#}", err);
        assert!(message.contains("line 3"), "{}", message);
        assert!(message.contains("at least 3 columns"), "{}", message);
    }

    #[test]
    fn test_invalid_positions() {
        assert!(parse_text("chr1\tx\t10\n").is_err());
        assert!(parse_text("chr1\t10\t10\n").is_err());
        assert!(parse_text("chr1\t20\t10\n").is_err());
    }

    #[test]
    fn test_gzip_bed() {
        let dir = std::env::temp_dir().join("segul-test-bed-gz");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("regions.bed.gz");
        let mut encoder =
            GzEncoder::new(std::fs::File::create(&path).unwrap(), Compression::fast());
        encoder.write_all(b"chr1\t0\t10\tg\n").unwrap();
        encoder.finish().unwrap();
        let records = BedParser::new(&path).parse().unwrap();
        assert_eq!(records.len(), 1);
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
