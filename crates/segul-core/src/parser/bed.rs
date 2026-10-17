//! Simple BED parser.

use std::{error::Error, io::BufRead, path::Path};

use serde::Deserialize;

use crate::helper::types::DnaStrand;

#[cfg(target_os = "windows")]
use super::CAR_RETURN;
use super::END_OF_LINE;

#[derive(Debug, Deserialize, Default)]
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
            score: None,
            strand: None,
            thick_start: None,
            thick_end: None,
            item_rgb: None,
            block_count: None,
            block_sizes: None,
            block_starts: None,
        }
    }

    pub fn from_vec_bytes(line: &[u8]) -> Result<Self, Box<dyn Error>> {
        let mut record = Self::default();
        let mut parts = line.split(|&b| b.is_ascii_whitespace());
        record.chrom = String::from_utf8_lossy(parts.next().unwrap()).to_string();
        record.chrom_start = String::from_utf8_lossy(parts.next().unwrap()).parse::<usize>()?;
        record.chrom_end = String::from_utf8_lossy(parts.next().unwrap()).parse::<usize>()?;
        if let Some(name) = parts.next() {
            record.name = Some(String::from_utf8_lossy(name).to_string());
        }
        if let Some(score) = parts.next() {
            record.score = Some(String::from_utf8_lossy(score).parse::<u16>()?);
        }
        if let Some(strand) = parts.next() {
            record.strand = Some(
                String::from_utf8_lossy(strand)
                    .parse::<DnaStrand>()
                    .unwrap_or_default(),
            );
        }
        if let Some(thick_start) = parts.next() {
            record.thick_start = Some(
                String::from_utf8_lossy(thick_start)
                    .parse::<usize>()
                    .unwrap_or_default(),
            );
        }
        if let Some(thick_end) = parts.next() {
            record.thick_end = Some(
                String::from_utf8_lossy(thick_end)
                    .parse::<usize>()
                    .unwrap_or_default(),
            );
        }
        if let Some(item_rgb) = parts.next() {
            record.item_rgb = Some(String::from_utf8_lossy(item_rgb).to_string());
        }
        if let Some(block_count) = parts.next() {
            record.block_count = Some(
                String::from_utf8_lossy(block_count)
                    .parse::<usize>()
                    .unwrap_or_default(),
            );
        }
        if let Some(block_sizes) = parts.next() {
            record.block_sizes = Some(String::from_utf8_lossy(block_sizes).to_string());
        }
        if let Some(block_starts) = parts.next() {
            record.block_starts = Some(String::from_utf8_lossy(block_starts).to_string());
        }

        if record.chrom.is_empty() || record.chrom_start >= record.chrom_end {
            return Err("Invalid BED record".into());
        }
        Ok(record)
    }
}

pub struct BedParser<'a> {
    pub file: &'a Path,
    pub has_header: bool,
}

impl<'a> BedParser<'a> {
    pub fn new(file: &'a Path) -> Self {
        Self {
            file,
            has_header: false,
        }
    }

    pub fn parse(&mut self) -> Result<Vec<BedRecord>, Box<dyn Error>> {
        self.parse_with_header()
    }

    fn parse_with_header(&self) -> Result<Vec<BedRecord>, Box<dyn Error>> {
        let file = std::fs::File::open(self.file)?;
        let mut buf = std::io::BufReader::new(file);
        let mut bed = Vec::new();
        let mut line = Vec::new();
        loop {
            line.clear();
            let bytes_read = buf.read_until(END_OF_LINE, &mut line)?;
            if bytes_read == 0 {
                break;
            }

            if line[0] == b'#' {
                continue;
            }
            if line.starts_with(b"track") || line.starts_with(b"browser") {
                continue;
            }

            #[cfg(target_os = "windows")]
            if line.contains(&CAR_RETURN) {
                // Remove the carriage return
                // and leave only the line feed
                line.retain(|&c| c != CAR_RETURN);
            }

            let record = BedRecord::from_vec_bytes(&line)
                .map_err(|e| format!("Failed to parse BED record: {}", e))?;
            bed.push(record);
        }

        Ok(bed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_bed_parser() {
        let bed_file = PathBuf::from("tests/files/bed/common.bed");
        let mut parser = BedParser::new(&bed_file);
        let records = parser.parse().unwrap();
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
        let mut parser = BedParser::new(&bed_file);
        let records = parser.parse().unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[0].chrom_start, 1000);
        assert_eq!(records[0].chrom_end, 5000);
        assert_eq!(records[0].name, Some("gene1".to_string()));
        assert_eq!(records[0].score, Some(0));
        assert_eq!(records[0].strand, Some(DnaStrand::Forward));
    }
}
