//! Simple BED parser.

use std::{error::Error, path::Path};

use csv::ReaderBuilder;
use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct BetRecord {
    pub chrom: String,
    pub chrom_start: usize,
    pub chrom_end: usize,
    pub name: Option<String>,
    pub score: Option<String>,
    pub strand: Option<String>,
    pub thick_start: Option<usize>,
    pub thick_end: Option<usize>,
    pub item_rgb: Option<String>,
    pub block_count: Option<usize>,
    pub block_sizes: Option<String>,
    pub block_starts: Option<String>,
}

impl BetRecord {
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
}

pub struct BedParser<'a> {
    pub file: &'a Path,
    pub has_header: bool,
}

impl<'a> BedParser<'a> {
    pub fn new(file: &'a Path, has_header: bool) -> Self {
        Self { file, has_header }
    }

    pub fn parser(&self) -> Result<Vec<BetRecord>, Box<dyn Error>> {
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(self.has_header)
            .from_path(self.file)?;

        let mut bed = Vec::new();
        for result in reader.deserialize() {
            let record: BetRecord = result?;
            bed.push(record);
        }
        Ok(bed)
    }
}
