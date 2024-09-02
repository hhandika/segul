//! BED parser. Partially implemented.

use std::{error::Error, path::Path};

use csv::ReaderBuilder;
use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct BetRecord {
    pub chrom: String,
    pub chrom_start: usize,
    pub chrom_end: usize,
    pub name: Option<String>,
}

impl BetRecord {
    pub fn new(chrom: String, chrom_start: usize, chrom_end: usize, name: Option<String>) -> Self {
        Self {
            chrom,
            chrom_start,
            chrom_end,
            name,
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
