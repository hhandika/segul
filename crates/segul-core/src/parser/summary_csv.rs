//! Reader for the summary CSV files written by segul.
//!
//! Header names are trimmed because some writers emit a space after the
//! comma. Empty or non-numeric cells read as NaN in numeric columns, which
//! covers undefined ratios in the MAF and distance summaries.
use std::fs::File;
use std::io::Read;
use std::path::Path;

use anyhow::{Context, Result, bail};
use csv::{ReaderBuilder, Trim};
use zip::ZipArchive;

/// A summary table held in memory as strings.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SummaryCsv {
    headers: Vec<String>,
    records: Vec<Vec<String>>,
}

impl SummaryCsv {
    /// Read a comma-delimited summary file.
    pub fn from_path(path: &Path) -> Result<Self> {
        let file =
            File::open(path).with_context(|| format!("Failed opening {}", path.display()))?;
        Self::from_reader(file).with_context(|| format!("Failed reading {}", path.display()))
    }

    /// Read every CSV entry in a zip archive (e.g. `read-pos-summary.zip`).
    /// Returns the entry file stem with its table, in archive order.
    pub fn from_zip(path: &Path) -> Result<Vec<(String, Self)>> {
        let file =
            File::open(path).with_context(|| format!("Failed opening {}", path.display()))?;
        let mut archive = ZipArchive::new(file)
            .with_context(|| format!("Failed reading zip archive {}", path.display()))?;
        let mut tables = Vec::with_capacity(archive.len());
        for i in 0..archive.len() {
            let mut entry = archive.by_index(i)?;
            let entry_name = entry.name()?.to_string();
            if !entry.is_file() || !entry_name.ends_with(".csv") {
                continue;
            }
            let name = Path::new(&entry_name)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or_default()
                .to_string();
            let mut text = String::new();
            entry
                .read_to_string(&mut text)
                .with_context(|| format!("Failed reading {} in {}", name, path.display()))?;
            tables.push((name, Self::from_reader(text.as_bytes())?));
        }
        Ok(tables)
    }

    pub fn from_reader<R: Read>(reader: R) -> Result<Self> {
        let mut reader = ReaderBuilder::new()
            .trim(Trim::All)
            .flexible(true)
            .from_reader(reader);
        let headers = reader.headers()?.iter().map(|h| h.to_string()).collect();
        let records = reader
            .records()
            .map(|r| r.map(|r| r.iter().map(|c| c.to_string()).collect()))
            .collect::<Result<Vec<Vec<String>>, _>>()?;
        Ok(Self { headers, records })
    }

    /// Number of data rows.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn headers(&self) -> &[String] {
        &self.headers
    }

    pub fn has_column(&self, name: &str) -> bool {
        self.headers.iter().any(|h| h == name)
    }

    /// Column values as strings. Missing cells read as empty strings.
    pub fn column_str(&self, name: &str) -> Result<Vec<String>> {
        let idx = self.column_index(name)?;
        Ok(self.column_at(idx))
    }

    /// Column values as numbers. Empty or invalid cells read as NaN.
    pub fn column_f64(&self, name: &str) -> Result<Vec<f64>> {
        let idx = self.column_index(name)?;
        Ok(self
            .records
            .iter()
            .map(|r| {
                r.get(idx)
                    .and_then(|c| c.parse::<f64>().ok())
                    .unwrap_or(f64::NAN)
            })
            .collect())
    }

    /// Column values by position, as strings.
    pub fn column_at(&self, idx: usize) -> Vec<String> {
        self.records
            .iter()
            .map(|r| r.get(idx).cloned().unwrap_or_default())
            .collect()
    }

    fn column_index(&self, name: &str) -> Result<usize> {
        match self.headers.iter().position(|h| h == name) {
            Some(idx) => Ok(idx),
            None => bail!(
                "Column '{}' not found. Available columns: {}",
                name,
                self.headers.join(", ")
            ),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_trimmed_headers() {
        let csv = SummaryCsv::from_path(Path::new("tests/files/plot/taxon_summary.csv")).unwrap();
        assert!(csv.has_column("total_chars"));
        assert_eq!(csv.len(), 3);
        assert_eq!(csv.column_f64("locus_counts").unwrap(), vec![4.0, 2.0, 1.0]);
    }

    #[test]
    fn test_empty_cells_are_nan() {
        let csv =
            SummaryCsv::from_path(Path::new("tests/files/plot/distance_summary.csv")).unwrap();
        let mean = csv.column_f64("mean").unwrap();
        assert!(mean[0].is_nan());
        assert_eq!(mean[1], 0.125);
    }

    #[test]
    fn test_missing_column() {
        let csv = SummaryCsv::from_reader("a,b\n1,2\n".as_bytes()).unwrap();
        assert!(csv.column_f64("c").is_err());
        assert_eq!(csv.column_str("b").unwrap(), vec!["2"]);
    }

    #[test]
    fn test_from_zip() {
        let tables =
            SummaryCsv::from_zip(Path::new("tests/files/plot/read-pos-summary.zip")).unwrap();
        assert_eq!(tables.len(), 4);
        let (name, table) = tables.iter().find(|(n, _)| n == "read_1").unwrap();
        assert_eq!(name, "read_1");
        assert!(table.has_column("mean_QScore"));
        assert!(!table.is_empty());
    }
}
