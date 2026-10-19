//! Filter MAF blocks before conversion.
use std::path::Path;

use anyhow::{Context, Result};

use super::block::RefBlock;

/// Criteria a MAF block must meet to be converted.
/// Blocks are checked in this order: species, length, score.
#[derive(Debug, Default, Clone)]
pub struct BlockFilter {
    /// Minimum number of reference bases (gaps excluded).
    pub min_ref_len: Option<usize>,
    /// Minimum number of species, counting the reference.
    pub min_species: Option<usize>,
    /// Minimum block score. Blocks without a score pass.
    pub min_score: Option<f64>,
    /// Write only these species. The reference is still used for coordinates.
    pub species: Option<Vec<String>>,
}

impl BlockFilter {
    /// Parse a comma-separated species list, ignoring blank entries.
    /// # Example
    /// ```
    /// use segul_core::core::maf::BlockFilter;
    ///
    /// assert_eq!(BlockFilter::parse_species_list(" hg18, mm4 ,,"), ["hg18", "mm4"]);
    /// ```
    pub fn parse_species_list(list: &str) -> Vec<String> {
        list.split(',')
            .map(|name| name.trim())
            .filter(|name| !name.is_empty())
            .map(|name| name.to_string())
            .collect()
    }

    /// Read one species per line, ignoring blank lines and lines starting with `#`.
    pub fn read_species_file(path: &Path) -> Result<Vec<String>> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("Failed reading species file {}", path.display()))?;
        Ok(text
            .lines()
            .map(|line| line.trim())
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .map(|line| line.to_string())
            .collect())
    }

    pub(crate) fn passes(&self, block: &RefBlock, ref_sample: &str) -> bool {
        let species = block.rows.len() + usize::from(!block.rows.contains_key(ref_sample));
        if self.min_species.is_some_and(|min| species < min) {
            return false;
        }
        if self
            .min_ref_len
            .is_some_and(|min| block.end - block.start < min)
        {
            return false;
        }
        match (self.min_score, block.score) {
            (Some(min), Some(score)) => score >= min,
            _ => true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use indexmap::IndexMap;

    fn block(samples: &[&str], len: usize, score: Option<f64>) -> RefBlock {
        let rows: IndexMap<String, Vec<u8>> = samples
            .iter()
            .map(|s| (s.to_string(), vec![b'A'; len]))
            .collect();
        RefBlock {
            source: "ref.chr1".to_string(),
            chrom: "chr1".to_string(),
            start: 0,
            end: len,
            src_size: 100,
            score,
            ref_text: vec![b'A'; len],
            rows,
            duplicates: 0,
        }
    }

    #[test]
    fn test_each_filter() {
        let b = block(&["ref", "a"], 10, Some(5.0));
        assert!(BlockFilter::default().passes(&b, "ref"));
        let min_species = BlockFilter {
            min_species: Some(3),
            ..Default::default()
        };
        assert!(!min_species.passes(&b, "ref"));
        let min_len = BlockFilter {
            min_ref_len: Some(11),
            ..Default::default()
        };
        assert!(!min_len.passes(&b, "ref"));
        let min_score = BlockFilter {
            min_score: Some(6.0),
            ..Default::default()
        };
        assert!(!min_score.passes(&b, "ref"));
        assert!(min_score.passes(&block(&["ref"], 10, None), "ref"));
    }

    #[test]
    fn test_reference_counts_when_not_written() {
        let b = block(&["a", "b"], 10, None);
        let filter = BlockFilter {
            min_species: Some(3),
            ..Default::default()
        };
        assert!(filter.passes(&b, "ref"));
    }
}
