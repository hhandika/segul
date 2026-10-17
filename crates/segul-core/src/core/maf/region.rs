//! Index BED regions for overlap queries.
use std::collections::{HashMap, HashSet};

use crate::{helper::types::DnaStrand, parser::bed::BedRecord};

/// A BED region on the forward strand. Coordinates are zero-based and half-open.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Region {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    pub strand: DnaStrand,
    /// Output name, e.g. `gene1-chr1-100-200`.
    pub label: String,
}

/// Regions grouped by chromosome and sorted by start position.
pub(crate) struct RegionIndex {
    regions: Vec<Region>,
    // Chromosome -> (first index, last index exclusive, longest region length).
    by_chrom: HashMap<String, (usize, usize, usize)>,
}

impl RegionIndex {
    /// Build the index. Exact duplicate records are removed, and records that
    /// produce an existing name are skipped, both with a warning.
    pub(crate) fn from_bed(records: &[BedRecord]) -> Self {
        let mut regions: Vec<Region> = records
            .iter()
            .map(|record| Region {
                chrom: record.chrom.clone(),
                start: record.chrom_start,
                end: record.chrom_end,
                strand: record.strand.unwrap_or(DnaStrand::Missing),
                label: format_label(record),
            })
            .collect();
        regions.sort_by(|a, b| {
            (&a.chrom, a.start, a.end, &a.label).cmp(&(&b.chrom, b.start, b.end, &b.label))
        });
        let before = regions.len();
        regions.dedup();
        if regions.len() < before {
            log::warn!("Skipped {} duplicate BED records", before - regions.len());
        }
        let mut labels = HashSet::new();
        regions.retain(|region| {
            let is_new = labels.insert(region.label.clone());
            if !is_new {
                log::warn!("Skipped BED record with a duplicate name: {}", region.label);
            }
            is_new
        });
        let mut by_chrom: HashMap<String, (usize, usize, usize)> = HashMap::new();
        for (idx, region) in regions.iter().enumerate() {
            let length = region.end - region.start;
            let entry = by_chrom
                .entry(region.chrom.clone())
                .or_insert((idx, idx, 0));
            entry.1 = idx + 1;
            entry.2 = entry.2.max(length);
        }
        Self { regions, by_chrom }
    }

    pub(crate) fn regions(&self) -> &[Region] {
        &self.regions
    }

    /// Ids of regions that overlap `[start, end)` on a chromosome.
    /// `fallback` is tried when `chrom` has no regions, e.g. the full MAF source.
    pub(crate) fn overlapping(
        &self,
        chrom: &str,
        fallback: &str,
        start: usize,
        end: usize,
    ) -> Vec<usize> {
        let Some(&(first, last, max_len)) = self
            .by_chrom
            .get(chrom)
            .or_else(|| self.by_chrom.get(fallback))
        else {
            return Vec::new();
        };
        let regions = &self.regions[first..last];
        // Regions starting this far before the query cannot reach it.
        let lo = regions.partition_point(|r| r.start + max_len <= start);
        let hi = regions.partition_point(|r| r.start < end);
        (lo..hi.max(lo))
            .filter(|&i| regions[i].end > start)
            .map(|i| first + i)
            .collect()
    }
}

/// Output name for a BED record: `{name}-{chrom}-{start}-{end}`, or
/// `{chrom}-{start}-{end}` without a name. Characters that are unsafe in
/// file names, including dots, become underscores.
pub(crate) fn format_label(record: &BedRecord) -> String {
    let label = match &record.name {
        Some(name) => format!(
            "{}-{}-{}-{}",
            name, record.chrom, record.chrom_start, record.chrom_end
        ),
        None => format!(
            "{}-{}-{}",
            record.chrom, record.chrom_start, record.chrom_end
        ),
    };
    sanitize_name(&label)
}

/// Replace characters that are unsafe in file names, and dots, with underscores.
pub(crate) fn sanitize_name(name: &str) -> String {
    name.chars()
        .map(|c| match c {
            '.' | '/' | '\\' | ':' | '*' | '?' | '"' | '<' | '>' | '|' => '_',
            c if c.is_whitespace() => '_',
            c => c,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(chrom: &str, start: usize, end: usize, name: &str) -> BedRecord {
        BedRecord::new(chrom.to_string(), start, end, Some(name.to_string()))
    }

    #[test]
    fn test_half_open_overlap() {
        let index = RegionIndex::from_bed(&[record("chr1", 10, 20, "a")]);
        assert!(index.overlapping("chr1", "x", 0, 10).is_empty());
        assert!(index.overlapping("chr1", "x", 20, 30).is_empty());
        assert_eq!(index.overlapping("chr1", "x", 19, 30), [0]);
        assert_eq!(index.overlapping("chr1", "x", 0, 11), [0]);
        assert!(index.overlapping("chr2", "x", 10, 20).is_empty());
    }

    #[test]
    fn test_long_early_region() {
        let mut records = vec![record("chr1", 0, 1000, "long")];
        records.extend((0..50).map(|i| record("chr1", 10 + i * 10, 15 + i * 10, &format!("s{i}"))));
        let index = RegionIndex::from_bed(&records);
        let labels = |hits: Vec<usize>| -> Vec<String> {
            hits.iter()
                .map(|&i| index.regions()[i].label.clone())
                .collect()
        };
        // Only the long region reaches past the short ones.
        assert_eq!(
            labels(index.overlapping("chr1", "x", 900, 910)),
            ["long-chr1-0-1000"]
        );
        // Short region s48 is [490, 495).
        assert_eq!(
            labels(index.overlapping("chr1", "x", 491, 492)),
            ["long-chr1-0-1000", "s48-chr1-490-495"]
        );
    }

    #[test]
    fn test_one_block_matches_many_regions_and_fallback() {
        let index = RegionIndex::from_bed(&[
            record("chr7", 0, 10, "a"),
            record("chr7", 5, 15, "b"),
            record("hg18.chrX", 0, 10, "c"),
        ]);
        assert_eq!(index.overlapping("chr7", "hg18.chr7", 8, 9).len(), 2);
        assert_eq!(index.overlapping("chrX", "hg18.chrX", 0, 1).len(), 1);
    }

    #[test]
    fn test_duplicates_and_labels() {
        let index = RegionIndex::from_bed(&[
            record("chr1", 0, 10, "g.1"),
            record("chr1", 0, 10, "g.1"),
            BedRecord::new("chr1".to_string(), 20, 30, None),
        ]);
        let labels: Vec<&str> = index.regions().iter().map(|r| r.label.as_str()).collect();
        assert_eq!(labels, ["g_1-chr1-0-10", "chr1-20-30"]);
        assert_eq!(sanitize_name("a b/c:d"), "a_b_c_d");
    }
}
