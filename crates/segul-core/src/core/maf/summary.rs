//! Summarize MAF files: blocks, species, reference coverage, and gaps.
//!
//! Writes two CSV files to the output directory:
//! - `{prefix}_per_file.csv`: one row per MAF file.
//! - `{prefix}_per_species.csv`: one row per species in each file. With more than
//!   one input file, rows with `path` set to `all` combine every file.
use std::{
    collections::HashMap,
    fs,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result};
use colored::Colorize;
use indexmap::IndexMap;
use rayon::prelude::*;

use super::block::{self, RefBlock, SampleNaming, SpeciesResolver};
use crate::{
    helper::{files, utils},
    parser::maf::{MafParagraph, MafReader},
};

/// Default prefix of the summary files.
pub const DEFAULT_SUMMARY_PREFIX: &str = "maf_summary";
// Path label for rows that combine every input file.
const ALL_FILES: &str = "all";

#[derive(Debug, Default, Clone)]
struct FileStats {
    path: String,
    blocks: usize,
    blocks_no_ref: usize,
    ref_chroms: HashMap<String, usize>,
    ref_bases: usize,
    aligned_columns: usize,
    min_block_columns: Option<usize>,
    max_block_columns: usize,
    species_in_blocks: usize,
    score_sum: f64,
    scored_blocks: usize,
    gaps: usize,
    cells: usize,
    duplicate_rows: usize,
    species: IndexMap<String, SpeciesStats>,
}

#[derive(Debug, Default, Clone)]
struct SpeciesStats {
    blocks_present: usize,
    aligned_columns: usize,
    bases: usize,
    gaps: usize,
    n_count: usize,
    ref_bases_covered: usize,
}

impl SpeciesStats {
    fn add(&mut self, other: &SpeciesStats) {
        self.blocks_present += other.blocks_present;
        self.aligned_columns += other.aligned_columns;
        self.bases += other.bases;
        self.gaps += other.gaps;
        self.n_count += other.n_count;
        self.ref_bases_covered += other.ref_bases_covered;
    }
}

impl FileStats {
    fn add_block(&mut self, block: &RefBlock) {
        let ncols = block.ncols();
        self.blocks += 1;
        self.ref_chroms.insert(block.source.clone(), block.src_size);
        self.ref_bases += block.end - block.start;
        self.aligned_columns += ncols;
        self.min_block_columns = Some(self.min_block_columns.map_or(ncols, |m| m.min(ncols)));
        self.max_block_columns = self.max_block_columns.max(ncols);
        self.species_in_blocks += block.rows.len();
        if let Some(score) = block.score {
            self.score_sum += score;
            self.scored_blocks += 1;
        }
        self.duplicate_rows += block.duplicates;
        for (sample, row) in &block.rows {
            let gaps = row.iter().filter(|&&b| b == b'-').count();
            let covered = row
                .iter()
                .zip(&block.ref_text)
                .filter(|&(&b, &r)| b != b'-' && r != b'-')
                .count();
            let stats = self.species.entry(sample.clone()).or_default();
            stats.blocks_present += 1;
            stats.aligned_columns += ncols;
            stats.bases += ncols - gaps;
            stats.gaps += gaps;
            stats.n_count += row
                .iter()
                .filter(|b| matches!(b, b'N' | b'n' | b'?'))
                .count();
            stats.ref_bases_covered += covered;
            self.gaps += gaps;
            self.cells += ncols;
        }
    }

    fn ref_chrom_size_total(&self) -> usize {
        self.ref_chroms.values().sum()
    }
}

/// Summarize MAF files.
///
/// # Example
/// ```
/// use std::path::PathBuf;
/// use segul_core::core::maf::summary::MafSummary;
///
/// let output = tempfile::TempDir::with_prefix("maf-summary").unwrap();
/// let files = vec![PathBuf::from("tests/files/maf/simple.maf")];
/// MafSummary::new(&files, output.path(), None).summarize().unwrap();
/// assert!(output.path().join("maf_summary_per_file.csv").exists());
/// assert!(output.path().join("maf_summary_per_species.csv").exists());
/// ```
pub struct MafSummary<'a> {
    input_files: &'a [PathBuf],
    output_dir: &'a Path,
    prefix: Option<&'a str>,
    naming: SampleNaming,
    ref_species: Option<String>,
    species: Option<Vec<String>>,
}

impl<'a> MafSummary<'a> {
    pub fn new(input_files: &'a [PathBuf], output_dir: &'a Path, prefix: Option<&'a str>) -> Self {
        Self {
            input_files,
            output_dir,
            prefix,
            naming: SampleNaming::default(),
            ref_species: None,
            species: None,
        }
    }

    pub fn with_naming(mut self, naming: SampleNaming) -> Self {
        self.naming = naming;
        self
    }

    /// Reference species. Defaults to the species of the first sequence in the first block.
    pub fn with_ref_species(mut self, ref_species: Option<String>) -> Self {
        self.ref_species = ref_species;
        self
    }

    /// Report only these species. Also names sources whose species contain dots.
    pub fn with_species(mut self, species: Option<Vec<String>>) -> Self {
        self.species = species;
        self
    }

    /// Write the summary files. Returns their paths: per file, then per species.
    pub fn summarize(&self) -> Result<(PathBuf, PathBuf)> {
        let mut input_files = self.input_files.to_vec();
        input_files.sort();
        let resolver = SpeciesResolver::new(self.naming, self.species.as_deref());
        let ref_sample = match &self.ref_species {
            Some(species) => species.clone(),
            None => block::first_sample(&input_files, &resolver)?,
        };
        let spin = utils::set_spinner();
        spin.set_message("Summarizing MAF files...");
        let stats = input_files
            .par_iter()
            .map(|file| self.summarize_file(file, &resolver, &ref_sample))
            .collect::<Result<Vec<FileStats>>>()?;
        fs::create_dir_all(self.output_dir).with_context(|| {
            format!(
                "Failed creating output directory {}",
                self.output_dir.display()
            )
        })?;
        let prefix = self.prefix.unwrap_or(DEFAULT_SUMMARY_PREFIX);
        let per_file = self.output_dir.join(format!("{}_per_file.csv", prefix));
        let per_species = self.output_dir.join(format!("{}_per_species.csv", prefix));
        self.write_per_file(&per_file, &stats, &ref_sample)?;
        self.write_per_species(&per_species, &stats, &ref_sample)?;
        spin.finish_with_message("Finished summarizing MAF files!\n");
        self.print_output_info(&stats, &ref_sample, &per_file, &per_species);
        Ok((per_file, per_species))
    }

    fn summarize_file(
        &self,
        file: &Path,
        resolver: &SpeciesResolver,
        ref_sample: &str,
    ) -> Result<FileStats> {
        let mut stats = FileStats {
            path: file.display().to_string(),
            ..Default::default()
        };
        let reader = MafReader::new(files::open_reader(file)?);
        for paragraph in reader {
            let paragraph =
                paragraph.with_context(|| format!("Failed reading {}", file.display()))?;
            let MafParagraph::Alignment(alignment) = paragraph else {
                continue;
            };
            let block_number = stats.blocks + stats.blocks_no_ref + 1;
            let block = RefBlock::from_alignment(alignment, resolver, ref_sample)
                .with_context(|| format!("Invalid block {} in {}", block_number, file.display()))?;
            match block {
                Some(block) => stats.add_block(&block),
                None => stats.blocks_no_ref += 1,
            }
        }
        Ok(stats)
    }

    fn write_per_file(&self, path: &Path, stats: &[FileStats], ref_sample: &str) -> Result<()> {
        let mut writer = csv::Writer::from_path(path)
            .with_context(|| format!("Failed creating {}", path.display()))?;
        writer.write_record([
            "path",
            "ref_species",
            "blocks",
            "blocks_no_ref",
            "species",
            "ref_chromosomes",
            "ref_chrom_size_total",
            "ref_bases",
            "ref_genome_coverage",
            "aligned_columns",
            "mean_block_columns",
            "min_block_columns",
            "max_block_columns",
            "mean_block_ref_len",
            "mean_species_per_block",
            "mean_score",
            "gap_proportion",
            "duplicate_rows",
        ])?;
        for file in stats {
            let mean_score = if file.scored_blocks > 0 {
                format!("{:.4}", file.score_sum / file.scored_blocks as f64)
            } else {
                String::new()
            };
            writer.write_record([
                file.path.clone(),
                ref_sample.to_string(),
                file.blocks.to_string(),
                file.blocks_no_ref.to_string(),
                file.species.len().to_string(),
                file.ref_chroms.len().to_string(),
                file.ref_chrom_size_total().to_string(),
                file.ref_bases.to_string(),
                ratio(file.ref_bases, file.ref_chrom_size_total()),
                file.aligned_columns.to_string(),
                ratio(file.aligned_columns, file.blocks),
                file.min_block_columns.unwrap_or(0).to_string(),
                file.max_block_columns.to_string(),
                ratio(file.ref_bases, file.blocks),
                ratio(file.species_in_blocks, file.blocks),
                mean_score,
                ratio(file.gaps, file.cells),
                file.duplicate_rows.to_string(),
            ])?;
        }
        writer.flush()?;
        Ok(())
    }

    fn write_per_species(&self, path: &Path, stats: &[FileStats], ref_sample: &str) -> Result<()> {
        let mut writer = csv::Writer::from_path(path)
            .with_context(|| format!("Failed creating {}", path.display()))?;
        writer.write_record([
            "path",
            "name",
            "blocks_present",
            "block_proportion",
            "aligned_columns",
            "bases",
            "gaps",
            "gap_proportion",
            "n_count",
            "ref_bases_covered",
            "ref_coverage",
        ])?;
        let mut rows: Vec<(String, usize, usize, IndexMap<String, SpeciesStats>)> = stats
            .iter()
            .map(|f| (f.path.clone(), f.blocks, f.ref_bases, f.species.clone()))
            .collect();
        if stats.len() > 1 {
            let mut combined: IndexMap<String, SpeciesStats> = IndexMap::new();
            for file in stats {
                for (name, species) in &file.species {
                    combined.entry(name.clone()).or_default().add(species);
                }
            }
            let blocks = stats.iter().map(|f| f.blocks).sum();
            let ref_bases = stats.iter().map(|f| f.ref_bases).sum();
            rows.push((ALL_FILES.to_string(), blocks, ref_bases, combined));
        }
        for (path, blocks, ref_bases, species) in rows {
            for (name, s) in sorted_species(&species, ref_sample) {
                writer.write_record([
                    path.clone(),
                    name.to_string(),
                    s.blocks_present.to_string(),
                    ratio(s.blocks_present, blocks),
                    s.aligned_columns.to_string(),
                    s.bases.to_string(),
                    s.gaps.to_string(),
                    ratio(s.gaps, s.aligned_columns),
                    s.n_count.to_string(),
                    s.ref_bases_covered.to_string(),
                    ratio(s.ref_bases_covered, ref_bases),
                ])?;
            }
        }
        writer.flush()?;
        Ok(())
    }

    fn print_output_info(
        &self,
        stats: &[FileStats],
        ref_sample: &str,
        per_file: &Path,
        per_species: &Path,
    ) {
        let blocks: usize = stats.iter().map(|f| f.blocks).sum();
        let no_ref: usize = stats.iter().map(|f| f.blocks_no_ref).sum();
        log::info!("{}", "Summary".yellow());
        log::info!("{:18}: {}", "Reference", ref_sample);
        log::info!("{:18}: {}", "Files", utils::fmt_num(&stats.len()));
        log::info!("{:18}: {}", "Blocks", utils::fmt_num(&blocks));
        log::info!("{:18}: {}\n", "No reference", utils::fmt_num(&no_ref));
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Per file", per_file.display());
        log::info!("{:18}: {}", "Per species", per_species.display());
    }
}

/// Reference first, then the other species alphabetically.
fn sorted_species<'a>(
    species: &'a IndexMap<String, SpeciesStats>,
    ref_sample: &str,
) -> Vec<(&'a str, &'a SpeciesStats)> {
    let mut sorted: Vec<(&str, &SpeciesStats)> =
        species.iter().map(|(k, v)| (k.as_str(), v)).collect();
    sorted.sort_by(|a, b| (a.0 != ref_sample, a.0).cmp(&(b.0 != ref_sample, b.0)));
    sorted
}

/// Format a ratio with four decimals, or in scientific notation when it is too
/// small to show, such as the coverage of a large genome. Empty when undefined.
fn ratio(numerator: usize, denominator: usize) -> String {
    if denominator == 0 {
        return String::new();
    }
    let value = numerator as f64 / denominator as f64;
    if value > 0.0 && value < 0.0001 {
        format!("{:.4e}", value)
    } else {
        format!("{:.4}", value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    fn read_csv(path: &Path) -> Vec<HashMap<String, String>> {
        let mut reader = csv::Reader::from_path(path).unwrap();
        let headers = reader.headers().unwrap().clone();
        reader
            .records()
            .map(|r| {
                let r = r.unwrap();
                headers
                    .iter()
                    .map(|h| h.to_string())
                    .zip(r.iter().map(|v| v.to_string()))
                    .collect()
            })
            .collect()
    }

    #[test]
    fn test_simple_summary() {
        let output = TempDir::with_prefix("maf-summary").unwrap();
        let files = vec![PathBuf::from("tests/files/maf/simple.maf")];
        let (per_file, per_species) = MafSummary::new(&files, output.path(), Some("simple"))
            .summarize()
            .unwrap();
        assert!(per_file.ends_with("simple_per_file.csv"));
        let file = &read_csv(&per_file)[0];
        assert_eq!(file["ref_species"], "hg18");
        assert_eq!(file["blocks"], "3");
        assert_eq!(file["species"], "5");
        assert_eq!(file["aligned_columns"], "61");
        assert_eq!(file["ref_bases"], "57");
        assert_eq!(file["ref_chrom_size_total"], "158545518");
        assert_eq!(file["min_block_columns"], "6");
        assert_eq!(file["max_block_columns"], "42");
        let species = read_csv(&per_species);
        let names: Vec<&str> = species.iter().map(|r| r["name"].as_str()).collect();
        assert_eq!(names, ["hg18", "baboon", "mm4", "panTro1", "rn3"]);
        let row = |name: &str| species.iter().find(|r| r["name"] == name).unwrap();
        assert_eq!(row("rn3")["blocks_present"], "2");
        assert_eq!(row("rn3")["bases"], "46");
        assert_eq!(row("rn3")["ref_bases_covered"], "43");
        assert_eq!(row("mm4")["ref_bases_covered"], "56");
        assert_eq!(row("hg18")["gaps"], "4");
        assert_eq!(row("hg18")["ref_coverage"], "1.0000");
        assert_eq!(file["ref_genome_coverage"], "3.5952e-7");
    }

    #[test]
    fn test_combined_rows_for_many_files() {
        let output = TempDir::with_prefix("maf-summary-many").unwrap();
        let files = vec![
            PathBuf::from("tests/files/maf/simple.maf"),
            PathBuf::from("tests/files/maf/simple.maf"),
        ];
        let (_, per_species) = MafSummary::new(&files, output.path(), None)
            .with_species(Some(vec!["hg18".to_string(), "rn3".to_string()]))
            .summarize()
            .unwrap();
        let species = read_csv(&per_species);
        let all: Vec<_> = species.iter().filter(|r| r["path"] == ALL_FILES).collect();
        assert_eq!(all.len(), 2);
        assert_eq!(all[1]["name"], "rn3");
        assert_eq!(all[1]["blocks_present"], "4");
    }
}
