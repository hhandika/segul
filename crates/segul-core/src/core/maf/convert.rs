//! Convert MAF blocks into per-locus alignments.
//!
//! Each block is oriented to the forward strand of the reference sequence and its
//! sources are renamed to sample names (`hg18.chr7` becomes `hg18`). Then:
//!
//! - With a BED file, a block is assigned to every region it overlaps on the same
//!   chromosome, clipped to the region, and merged with other blocks of the region
//!   across all input files. Samples missing from a block are padded with gaps.
//! - Without a BED file, each block is written to its own file, or blocks are merged
//!   per reference chromosome.
use std::{
    collections::{BTreeMap, HashMap},
    fmt::Display,
    fs,
    path::{Path, PathBuf},
    str::FromStr,
};

use anyhow::{Context, Result};
use colored::Colorize;
use rayon::prelude::*;

use super::{
    block::{self, RefBlock, SampleNaming, SpeciesResolver},
    clip::clip_columns,
    filter::BlockFilter,
    locus::{LocusPiece, LocusStats, build_locus},
    region::{RegionIndex, sanitize_name},
};
use crate::{
    core::OutputPrint,
    helper::{
        files,
        types::{DnaStrand, Header, OutputFmt, SeqMatrix},
        utils,
    },
    parser::{
        bed::BedParser,
        maf::{MafParagraph, MafReader},
    },
    writer::sequences::SeqWriter,
};

/// Directory for blocks that match no BED region.
pub const UNMATCHED_DIR: &str = "unmatched";
/// Default name of the locus report, without the `.csv` extension.
pub const DEFAULT_REPORT_PREFIX: &str = "maf_convert";
// Blocks written together when each block is its own file.
const WRITE_BATCH: usize = 1024;

/// How blocks are grouped into output files when no BED file is given.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum MafGrouping {
    /// One file per block, named `{chrom}-{start}-{end}`.
    #[default]
    Block,
    /// One file per reference chromosome, with blocks merged in coordinate order.
    Chrom,
}

impl FromStr for MafGrouping {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "block" => Ok(Self::Block),
            "chrom" => Ok(Self::Chrom),
            _ => Err(format!("{} is not a valid MAF grouping", s)),
        }
    }
}

impl Display for MafGrouping {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Block => write!(f, "block"),
            Self::Chrom => write!(f, "chrom"),
        }
    }
}

/// Counts from a MAF conversion.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct MafConvertReport {
    pub files: usize,
    pub blocks_read: usize,
    /// Blocks without a row from the reference species.
    pub blocks_no_ref: usize,
    pub blocks_filtered: usize,
    /// Blocks that overlap no BED region.
    pub blocks_unmatched: usize,
    /// Rows dropped because a sample appeared more than once in a block.
    pub duplicate_rows: usize,
    /// Block pieces dropped because they overlap another piece of the same locus.
    pub overlapping_pieces: usize,
    pub loci_written: usize,
    /// BED regions without any aligned data.
    pub loci_empty: usize,
}

#[derive(Default)]
struct BlockCounts {
    read: usize,
    no_ref: usize,
    filtered: usize,
    duplicates: usize,
}

impl MafConvertReport {
    fn add_counts(&mut self, counts: &BlockCounts) {
        self.blocks_read += counts.read;
        self.blocks_no_ref += counts.no_ref;
        self.blocks_filtered += counts.filtered;
        self.duplicate_rows += counts.duplicates;
    }
}

// Pieces, counts, unmatched count, and unmatched blocks from one file.
type FileRegions<K> = (
    HashMap<K, Vec<LocusPiece>>,
    BlockCounts,
    usize,
    Vec<(usize, usize, RefBlock)>,
);

/// Convert MAF files to FASTA, NEXUS, or PHYLIP alignments.
///
/// # Example
/// ```
/// use std::path::{Path, PathBuf};
/// use segul_core::core::maf::convert::MafConverter;
/// use segul_core::helper::types::OutputFmt;
///
/// let output = tempfile::TempDir::with_prefix("maf").unwrap();
/// let files = vec![PathBuf::from("tests/files/maf/simple.maf")];
/// let report = MafConverter::new(&files, output.path(), &OutputFmt::Fasta)
///     .convert()
///     .unwrap();
/// assert_eq!(report.loci_written, 3);
/// assert!(output.path().join("chr7-27578828-27578866.fas").exists());
/// ```
pub struct MafConverter<'a> {
    input_files: &'a [PathBuf],
    output_dir: &'a Path,
    output_fmt: &'a OutputFmt,
    bed: Option<&'a Path>,
    grouping: MafGrouping,
    clip: bool,
    strand_aware: bool,
    naming: SampleNaming,
    ref_species: Option<String>,
    filter: BlockFilter,
    write_unmatched: bool,
    prefix: Option<String>,
}

impl OutputPrint for MafConverter<'_> {}

impl<'a> MafConverter<'a> {
    pub fn new(
        input_files: &'a [PathBuf],
        output_dir: &'a Path,
        output_fmt: &'a OutputFmt,
    ) -> Self {
        Self {
            input_files,
            output_dir,
            output_fmt,
            bed: None,
            grouping: MafGrouping::default(),
            clip: true,
            strand_aware: false,
            naming: SampleNaming::default(),
            ref_species: None,
            filter: BlockFilter::default(),
            write_unmatched: false,
            prefix: None,
        }
    }

    /// Assign blocks to the regions in a BED file.
    pub fn with_bed(mut self, bed: &'a Path) -> Self {
        self.bed = Some(bed);
        self
    }

    /// Group blocks without a BED file. Ignored when a BED file is set.
    pub fn with_grouping(mut self, grouping: MafGrouping) -> Self {
        self.grouping = grouping;
        self
    }

    /// Clip blocks to BED regions (default). When off, overlapping blocks are kept whole.
    pub fn with_clip(mut self, clip: bool) -> Self {
        self.clip = clip;
        self
    }

    /// Reverse complement loci whose BED strand is `-`.
    pub fn with_strand_aware(mut self, strand_aware: bool) -> Self {
        self.strand_aware = strand_aware;
        self
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

    pub fn with_filter(mut self, filter: BlockFilter) -> Self {
        self.filter = filter;
        self
    }

    /// Write blocks that match no BED region to the `unmatched` directory.
    pub fn with_unmatched(mut self, write_unmatched: bool) -> Self {
        self.write_unmatched = write_unmatched;
        self
    }

    /// Name of the locus report written with a BED file, without the extension.
    pub fn with_prefix(mut self, prefix: Option<String>) -> Self {
        self.prefix = prefix;
        self
    }

    pub fn convert(&self) -> Result<MafConvertReport> {
        let mut input_files = self.input_files.to_vec();
        input_files.sort();
        let resolver = SpeciesResolver::new(self.naming, self.filter.species.as_deref());
        let ref_sample = match &self.ref_species {
            Some(species) => species.clone(),
            None => block::first_sample(&input_files, &resolver)?,
        };
        fs::create_dir_all(self.output_dir).with_context(|| {
            format!(
                "Failed creating output directory {}",
                self.output_dir.display()
            )
        })?;
        let spin = utils::set_spinner();
        spin.set_message("Converting MAF files...");
        let mut report = match (self.bed, self.grouping) {
            (Some(bed), _) => self.convert_regions(bed, &input_files, &resolver, &ref_sample)?,
            (None, MafGrouping::Block) => {
                self.convert_blocks(&input_files, &resolver, &ref_sample)?
            }
            (None, MafGrouping::Chrom) => {
                self.convert_chroms(&input_files, &resolver, &ref_sample)?
            }
        };
        report.files = input_files.len();
        spin.finish_with_message("Finished converting MAF files!\n");
        self.print_output_info(&ref_sample, &report);
        Ok(report)
    }

    fn convert_regions(
        &self,
        bed: &Path,
        input_files: &[PathBuf],
        resolver: &SpeciesResolver,
        ref_sample: &str,
    ) -> Result<MafConvertReport> {
        let index = RegionIndex::from_bed(&BedParser::new(bed).parse()?);
        let results = input_files
            .par_iter()
            .enumerate()
            .map(|(file_idx, file)| {
                self.collect_region_pieces(file_idx, file, &index, resolver, ref_sample)
            })
            .collect::<Result<Vec<FileRegions<usize>>>>()?;
        let mut report = MafConvertReport::default();
        let mut by_region: BTreeMap<usize, Vec<LocusPiece>> = BTreeMap::new();
        let mut unmatched = Vec::new();
        for (pieces, counts, unmatched_count, blocks) in results {
            report.add_counts(&counts);
            report.blocks_unmatched += unmatched_count;
            unmatched.extend(blocks);
            for (id, region_pieces) in pieces {
                by_region.entry(id).or_default().extend(region_pieces);
            }
        }
        let regions = index.regions();
        let written = by_region
            .into_par_iter()
            .map(|(id, pieces)| {
                let region = &regions[id];
                let reverse = self.strand_aware && region.strand == DnaStrand::Reverse;
                let (matrix, stats) = build_locus(pieces, ref_sample, reverse);
                let is_written = match matrix {
                    Some(matrix) => {
                        self.write_matrix(self.output_dir, &region.label, &matrix)?;
                        true
                    }
                    None => false,
                };
                Ok((id, (stats, is_written)))
            })
            .collect::<Result<HashMap<usize, (LocusStats, bool)>>>()?;
        report.overlapping_pieces = written.values().map(|(s, _)| s.overlapping).sum();
        report.loci_written = written.values().filter(|(_, w)| *w).count();
        report.loci_empty = regions.len() - report.loci_written;
        self.write_report(&index, &written)?;
        if self.write_unmatched && !unmatched.is_empty() {
            unmatched.sort_by_key(|(file_idx, block_idx, _)| (*file_idx, *block_idx));
            let blocks: Vec<RefBlock> = unmatched.into_iter().map(|(_, _, b)| b).collect();
            self.write_named_blocks(&self.output_dir.join(UNMATCHED_DIR), blocks)?;
        }
        Ok(report)
    }

    fn collect_region_pieces(
        &self,
        file_idx: usize,
        file: &Path,
        index: &RegionIndex,
        resolver: &SpeciesResolver,
        ref_sample: &str,
    ) -> Result<FileRegions<usize>> {
        let mut pieces: HashMap<usize, Vec<LocusPiece>> = HashMap::new();
        let mut unmatched_count = 0;
        let mut unmatched = Vec::new();
        let counts = self.read_blocks(file, resolver, ref_sample, |block_idx, block| {
            let hits = index.overlapping(&block.chrom, &block.source, block.start, block.end);
            if hits.is_empty() {
                unmatched_count += 1;
                if self.write_unmatched {
                    unmatched.push((file_idx, block_idx, block));
                }
                return Ok(());
            }
            for id in hits {
                let region = &index.regions()[id];
                let piece = if self.clip {
                    let Some(range) =
                        clip_columns(&block.ref_text, block.start, region.start, region.end)
                    else {
                        continue;
                    };
                    LocusPiece {
                        ref_start: block.start.max(region.start),
                        ref_end: block.end.min(region.end),
                        file_idx,
                        block_idx,
                        width: range.len(),
                        rows: block
                            .rows
                            .iter()
                            .map(|(sample, row)| (sample.clone(), row[range.clone()].to_vec()))
                            .collect(),
                    }
                } else {
                    Self::whole_piece(&block, file_idx, block_idx)
                };
                pieces.entry(id).or_default().push(piece);
            }
            Ok(())
        })?;
        Ok((pieces, counts, unmatched_count, unmatched))
    }

    fn convert_chroms(
        &self,
        input_files: &[PathBuf],
        resolver: &SpeciesResolver,
        ref_sample: &str,
    ) -> Result<MafConvertReport> {
        let results = input_files
            .par_iter()
            .enumerate()
            .map(|(file_idx, file)| {
                let mut pieces: HashMap<String, Vec<LocusPiece>> = HashMap::new();
                let counts = self.read_blocks(file, resolver, ref_sample, |block_idx, block| {
                    let piece = Self::whole_piece(&block, file_idx, block_idx);
                    pieces.entry(block.chrom).or_default().push(piece);
                    Ok(())
                })?;
                Ok((pieces, counts, 0, Vec::new()))
            })
            .collect::<Result<Vec<FileRegions<String>>>>()?;
        let mut report = MafConvertReport::default();
        let mut by_chrom: BTreeMap<String, Vec<LocusPiece>> = BTreeMap::new();
        for (pieces, counts, _, _) in results {
            report.add_counts(&counts);
            for (chrom, chrom_pieces) in pieces {
                by_chrom.entry(chrom).or_default().extend(chrom_pieces);
            }
        }
        let stats = by_chrom
            .into_par_iter()
            .map(|(chrom, pieces)| {
                let (matrix, stats) = build_locus(pieces, ref_sample, false);
                if let Some(matrix) = matrix {
                    self.write_matrix(self.output_dir, &sanitize_name(&chrom), &matrix)?;
                }
                Ok(stats)
            })
            .collect::<Result<Vec<LocusStats>>>()?;
        report.overlapping_pieces = stats.iter().map(|s| s.overlapping).sum();
        report.loci_written = stats
            .iter()
            .filter(|s| s.columns > 0 && s.samples > 0)
            .count();
        Ok(report)
    }

    fn convert_blocks(
        &self,
        input_files: &[PathBuf],
        resolver: &SpeciesResolver,
        ref_sample: &str,
    ) -> Result<MafConvertReport> {
        let mut report = MafConvertReport::default();
        let mut names: HashMap<String, usize> = HashMap::new();
        let mut batch: Vec<(String, RefBlock)> = Vec::new();
        // Files are read in order so repeated names get the same suffix on every run.
        for file in input_files {
            let counts = self.read_blocks(file, resolver, ref_sample, |_, block| {
                let name = unique_name(&mut names, block_name(&block));
                batch.push((name, block));
                if batch.len() >= WRITE_BATCH {
                    report.loci_written += self.write_batch(self.output_dir, &mut batch)?;
                }
                Ok(())
            })?;
            report.add_counts(&counts);
        }
        report.loci_written += self.write_batch(self.output_dir, &mut batch)?;
        Ok(report)
    }

    /// Stream the alignment blocks of a file, normalized and filtered.
    fn read_blocks<F>(
        &self,
        file: &Path,
        resolver: &SpeciesResolver,
        ref_sample: &str,
        mut on_block: F,
    ) -> Result<BlockCounts>
    where
        F: FnMut(usize, RefBlock) -> Result<()>,
    {
        let mut counts = BlockCounts::default();
        let reader = MafReader::new(files::open_reader(file)?);
        for paragraph in reader {
            let paragraph =
                paragraph.with_context(|| format!("Failed reading {}", file.display()))?;
            let MafParagraph::Alignment(alignment) = paragraph else {
                continue;
            };
            let block_idx = counts.read;
            counts.read += 1;
            let block =
                RefBlock::from_alignment(alignment, resolver, ref_sample).with_context(|| {
                    format!("Invalid block {} in {}", block_idx + 1, file.display())
                })?;
            match block {
                None => counts.no_ref += 1,
                Some(block) => {
                    counts.duplicates += block.duplicates;
                    if self.filter.passes(&block, ref_sample) {
                        on_block(block_idx, block)?;
                    } else {
                        counts.filtered += 1;
                    }
                }
            }
        }
        Ok(counts)
    }

    fn whole_piece(block: &RefBlock, file_idx: usize, block_idx: usize) -> LocusPiece {
        LocusPiece {
            ref_start: block.start,
            ref_end: block.end,
            file_idx,
            block_idx,
            width: block.ncols(),
            rows: block.rows.clone(),
        }
    }

    fn write_named_blocks(&self, dir: &Path, blocks: Vec<RefBlock>) -> Result<usize> {
        let mut names: HashMap<String, usize> = HashMap::new();
        let mut batch: Vec<(String, RefBlock)> = blocks
            .into_iter()
            .map(|block| (unique_name(&mut names, block_name(&block)), block))
            .collect();
        self.write_batch(dir, &mut batch)
    }

    /// Write and clear a batch of blocks, one file each. Returns the number written.
    fn write_batch(&self, dir: &Path, batch: &mut Vec<(String, RefBlock)>) -> Result<usize> {
        let written = batch
            .par_iter()
            .map(|(name, block)| {
                if block.rows.is_empty() {
                    return Ok(0);
                }
                let matrix: SeqMatrix = block
                    .rows
                    .iter()
                    .map(|(sample, row)| {
                        (sample.clone(), String::from_utf8_lossy(row).into_owned())
                    })
                    .collect();
                self.write_matrix(dir, name, &matrix)?;
                Ok(1)
            })
            .collect::<Result<Vec<usize>>>()?
            .into_iter()
            .sum();
        batch.clear();
        Ok(written)
    }

    fn write_matrix(&self, dir: &Path, stem: &str, matrix: &SeqMatrix) -> Result<()> {
        let output = files::create_output_fname_from_stem(dir, stem, self.output_fmt);
        let mut header = Header::new();
        header.from_seq_matrix(matrix, true);
        SeqWriter::new(&output, matrix, &header)
            .write_sequence(self.output_fmt)
            .with_context(|| format!("Failed writing {}", output.display()))
    }

    fn write_report(
        &self,
        index: &RegionIndex,
        loci: &HashMap<usize, (LocusStats, bool)>,
    ) -> Result<()> {
        let prefix = self.prefix.as_deref().unwrap_or(DEFAULT_REPORT_PREFIX);
        let path = self.output_dir.join(format!("{}.csv", prefix));
        let mut writer = csv::Writer::from_path(&path)
            .with_context(|| format!("Failed creating {}", path.display()))?;
        writer.write_record([
            "locus",
            "chrom",
            "start",
            "end",
            "strand",
            "blocks",
            "ref_bases",
            "aligned_columns",
            "samples",
            "written",
        ])?;
        let empty = (LocusStats::default(), false);
        for (id, region) in index.regions().iter().enumerate() {
            let (stats, written) = loci.get(&id).unwrap_or(&empty);
            writer.write_record([
                region.label.clone(),
                region.chrom.clone(),
                region.start.to_string(),
                region.end.to_string(),
                region.strand.to_string(),
                stats.blocks.to_string(),
                stats.ref_bases.to_string(),
                stats.columns.to_string(),
                stats.samples.to_string(),
                written.to_string(),
            ])?;
        }
        writer.flush()?;
        Ok(())
    }

    fn print_output_info(&self, ref_sample: &str, report: &MafConvertReport) {
        log::info!("{}", "Blocks".yellow());
        log::info!("{:18}: {}", "Reference", ref_sample);
        log::info!("{:18}: {}", "Read", utils::fmt_num(&report.blocks_read));
        log::info!(
            "{:18}: {}",
            "No reference",
            utils::fmt_num(&report.blocks_no_ref)
        );
        log::info!(
            "{:18}: {}",
            "Filtered",
            utils::fmt_num(&report.blocks_filtered)
        );
        if self.bed.is_some() {
            log::info!(
                "{:18}: {}",
                "Unmatched",
                utils::fmt_num(&report.blocks_unmatched)
            );
        }
        if report.duplicate_rows > 0 {
            log::warn!("{:18}: {}", "Duplicate rows", report.duplicate_rows);
        }
        if report.overlapping_pieces > 0 {
            log::warn!("{:18}: {}", "Overlaps dropped", report.overlapping_pieces);
        }
        log::info!("");
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Output dir", self.output_dir.display());
        log::info!(
            "{:18}: {}",
            "Files written",
            utils::fmt_num(&report.loci_written)
        );
        if self.bed.is_some() {
            log::info!(
                "{:18}: {}",
                "Empty loci",
                utils::fmt_num(&report.loci_empty)
            );
        }
        self.print_output_fmt(self.output_fmt);
    }
}

fn block_name(block: &RefBlock) -> String {
    sanitize_name(&format!("{}-{}-{}", block.chrom, block.start, block.end))
}

/// Return `base` the first time, then `base-2`, `base-3`, and so on.
fn unique_name(names: &mut HashMap<String, usize>, base: String) -> String {
    let count = names.entry(base.clone()).or_insert(0);
    *count += 1;
    if *count == 1 {
        base
    } else {
        format!("{}-{}", base, count)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    use flate2::{Compression, write::GzEncoder};
    use tempfile::TempDir;

    const SIMPLE_MAF: &str = "tests/files/maf/simple.maf";
    const SIMPLE_REGIONS: &str = "tests/files/bed/simple_regions.bed";

    fn read_fasta(path: &Path) -> Vec<(String, String)> {
        let text = fs::read_to_string(path).unwrap();
        let mut records = Vec::new();
        for entry in text.split('>').skip(1) {
            let mut lines = entry.lines();
            let id = lines.next().unwrap().trim().to_string();
            records.push((id, lines.collect::<String>()));
        }
        records
    }

    fn output_names(dir: &Path) -> Vec<String> {
        let mut names: Vec<String> = fs::read_dir(dir)
            .unwrap()
            .filter_map(|e| e.ok())
            .filter(|e| e.path().is_file())
            .map(|e| e.file_name().to_string_lossy().to_string())
            .collect();
        names.sort();
        names
    }

    fn sequence_of(records: &[(String, String)], id: &str) -> String {
        records
            .iter()
            .find(|(name, _)| name == id)
            .unwrap()
            .1
            .clone()
    }

    #[test]
    fn test_blocks_without_bed() {
        let output = TempDir::with_prefix("maf-blocks").unwrap();
        let files = vec![PathBuf::from(SIMPLE_MAF)];
        let report = MafConverter::new(&files, output.path(), &OutputFmt::Fasta)
            .convert()
            .unwrap();
        assert_eq!(report.blocks_read, 3);
        assert_eq!(
            output_names(output.path()),
            [
                "chr7-27578828-27578866.fas",
                "chr7-27699739-27699745.fas",
                "chr7-27707221-27707234.fas"
            ]
        );
        let records = read_fasta(&output.path().join("chr7-27578828-27578866.fas"));
        let ids: Vec<&str> = records.iter().map(|(id, _)| id.as_str()).collect();
        assert_eq!(ids, ["hg18", "panTro1", "baboon", "mm4", "rn3"]);
    }

    #[test]
    fn test_repeated_blocks_get_suffixes() {
        let output = TempDir::with_prefix("maf-twice").unwrap();
        let files = vec![PathBuf::from(SIMPLE_MAF), PathBuf::from(SIMPLE_MAF)];
        let report = MafConverter::new(&files, output.path(), &OutputFmt::Fasta)
            .convert()
            .unwrap();
        assert_eq!(report.loci_written, 6);
        assert!(output.path().join("chr7-27578828-27578866-2.fas").exists());
    }

    #[test]
    fn test_group_by_chrom() {
        let output = TempDir::with_prefix("maf-chrom").unwrap();
        let files = vec![PathBuf::from(SIMPLE_MAF)];
        MafConverter::new(&files, output.path(), &OutputFmt::Fasta)
            .with_grouping(MafGrouping::Chrom)
            .convert()
            .unwrap();
        let records = read_fasta(&output.path().join("chr7.fas"));
        assert!(records.iter().all(|(_, seq)| seq.len() == 61));
        assert!(sequence_of(&records, "rn3").ends_with(&"-".repeat(13)));
    }

    #[test]
    fn test_bed_regions() {
        let output = TempDir::with_prefix("maf-bed").unwrap();
        let files = vec![PathBuf::from(SIMPLE_MAF)];
        let bed = PathBuf::from(SIMPLE_REGIONS);
        let report = MafConverter::new(&files, output.path(), &OutputFmt::Fasta)
            .with_bed(&bed)
            .convert()
            .unwrap();
        assert_eq!(report.loci_written, 4);
        assert_eq!(report.loci_empty, 2);
        assert_eq!(
            output_names(output.path()),
            [
                "clipA-chr7-27578830-27578834.fas",
                "gapIns-chr7-27578848-27578852.fas",
                "locA-chr7-27578828-27578866.fas",
                "locMerge-chr7-27699739-27707234.fas",
                "maf_convert.csv",
            ]
        );
        let clip = read_fasta(&output.path().join("clipA-chr7-27578830-27578834.fas"));
        assert_eq!(sequence_of(&clip, "hg18"), "A-GGG");
        assert_eq!(sequence_of(&clip, "mm4"), "ATGGG");
        assert_eq!(sequence_of(&clip, "rn3"), "A-GGG");
        let insertion = read_fasta(&output.path().join("gapIns-chr7-27578848-27578852.fas"));
        assert_eq!(sequence_of(&insertion, "hg18"), "GA---AT");
        assert_eq!(sequence_of(&insertion, "rn3"), "GAGTTGT");
        let merged = read_fasta(&output.path().join("locMerge-chr7-27699739-27707234.fas"));
        assert_eq!(sequence_of(&merged, "hg18"), "TAAAGAgcagctgaaaaca");
        assert_eq!(
            sequence_of(&merged, "rn3"),
            format!("taagga{}", "-".repeat(13))
        );
        let report_csv = fs::read_to_string(output.path().join("maf_convert.csv")).unwrap();
        assert!(report_csv.starts_with("locus,chrom,start,end,strand,blocks"));
        assert!(report_csv.contains("noData-chr7-1-100,chr7,1,100,+,0,0,0,0,false"));
        assert!(report_csv.contains("otherChr-chr8-27578828-27578866"));
    }

    #[test]
    fn test_strand_aware_and_no_clip() {
        let output = TempDir::with_prefix("maf-strand").unwrap();
        let files = vec![PathBuf::from(SIMPLE_MAF)];
        let bed = PathBuf::from(SIMPLE_REGIONS);
        MafConverter::new(&files, output.path(), &OutputFmt::Fasta)
            .with_bed(&bed)
            .with_strand_aware(true)
            .with_clip(false)
            .convert()
            .unwrap();
        let merged = read_fasta(&output.path().join("locMerge-chr7-27699739-27707234.fas"));
        assert_eq!(sequence_of(&merged, "hg18"), "tgttttcagctgcTCTTTA");
        let clip = read_fasta(&output.path().join("clipA-chr7-27578830-27578834.fas"));
        assert_eq!(sequence_of(&clip, "hg18").len(), 42);
    }

    #[test]
    fn test_filters_and_unmatched() {
        let output = TempDir::with_prefix("maf-filter").unwrap();
        let files = vec![PathBuf::from(SIMPLE_MAF)];
        let filter = BlockFilter {
            min_ref_len: Some(10),
            species: Some(vec!["hg18".to_string(), "mm4".to_string()]),
            ..Default::default()
        };
        let report = MafConverter::new(&files, output.path(), &OutputFmt::Fasta)
            .with_filter(filter)
            .convert()
            .unwrap();
        assert_eq!(report.blocks_filtered, 1);
        let records = read_fasta(&output.path().join("chr7-27707221-27707234.fas"));
        let ids: Vec<&str> = records.iter().map(|(id, _)| id.as_str()).collect();
        assert_eq!(ids, ["hg18", "mm4"]);

        let output = TempDir::with_prefix("maf-unmatched").unwrap();
        let bed = output.path().join("other.bed");
        fs::write(&bed, "chr1\t0\t10\tx\n").unwrap();
        let out_dir = output.path().join("out");
        let report = MafConverter::new(&files, &out_dir, &OutputFmt::Fasta)
            .with_bed(&bed)
            .with_unmatched(true)
            .convert()
            .unwrap();
        assert_eq!(report.blocks_unmatched, 3);
        assert_eq!(output_names(&out_dir.join(UNMATCHED_DIR)).len(), 3);
    }

    #[test]
    fn test_file_order_and_gzip_do_not_change_output() {
        let input = TempDir::with_prefix("maf-order-input").unwrap();
        let text = fs::read_to_string(SIMPLE_MAF).unwrap();
        let (first, rest) = text.split_at(text.find("a score=5062.0").unwrap());
        let file_a = input.path().join("a.maf");
        fs::write(&file_a, first).unwrap();
        let file_b = input.path().join("b.maf.gz");
        let mut encoder = GzEncoder::new(fs::File::create(&file_b).unwrap(), Compression::fast());
        encoder.write_all(rest.as_bytes()).unwrap();
        encoder.finish().unwrap();
        let bed = PathBuf::from(SIMPLE_REGIONS);
        let run = |files: Vec<PathBuf>, name: &str| {
            let dir = input.path().join(name);
            MafConverter::new(&files, &dir, &OutputFmt::Fasta)
                .with_bed(&bed)
                .with_ref_species(Some("hg18".to_string()))
                .convert()
                .unwrap();
            fs::read_to_string(dir.join("locMerge-chr7-27699739-27707234.fas")).unwrap()
        };
        let forward = run(vec![file_a.clone(), file_b.clone()], "forward");
        let backward = run(vec![file_b, file_a], "backward");
        assert_eq!(forward, backward);
        let whole = TempDir::with_prefix("maf-order-whole").unwrap();
        let files = vec![PathBuf::from(SIMPLE_MAF)];
        MafConverter::new(&files, whole.path(), &OutputFmt::Fasta)
            .with_bed(&bed)
            .convert()
            .unwrap();
        let expected =
            fs::read_to_string(whole.path().join("locMerge-chr7-27699739-27707234.fas")).unwrap();
        assert_eq!(forward, expected);
    }

    #[test]
    fn test_nexus_output() {
        let output = TempDir::with_prefix("maf-nexus").unwrap();
        let files = vec![PathBuf::from(SIMPLE_MAF)];
        MafConverter::new(&files, output.path(), &OutputFmt::Nexus)
            .convert()
            .unwrap();
        let text = fs::read_to_string(output.path().join("chr7-27578828-27578866.nex")).unwrap();
        assert!(text.to_uppercase().starts_with("#NEXUS"));
    }
}
