//! Read alignments for interactive viewers.
//!
//! Residues are read on request in small tiles, so viewers only touch the
//! part of the alignment in view. FASTA files with regular line widths are
//! memory-mapped and indexed by byte offset, like `samtools faidx`. Other
//! files are parsed by SEGUL and kept in memory.
//!
//! # Example
//!
//! ```
//! use std::path::Path;
//! use segul_core::core::align::view::AlignmentView;
//! use segul_core::helper::types::{DataType, InputFmt};
//!
//! let path = Path::new("tests/files/simple.fas");
//! let view = AlignmentView::open(path, &InputFmt::Auto, &DataType::Ignore).unwrap();
//! let summary = view.summary();
//! assert_eq!(summary.names, vec!["ABCD", "CDEF"]);
//! assert!(summary.is_dna);
//! // The first three sites of both sequences, row after row.
//! assert_eq!(view.tile(0, 2, 0, 3), b"AGTATG".to_vec());
//! ```

use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::File;
use std::ops::Range;
use std::panic::{self, AssertUnwindSafe};
use std::path::Path;
use std::sync::{Arc, Mutex, RwLock};

use memchr::memchr;
use memmap2::Mmap;
use rayon::prelude::*;

use crate::helper::sequence::SeqParser;
use crate::helper::types::{DataType, InputFmt};

const DNA_CHARS: &[u8] = b"ACGTUNRYSWKMBDHV";
// Share of residues that must be nucleotide codes to call a file DNA.
const DNA_THRESHOLD: f64 = 0.9;
// Residues read from the start of each row to guess the data type.
const DNA_SAMPLE_PER_ROW: usize = 4096;
const DNA_SAMPLE_TOTAL: usize = 1 << 20;
// Sites between stored ungapped counts. Lookups scan at most this many.
const UNGAPPED_STEP: usize = 4096;
// Sites counted together for chunk statistics. Blocks are counted in
// parallel, only when a chunk over them is asked for, and then kept.
const STAT_BLOCK: usize = 4096;
const GAP: u8 = b'-';
const MISSING: u8 = b'?';
const DNA_UNKNOWN: u8 = b'N';
const AA_UNKNOWN: u8 = b'X';

/// Names and shape of an open alignment.
#[derive(Debug, Clone, PartialEq)]
pub struct AlignmentSummary {
    pub names: Vec<String>,
    /// Residues in each row.
    pub lengths: Vec<usize>,
    /// Longest sequence length.
    pub length: usize,
    /// Whether every sequence has the same length.
    pub is_aligned: bool,
    /// Whether the sequences are DNA rather than amino acids.
    pub is_dna: bool,
}

/// Gap and variable-site share of a run of sites.
#[derive(Debug, Clone, PartialEq)]
pub struct ChunkStat {
    /// Zero-based first site.
    pub start: usize,
    /// Exclusive last site.
    pub end: usize,
    /// Share of cells that are gaps, missing data, or past a sequence's end.
    pub gap_fraction: f64,
    /// Share of sites with more than one residue, ignoring gaps and unknowns.
    pub variable_fraction: f64,
}

/// An open sequence file. Residues are read on request.
pub struct AlignmentView {
    names: Vec<String>,
    lengths: Vec<usize>,
    length: usize,
    is_aligned: bool,
    is_dna: bool,
    store: Store,
    // Per-site counts of each STAT_BLOCK sites counted so far.
    site_blocks: RwLock<HashMap<usize, Arc<SiteBlock>>>,
    // Residues before every UNGAPPED_STEP sites, per row, built on first use.
    ungapped: Mutex<HashMap<usize, Vec<usize>>>,
}

impl AlignmentView {
    /// Opens a FASTA, NEXUS, or PHYLIP file. `InputFmt::Auto` infers the
    /// format from the extension. `DataType::Ignore` guesses DNA or amino
    /// acids from the residues; `Dna` and `Aa` set it.
    pub fn open(path: &Path, input_fmt: &InputFmt, datatype: &DataType) -> anyhow::Result<Self> {
        let input_fmt = match input_fmt {
            InputFmt::Auto => Self::format_of(path)?,
            other => *other,
        };
        let indexed = match input_fmt {
            InputFmt::Fasta => FastaIndex::build(path)?,
            _ => None,
        };
        let (names, store) = match indexed {
            Some(index) => (index.names, Store::Mapped(index.map, index.records)),
            None => Self::parse(path, &input_fmt)?,
        };
        let mut view = Self::new(names, store);
        match datatype {
            DataType::Dna => view.is_dna = true,
            DataType::Aa => view.is_dna = false,
            DataType::Ignore => (),
        }
        Ok(view)
    }

    /// Names and shape of the alignment.
    pub fn summary(&self) -> AlignmentSummary {
        AlignmentSummary {
            names: self.names.clone(),
            lengths: self.lengths.clone(),
            length: self.length,
            is_aligned: self.is_aligned,
            is_dna: self.is_dna,
        }
    }

    /// Number of sequences.
    pub fn rows(&self) -> usize {
        self.names.len()
    }

    /// Longest sequence length.
    pub fn length(&self) -> usize {
        self.length
    }

    /// Residues in `row`.
    pub fn row_len(&self, row: usize) -> usize {
        self.lengths[row]
    }

    /// Name of the sequence at `row`.
    pub fn name(&self, row: usize) -> &str {
        &self.names[row]
    }

    /// All sequence names, in file order.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Whether every sequence has the same length.
    pub fn is_aligned(&self) -> bool {
        self.is_aligned
    }

    /// Whether the sequences are DNA rather than amino acids.
    pub fn is_dna(&self) -> bool {
        self.is_dna
    }

    /// Upper-case residue at `site` of `row`, or none past the row's end.
    pub fn residue(&self, row: usize, site: usize) -> Option<u8> {
        if row < self.rows() && site < self.lengths[row] {
            Some(self.store.residue(row, site))
        } else {
            None
        }
    }

    /// Residues of `row_count` rows from `row_start`, `site_count` sites each
    /// from `site_start`, row after row. Cells past a row's end are zero.
    pub fn tile(
        &self,
        row_start: usize,
        row_count: usize,
        site_start: usize,
        site_count: usize,
    ) -> Vec<u8> {
        let rows = row_start.min(self.rows())..(row_start + row_count).min(self.rows());
        let mut out = vec![0u8; rows.len() * site_count];
        for (row, cells) in rows.zip(out.chunks_exact_mut(site_count.max(1))) {
            self.store.read(row, site_start, cells);
        }
        out
    }

    /// Residues of `rows` at `width` evenly spaced sites, row after row, for
    /// an overview. Cells past a row's end are zero.
    pub fn sample(&self, rows: &[usize], width: usize) -> Vec<u8> {
        let mut out = vec![0u8; rows.len() * width];
        for (&row, cells) in rows.iter().zip(out.chunks_exact_mut(width.max(1))) {
            if row >= self.rows() {
                continue;
            }
            for (x, cell) in cells.iter_mut().enumerate() {
                let site = x * self.length / width;
                if site < self.lengths[row] {
                    *cell = self.store.residue(row, site);
                }
            }
        }
        out
    }

    /// Gap and variable-site shares of up to `count` chunks of `chunk_size`
    /// sites, from chunk `first`. Only the sites of those chunks are read;
    /// their counts are kept, so later calls and other chunk sizes reuse them.
    pub fn chunk_stats(&self, chunk_size: usize, first: usize, count: usize) -> Vec<ChunkStat> {
        let size = chunk_size.max(1);
        let start = first.saturating_mul(size).min(self.length);
        let end = first
            .saturating_add(count)
            .saturating_mul(size)
            .min(self.length);
        if start >= end {
            return Vec::new();
        }
        let blocks = self.site_blocks(start..end);
        let first_block = start / STAT_BLOCK;
        let rows = self.rows();
        (start..end)
            .step_by(size)
            .map(|chunk_start| {
                let chunk_end = (chunk_start + size).min(end);
                let mut gaps = 0u64;
                let mut variable = 0usize;
                let mut site = chunk_start;
                while site < chunk_end {
                    let block = &blocks[site / STAT_BLOCK - first_block];
                    let block_start = site / STAT_BLOCK * STAT_BLOCK;
                    let to = chunk_end.min(block_start + STAT_BLOCK);
                    let range = site - block_start..to - block_start;
                    gaps += block.gaps[range.clone()]
                        .iter()
                        .map(|&g| g as u64)
                        .sum::<u64>();
                    variable += block.variable[range].iter().filter(|&&v| v).count();
                    site = to;
                }
                let sites = chunk_end - chunk_start;
                let cells = sites * rows;
                ChunkStat {
                    start: chunk_start,
                    end: chunk_end,
                    gap_fraction: if cells == 0 {
                        0.0
                    } else {
                        gaps as f64 / cells as f64
                    },
                    variable_fraction: variable as f64 / sites as f64,
                }
            })
            .collect()
    }

    /// One-based position of `site` in `row` without gaps, or none on a gap.
    pub fn ungapped_position(&self, row: usize, site: usize) -> Option<usize> {
        if row >= self.rows() || site >= self.lengths[row] {
            return None;
        }
        if is_gap(self.store.residue(row, site)) {
            return None;
        }
        let step = site / UNGAPPED_STEP;
        let before = {
            let mut cache = self.ungapped.lock().ok()?;
            let counts = cache
                .entry(row)
                .or_insert_with(|| self.ungapped_counts(row));
            counts[step]
        };
        let start = step * UNGAPPED_STEP;
        let mut cells = vec![0u8; site + 1 - start];
        self.store.read(row, start, &mut cells);
        let residues = cells.iter().filter(|&&c| !is_gap(c)).count();
        Some(before + residues)
    }

    /// Rows in FASTA format, in alignment order.
    pub fn fasta(&self, rows: &[usize]) -> String {
        let mut rows: Vec<usize> = rows.iter().copied().filter(|&r| r < self.rows()).collect();
        rows.sort_unstable();
        rows.dedup();
        let mut out = String::new();
        for row in rows {
            let mut cells = vec![0u8; self.lengths[row]];
            self.store.read(row, 0, &mut cells);
            out.push('>');
            out.push_str(&self.names[row]);
            out.push('\n');
            out.push_str(&String::from_utf8_lossy(&cells));
            out.push('\n');
        }
        out
    }

    fn new(names: Vec<String>, store: Store) -> Self {
        let lengths: Vec<usize> = (0..names.len()).map(|row| store.len(row)).collect();
        let length = lengths.iter().copied().max().unwrap_or(0);
        let is_aligned = lengths.iter().all(|&l| l == length);
        let is_dna = Self::guess_dna(&store, &lengths);
        Self {
            names,
            lengths,
            length,
            is_aligned,
            is_dna,
            store,
            site_blocks: RwLock::new(HashMap::new()),
            ungapped: Mutex::new(HashMap::new()),
        }
    }

    fn format_of(path: &Path) -> anyhow::Result<InputFmt> {
        let ext = path
            .extension()
            .and_then(OsStr::to_str)
            .map(str::to_ascii_lowercase)
            .unwrap_or_default();
        match ext.as_str() {
            "fa" | "fasta" | "fna" | "fsa" | "fas" => Ok(InputFmt::Fasta),
            "nex" | "nxs" | "nexus" => Ok(InputFmt::Nexus),
            "phy" | "phylip" | "ph" => Ok(InputFmt::Phylip),
            _ => anyhow::bail!(
                "Failed to infer the format of {}. Use a FASTA, NEXUS, or PHYLIP extension",
                path.display()
            ),
        }
    }

    fn parse(path: &Path, input_fmt: &InputFmt) -> anyhow::Result<(Vec<String>, Store)> {
        if !path.is_file() {
            anyhow::bail!("Failed to open {}: no such file", path.display());
        }
        // SEGUL parsers panic on malformed input; report it as an error instead.
        let (matrix, _) = panic::catch_unwind(AssertUnwindSafe(|| {
            SeqParser::new(path, &DataType::Ignore).parse(input_fmt)
        }))
        .map_err(|_| anyhow::anyhow!("Failed to read {} as a sequence file", path.display()))?;
        let mut names = Vec::with_capacity(matrix.len());
        let mut rows = Vec::with_capacity(matrix.len());
        for (name, seq) in matrix {
            let mut seq = seq.into_bytes();
            seq.make_ascii_uppercase();
            names.push(name);
            rows.push(seq);
        }
        Ok((names, Store::Memory(rows)))
    }

    // Guesses the data type from the start of each row.
    fn guess_dna(store: &Store, lengths: &[usize]) -> bool {
        let per_row = (DNA_SAMPLE_TOTAL / lengths.len().max(1)).clamp(1, DNA_SAMPLE_PER_ROW);
        let mut residues = 0usize;
        let mut nucleotides = 0usize;
        let mut cells = Vec::new();
        for (row, &len) in lengths.iter().enumerate() {
            cells.resize(len.min(per_row), 0);
            store.read(row, 0, &mut cells);
            for &c in &cells {
                if c == GAP || c == MISSING || c == b'.' {
                    continue;
                }
                residues += 1;
                if DNA_CHARS.contains(&c) {
                    nucleotides += 1;
                }
            }
        }
        residues == 0 || nucleotides as f64 / residues as f64 >= DNA_THRESHOLD
    }

    // Counts of every block over `sites`, in order, counting missing ones in
    // parallel.
    fn site_blocks(&self, sites: Range<usize>) -> Vec<Arc<SiteBlock>> {
        let ids = sites.start / STAT_BLOCK..sites.end.div_ceil(STAT_BLOCK);
        let missing: Vec<usize> = match self.site_blocks.read() {
            Ok(cache) => ids.clone().filter(|id| !cache.contains_key(id)).collect(),
            Err(_) => ids.clone().collect(),
        };
        let counted: Vec<(usize, Arc<SiteBlock>)> = missing
            .into_par_iter()
            .map(|id| (id, Arc::new(self.count_block(id))))
            .collect();
        let mut cache = match self.site_blocks.write() {
            Ok(cache) => cache,
            Err(poisoned) => poisoned.into_inner(),
        };
        cache.extend(counted);
        ids.map(|id| Arc::clone(&cache[&id])).collect()
    }

    fn count_block(&self, id: usize) -> SiteBlock {
        let start = id * STAT_BLOCK;
        let width = self.length.saturating_sub(start).min(STAT_BLOCK);
        let unknown = if self.is_dna { DNA_UNKNOWN } else { AA_UNKNOWN };
        let mut gaps = vec![0u32; width];
        let mut first = vec![0u8; width];
        let mut variable = vec![false; width];
        let mut cells = vec![0u8; width];
        for (row, &len) in self.lengths.iter().enumerate() {
            let n = len.saturating_sub(start).min(width);
            // Sites past the row's end count as gaps.
            gaps[n..].iter_mut().for_each(|g| *g += 1);
            let cells = &mut cells[..n];
            self.store.read(row, start, cells);
            for (i, &code) in cells.iter().enumerate() {
                if is_gap(code) {
                    gaps[i] += 1;
                } else if code != unknown {
                    if first[i] == 0 {
                        first[i] = code;
                    } else if first[i] != code {
                        variable[i] = true;
                    }
                }
            }
        }
        SiteBlock { gaps, variable }
    }

    fn ungapped_counts(&self, row: usize) -> Vec<usize> {
        let len = self.lengths[row];
        let mut counts = Vec::with_capacity(len / UNGAPPED_STEP + 1);
        let mut cells = vec![0u8; UNGAPPED_STEP];
        let mut total = 0;
        let mut start = 0;
        loop {
            counts.push(total);
            if start >= len {
                break;
            }
            let end = (start + UNGAPPED_STEP).min(len);
            let cells = &mut cells[..end - start];
            self.store.read(row, start, cells);
            total += cells.iter().filter(|&&c| !is_gap(c)).count();
            start = end;
        }
        counts
    }
}

/// Whether `code` is a gap or missing data.
pub fn is_gap(code: u8) -> bool {
    code == GAP || code == MISSING
}

/// Per-site counts of a block of sites, behind the chunk statistics.
struct SiteBlock {
    /// Rows with a gap, missing data, or no residue at each site.
    gaps: Vec<u32>,
    /// Whether a site holds more than one residue.
    variable: Vec<bool>,
}

/// Where residues are read from.
enum Store {
    /// A memory-mapped FASTA file and where each sequence sits in it.
    Mapped(Mmap, Vec<FastaRecord>),
    /// Upper-case residues parsed by SEGUL, one row per taxon.
    Memory(Vec<Vec<u8>>),
}

impl Store {
    fn len(&self, row: usize) -> usize {
        match self {
            Store::Mapped(_, records) => records[row].len,
            Store::Memory(rows) => rows[row].len(),
        }
    }

    fn residue(&self, row: usize, site: usize) -> u8 {
        match self {
            Store::Mapped(map, records) => map[records[row].offset_of(site)].to_ascii_uppercase(),
            Store::Memory(rows) => rows[row][site],
        }
    }

    /// Copies upper-case residues of `row` from `site` into `out`, stopping
    /// at the row's end.
    fn read(&self, row: usize, site: usize, out: &mut [u8]) {
        match self {
            Store::Mapped(map, records) => records[row].read(map, site, out),
            Store::Memory(rows) => {
                let seq = &rows[row];
                if site < seq.len() {
                    let n = out.len().min(seq.len() - site);
                    out[..n].copy_from_slice(&seq[site..site + n]);
                }
            }
        }
    }
}

/// A sequence in a FASTA file, laid out like a `samtools faidx` record.
#[derive(Debug, PartialEq)]
struct FastaRecord {
    /// Byte offset of the first residue.
    offset: usize,
    /// Residues in the sequence.
    len: usize,
    /// Residues on every full line.
    line_bases: usize,
    /// Bytes on every full line, including the line ending.
    line_bytes: usize,
}

impl FastaRecord {
    fn offset_of(&self, site: usize) -> usize {
        self.offset + site / self.line_bases * self.line_bytes + site % self.line_bases
    }

    fn read(&self, map: &[u8], site: usize, out: &mut [u8]) {
        let end = (site + out.len()).min(self.len);
        let mut site = site;
        let mut written = 0;
        while site < end {
            let in_line = site % self.line_bases;
            let take = (self.line_bases - in_line).min(end - site);
            let from = self.offset_of(site);
            let cells = &mut out[written..written + take];
            cells.copy_from_slice(&map[from..from + take]);
            cells.make_ascii_uppercase();
            site += take;
            written += take;
        }
    }
}

/// Byte offsets of every sequence in a FASTA file.
struct FastaIndex {
    map: Mmap,
    names: Vec<String>,
    records: Vec<FastaRecord>,
}

impl FastaIndex {
    /// Indexes a FASTA file, or none when its lines are too irregular to
    /// index, such as wrapped lines of different widths.
    fn build(path: &Path) -> anyhow::Result<Option<FastaIndex>> {
        let file = File::open(path)
            .map_err(|e| anyhow::anyhow!("Failed to open {}: {e}", path.display()))?;
        if file.metadata()?.len() == 0 {
            anyhow::bail!("Failed to read {}: the file is empty", path.display());
        }
        // SAFETY: the map is only read. Changing the file while it is open
        // shows changed residues; viewers reopen files to refresh them.
        let map = unsafe { Mmap::map(&file)? };
        Ok(Self::scan(&map).map(|(names, records)| FastaIndex {
            map,
            names,
            records,
        }))
    }

    fn scan(bytes: &[u8]) -> Option<(Vec<String>, Vec<FastaRecord>)> {
        let mut names = Vec::new();
        let mut records = Vec::new();
        let mut lines = Lines { bytes, pos: 0 };
        let mut next = lines.next();
        // Blank lines before the first header.
        while next
            .as_ref()
            .is_some_and(|l| l.content.trim_ascii().is_empty())
        {
            next = lines.next();
        }
        while let Some(header) = next {
            let id = header.content.strip_prefix(b">")?;
            names.push(String::from_utf8_lossy(id.trim_ascii()).into_owned());
            let (record, after) = Self::record(&mut lines)?;
            records.push(record);
            next = after;
        }
        if records.is_empty() {
            return None;
        }
        Some((names, records))
    }

    // Reads the sequence lines after a header, returning the record and the
    // next header.
    fn record<'a>(lines: &mut Lines<'a>) -> Option<(FastaRecord, Option<Line<'a>>)> {
        let mut record = FastaRecord {
            offset: lines.pos,
            len: 0,
            line_bases: 0,
            line_bytes: 0,
        };
        // A line shorter than the first ends the sequence.
        let mut is_closed = false;
        for line in lines.by_ref() {
            let bases = line.content.len();
            if line.content.first() == Some(&b'>') {
                if record.line_bases == 0 {
                    record.line_bases = 1;
                }
                return Some((record, Some(line)));
            }
            if bases == 0 {
                is_closed = true;
                continue;
            }
            let is_padded = line.content[0].is_ascii_whitespace()
                || line.content[bases - 1].is_ascii_whitespace();
            if is_closed || is_padded {
                return None;
            }
            if record.line_bases == 0 {
                record.offset = line.start;
                record.line_bases = bases;
                record.line_bytes = line.bytes;
            } else {
                // The last line of a file may have no line ending.
                let ending = line.bytes - bases;
                let is_other_ending =
                    ending != 0 && ending != record.line_bytes - record.line_bases;
                if bases > record.line_bases || is_other_ending {
                    return None;
                }
            }
            if bases < record.line_bases {
                is_closed = true;
            }
            record.len += bases;
        }
        if record.line_bases == 0 {
            record.line_bases = 1;
        }
        Some((record, None))
    }
}

/// A line of a file, without its line ending.
struct Line<'a> {
    start: usize,
    content: &'a [u8],
    /// Bytes including the line ending.
    bytes: usize,
}

struct Lines<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> Iterator for Lines<'a> {
    type Item = Line<'a>;

    fn next(&mut self) -> Option<Line<'a>> {
        if self.pos >= self.bytes.len() {
            return None;
        }
        let start = self.pos;
        let rest = &self.bytes[start..];
        let (mut content, bytes) = match memchr(b'\n', rest) {
            Some(i) => (&rest[..i], i + 1),
            None => (rest, rest.len()),
        };
        if let Some(stripped) = content.strip_suffix(b"\r") {
            content = stripped;
        }
        self.pos += bytes;
        Some(Line {
            start,
            content,
            bytes,
        })
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use std::fs;

    use tempfile::TempDir;

    fn rows(view: &AlignmentView) -> Vec<Vec<u8>> {
        let width = view.length;
        let tile = view.tile(0, view.rows(), 0, width);
        tile.chunks(width.max(1))
            .zip(&view.lengths)
            .map(|(row, &len)| row[..len].to_vec())
            .collect()
    }

    // The view comes first so it drops before the directory; Windows cannot
    // remove a mapped file.
    fn open(name: &str, content: &str) -> (AlignmentView, TempDir) {
        let dir = TempDir::with_prefix("view").expect("Failed creating a temp dir");
        let path = dir.path().join(name);
        fs::write(&path, content).expect("Failed writing a temp file");
        let view = AlignmentView::open(&path, &InputFmt::Auto, &DataType::Ignore)
            .expect("Failed opening a temp file");
        (view, dir)
    }

    #[test]
    fn test_open_fasta() {
        let path = Path::new("tests/files/simple.fas");
        let view = AlignmentView::open(path, &InputFmt::Auto, &DataType::Ignore)
            .expect("Failed loading fixture");
        let summary = view.summary();
        assert!(matches!(view.store, Store::Mapped(..)));
        assert_eq!(summary.names, vec!["ABCD", "CDEF"]);
        assert_eq!(rows(&view)[0], b"AGTATG".to_vec());
        assert_eq!(summary.length, 6);
        assert!(summary.is_aligned);
        assert!(summary.is_dna);
    }

    #[test]
    fn test_open_nexus_fixture() {
        let path = Path::new("tests/files/alignments/gene_1.nex");
        let view = AlignmentView::open(path, &InputFmt::Nexus, &DataType::Dna)
            .expect("Failed loading fixture");
        assert!(matches!(view.store, Store::Memory(..)));
        assert!(view.rows() > 0);
        assert!(view.is_aligned());
    }

    #[test]
    fn test_index_matches_segul() {
        let cases = [
            ("single.fas", ">a\nacgt-a\n>b desc\nACGTTA\n"),
            ("wrapped.fas", ">a\nACG\nTAC\nG\n>b\nACG\nTAC\nGT\n\n"),
            ("crlf.fas", ">a\r\nACG\r\nTA\r\n>b\r\nAC-\r\nTA\r\n"),
            ("ragged.fas", ">a\nACGTACGT\n>b\nAC\n>c\n\n"),
            ("no_newline.fas", ">a\nACG\nTA\n>b\nACG\nTA"),
        ];
        for (name, content) in cases {
            let (mapped, dir) = open(name, content);
            assert!(
                matches!(mapped.store, Store::Mapped(..)),
                "{name} was not indexed"
            );
            let (names, store) = AlignmentView::parse(&dir.path().join(name), &InputFmt::Fasta)
                .expect("Failed parsing");
            let parsed = AlignmentView::new(names, store);
            assert_eq!(mapped.names, parsed.names, "{name}");
            assert_eq!(rows(&mapped), rows(&parsed), "{name}");
        }
    }

    #[test]
    fn test_irregular_fasta_falls_back() {
        let (view, _dir) = open("irregular.fas", ">a\nAC\nACGT\n>b\nACGTAC\n");
        assert!(matches!(view.store, Store::Memory(..)));
        assert_eq!(rows(&view)[0], b"ACACGT".to_vec());
    }

    #[test]
    fn test_upper_case_extension() {
        let (view, _dir) = open("upper.FAS", ">a\nACGT\n");
        assert_eq!(view.length(), 4);
    }

    #[test]
    fn test_open_phylip_and_nexus() {
        let (phylip, _dir) = open("seq.phy", "2 4\na ACGT\nb AC-T\n");
        assert_eq!(rows(&phylip), vec![b"ACGT".to_vec(), b"AC-T".to_vec()]);
        let (nexus, _dir) = open(
            "seq.nex",
            "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=4;\n\
             format datatype=dna missing=? gap=- interleave;\nmatrix\n\
             a AC\nb AC\n\na GT\nb -T\n;\nend;\n",
        );
        assert_eq!(rows(&nexus), vec![b"ACGT".to_vec(), b"AC-T".to_vec()]);
    }

    #[test]
    fn test_tile_pads_past_row_end() {
        let (view, _dir) = open("pad.fas", ">a\nACG\nTA\n>b\nAC\n");
        assert_eq!(view.tile(0, 2, 1, 4), b"CGTAC\0\0\0".to_vec());
        assert_eq!(view.tile(1, 5, 0, 2), b"AC".to_vec());
        assert!(view.tile(5, 2, 0, 2).is_empty());
    }

    #[test]
    fn test_residue() {
        let (view, _dir) = open("residue.fas", ">a\nacg\n>b\nAC\n");
        assert_eq!(view.residue(0, 1), Some(b'C'));
        assert_eq!(view.residue(1, 2), None);
        assert_eq!(view.residue(2, 0), None);
    }

    #[test]
    fn test_sample() {
        let (view, _dir) = open("sample.fas", ">a\nACGTACGT\n>b\nTT\n");
        assert_eq!(view.sample(&[1, 0], 4), b"T\0\0\0AGAG".to_vec());
    }

    #[test]
    fn test_chunk_stats() {
        let (view, _dir) = open("stats.fas", ">a\nACGT--\n>b\nATGTN-\n>c\nACGA\n");
        let stats = view.chunk_stats(2, 0, usize::MAX);
        assert_eq!(stats.len(), 3);
        assert_eq!((stats[0].start, stats[0].end), (0, 2));
        assert_eq!(stats[0].variable_fraction, 0.5);
        assert_eq!(stats[0].gap_fraction, 0.0);
        assert_eq!(stats[1].variable_fraction, 0.5);
        assert_eq!(stats[2].gap_fraction, 5.0 / 6.0);
        assert_eq!(stats[2].variable_fraction, 0.0);
    }

    #[test]
    fn test_chunk_stats_ignore_unknowns() {
        let (view, _dir) = open("unknown.fas", ">a\nACGT\n>b\nACGA\n>c\nA-N-\n");
        let stats = view.chunk_stats(2, 0, usize::MAX);
        assert_eq!(stats.len(), 2);
        // Sites 1-2: one gap in six cells, no variable site.
        assert!((stats[0].gap_fraction - 1.0 / 6.0).abs() < 1e-9);
        assert_eq!(stats[0].variable_fraction, 0.0);
        // Sites 3-4: N is ignored, so only site 4 (T, A) varies.
        assert!((stats[1].gap_fraction - 1.0 / 6.0).abs() < 1e-9);
        assert_eq!(stats[1].variable_fraction, 0.5);
    }

    #[test]
    fn test_chunk_stats_count_short_rows_as_gaps() {
        let (view, _dir) = open("short.fas", ">a\nACGT\n>b\nAC\n");
        let stats = view.chunk_stats(4, 0, usize::MAX);
        assert_eq!(stats.len(), 1);
        assert_eq!(stats[0].gap_fraction, 0.25);
        // Changing the chunk size reuses the per-site counts.
        assert_eq!(view.chunk_stats(1, 0, usize::MAX).len(), 4);
    }

    #[test]
    fn test_chunk_stats_window_across_blocks() {
        // Row b is variable at every 3rd site and stops before the end.
        let a = "A".repeat(3 * STAT_BLOCK);
        let b: String = (0..3 * STAT_BLOCK - 7)
            .map(|i| if i % 3 == 0 { 'C' } else { 'A' })
            .collect();
        let (view, _dir) = open("blocks.fas", &format!(">a\n{a}\n>b\n{b}\n"));
        let size = 1000;
        let all = view.chunk_stats(size, 0, usize::MAX);
        assert_eq!(all.len(), (3 * STAT_BLOCK).div_ceil(size));
        // A window of chunks matches the same chunks of the whole.
        let window = view.chunk_stats(size, 3, 6);
        assert_eq!(window, all[3..9].to_vec());
        // Chunk 4 (sites 4000-4999) crosses the first block boundary.
        assert_eq!(all[4].variable_fraction, 0.333);
        let last = all.last().unwrap();
        assert_eq!(last.end, 3 * STAT_BLOCK);
        let sites = (last.end - last.start) as f64;
        assert_eq!(last.gap_fraction, 7.0 / (2.0 * sites));
        assert!(view.chunk_stats(size, 999, 5).is_empty());
    }

    #[test]
    fn test_ungapped_position() {
        let (view, _dir) = open("ungapped.fas", ">a\nA--C?G\n");
        let positions: Vec<Option<usize>> = (0..7).map(|s| view.ungapped_position(0, s)).collect();
        assert_eq!(
            positions,
            vec![Some(1), None, None, Some(2), None, Some(3), None]
        );
    }

    #[test]
    fn test_ungapped_position_across_steps() {
        let seq = "A-".repeat(UNGAPPED_STEP + 10);
        let (view, _dir) = open("long.fas", &format!(">a\n{seq}\n"));
        let site = 2 * (UNGAPPED_STEP + 5);
        assert_eq!(view.ungapped_position(0, site), Some(UNGAPPED_STEP + 6));
    }

    #[test]
    fn test_fasta() {
        let (view, _dir) = open("copy.fas", ">a\nac-t\n>b\nGGCC\n>c\nTTAA\n");
        assert_eq!(view.fasta(&[2, 0]), ">a\nAC-T\n>c\nTTAA\n");
    }

    #[test]
    fn test_amino_acid_detection() {
        let path = Path::new("tests/files/simple_aa.fas");
        let view = AlignmentView::open(path, &InputFmt::Fasta, &DataType::Ignore)
            .expect("Failed loading fixture");
        assert!(!view.is_dna());
        assert!(view.is_aligned());
    }

    #[test]
    fn test_datatype_override() {
        let path = Path::new("tests/files/simple_aa.fas");
        let view = AlignmentView::open(path, &InputFmt::Fasta, &DataType::Dna)
            .expect("Failed loading fixture");
        assert!(view.is_dna());
    }

    #[test]
    fn test_unaligned_length() {
        let (view, _dir) = open("unaligned.fas", ">a\nACGT\n>b\nACGTAC\n");
        assert_eq!(view.length(), 6);
        assert!(!view.is_aligned());
    }

    #[test]
    fn test_unknown_extension() {
        let open = AlignmentView::open(Path::new("seq.txt"), &InputFmt::Auto, &DataType::Ignore);
        assert!(open.is_err());
    }

    #[test]
    fn test_missing_file() {
        let path = Path::new("tests/files/missing.nex");
        assert!(AlignmentView::open(path, &InputFmt::Auto, &DataType::Ignore).is_err());
        let path = Path::new("tests/files/missing.fas");
        assert!(AlignmentView::open(path, &InputFmt::Auto, &DataType::Ignore).is_err());
    }
}
