//! Data types for supported formats
use std::ffi::OsStr;
use std::path::Path;

use ahash::AHashMap as HashMap;
use indexmap::IndexMap;
use serde::{Deserialize, Serialize};

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum GenomicFmt {
    /// Fastq format with auto-detection
    /// to inter gzipped and non-gzipped files
    FastqAuto,
    /// Fastq format
    /// with gzipped files
    FastqGzip,
    /// Fastq format
    /// with non-gzipped files
    Fastq,
    /// Contig format with auto-detection
    /// to inter gzipped and non-gzipped files
    /// The contig file should be in FASTA format
    ContigAuto,
    /// Contig format
    /// with gzipped files
    ContigGzip,
    /// Contig format
    /// with non-gzipped files
    Contig,
    /// Multi-Alignment Format (MAF) format
    Maf,
}

impl std::fmt::Display for GenomicFmt {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::FastqAuto => write!(f, "fastq-auto"),
            Self::FastqGzip => write!(f, "fastq-gzip"),
            Self::Fastq => write!(f, "fastq"),
            Self::ContigAuto => write!(f, "contig-auto"),
            Self::ContigGzip => write!(f, "contig-gzip"),
            Self::Contig => write!(f, "contig"),
            Self::Maf => write!(f, "maf"),
        }
    }
}

impl std::str::FromStr for GenomicFmt {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "fastq-auto" => Ok(Self::FastqAuto),
            "fastq-gzip" => Ok(Self::FastqGzip),
            "fastq" => Ok(Self::Fastq),
            "contig-auto" => Ok(Self::ContigAuto),
            "contig-gzip" => Ok(Self::ContigGzip),
            "contig" => Ok(Self::Contig),
            "maf" => Ok(Self::Maf),
            _ => Err(format!("{} is not a valid format", s)),
        }
    }
}

/// Data types for high-throughput sequencing reads
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum SeqReadFmt {
    /// Infer format from file extension
    Auto,
    /// Fastq format
    Fastq,
    /// Gzip compressed fastq format
    Gzip,
}

impl std::fmt::Display for SeqReadFmt {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Auto => write!(f, "auto"),
            Self::Fastq => write!(f, "fastq"),
            Self::Gzip => write!(f, "gzip"),
        }
    }
}

impl std::str::FromStr for SeqReadFmt {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "auto" => Ok(Self::Auto),
            "fastq" => Ok(Self::Fastq),
            "gzip" => Ok(Self::Gzip),
            _ => Err(format!("{} is not a valid format", s)),
        }
    }
}

/// Quality score encoding of FASTQ reads.
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum PhredOffset {
    /// Detect the offset from the quality scores of the first reads
    Auto,
    /// Sanger and Illumina 1.8+ (ASCII 33)
    Phred33,
    /// Illumina 1.3 to 1.7 (ASCII 64)
    Phred64,
}

impl std::fmt::Display for PhredOffset {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Auto => write!(f, "auto"),
            Self::Phred33 => write!(f, "33"),
            Self::Phred64 => write!(f, "64"),
        }
    }
}

impl std::str::FromStr for PhredOffset {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "auto" => Ok(Self::Auto),
            "33" => Ok(Self::Phred33),
            "64" => Ok(Self::Phred64),
            _ => Err(format!("{} is not a valid phred offset", s)),
        }
    }
}

pub fn infer_raw_input_auto(input: &Path) -> SeqReadFmt {
    let ext: &str = input
        .extension()
        .and_then(OsStr::to_str)
        .expect("Failed parsing extension");
    match ext {
        "fq" | "fastq" => SeqReadFmt::Fastq,
        "gz" | "gzip" => SeqReadFmt::Gzip,
        _ => panic!(
            "The program cannot recognize the file extension. \
        Try to specify the input format."
        ),
    }
}

/// Data type for contig sequences
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum ContigFmt {
    /// Infer format from file extension
    Auto,
    /// Fasta format
    Fasta,
    /// Gzip compressed fasta format
    Gzip,
}

impl std::fmt::Display for ContigFmt {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Auto => write!(f, "auto"),
            Self::Fasta => write!(f, "fasta"),
            Self::Gzip => write!(f, "gzip"),
        }
    }
}

pub fn infer_contig_fmt_auto(input: &Path) -> ContigFmt {
    let ext: &str = input
        .extension()
        .and_then(OsStr::to_str)
        .expect("Failed parsing extension");
    match ext {
        "fa" | "fasta" | "fna" | "fsa" | "fas" => ContigFmt::Fasta,
        "gz" | "gzip" => ContigFmt::Gzip,
        _ => panic!(
            "The program cannot recognize the file extension. \
        Try to specify the input format."
        ),
    }
}

impl std::str::FromStr for ContigFmt {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "auto" => Ok(Self::Auto),
            "fasta" => Ok(Self::Fasta),
            "gzip" => Ok(Self::Gzip),
            _ => Err(format!("{} is not a valid format", s)),
        }
    }
}

/// Data types for input sequences
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum InputFmt {
    /// Infer format from file extension
    Auto,
    /// Fasta format
    Fasta,
    /// Nexus format
    Nexus,
    /// Phylip format
    Phylip,
}

impl std::fmt::Display for InputFmt {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Auto => write!(f, "auto"),
            Self::Fasta => write!(f, "fasta"),
            Self::Nexus => write!(f, "nexus"),
            Self::Phylip => write!(f, "phylip"),
        }
    }
}

impl std::str::FromStr for InputFmt {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "auto" => Ok(Self::Auto),
            "fasta" => Ok(Self::Fasta),
            "nexus" => Ok(Self::Nexus),
            "phylip" => Ok(Self::Phylip),
            _ => Err(format!("{} is not a valid format", s)),
        }
    }
}

/// Infer input format automatically based on the file extension.
/// Return the input format.
/// # Example
/// ```
/// use std::path::Path;
/// use segul_core::helper::types::InputFmt;
/// use segul_core::helper::types::infer_input_auto;
///
/// let file = Path::new("tests/files/simple.fas");
/// let input_fmt = infer_input_auto(&file);
/// assert_eq!(input_fmt, InputFmt::Fasta);
/// ```
pub fn infer_input_auto(input: &Path) -> InputFmt {
    let ext: &str = input
        .extension()
        .and_then(OsStr::to_str)
        .expect("Failed parsing extension");
    match ext {
        "fa" | "fasta" | "fna" | "fsa" | "fas" => InputFmt::Fasta,
        "nex" | "nxs" | "nexus" => InputFmt::Nexus,
        "phy" | "phylip" | "ph" => InputFmt::Phylip,
        _ => panic!(
            "The program cannot recognize the file extension. \
        Try to specify the input format."
        ),
    }
}

/// Data types for output sequences
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum OutputFmt {
    /// Fastq format
    Fasta,
    /// Nexus format
    Nexus,
    /// Phylip format
    Phylip,
    /// Interleaved Fasta format
    FastaInt,
    /// Interleaved Nexus format
    NexusInt,
    /// Interleaved Phylip format
    PhylipInt,
}

impl std::fmt::Display for OutputFmt {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Fasta => write!(f, "fasta"),
            Self::Nexus => write!(f, "nexus"),
            Self::Phylip => write!(f, "phylip"),
            Self::FastaInt => write!(f, "Interleaved fasta"),
            Self::NexusInt => write!(f, "Interleaved nexus"),
            Self::PhylipInt => write!(f, "Interleaved phylip"),
        }
    }
}

impl std::str::FromStr for OutputFmt {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "fasta" => Ok(Self::Fasta),
            "nexus" => Ok(Self::Nexus),
            "phylip" => Ok(Self::Phylip),
            "fasta-int" => Ok(Self::FastaInt),
            "nexus-int" => Ok(Self::NexusInt),
            "phylip-int" => Ok(Self::PhylipInt),
            _ => Err(format!("{} is not a valid format", s)),
        }
    }
}

/// Data types for alignment partitions
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum PartitionFmt {
    /// Nexus charset format embedded in a nexus alignment file
    Charset,
    /// Nexus charset format embedded in nexus file for codon model partitions
    CharsetCodon,
    /// Nexus format separate from nexus alignment files
    Nexus,
    /// Nexus format separate from nexus alignment files for codon model partitions
    NexusCodon,
    /// RAxML format partition file
    Raxml,
    /// RAxML format partition file for codon model partitions
    RaxmlCodon,
}

impl std::fmt::Display for PartitionFmt {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Charset => write!(f, "charset"),
            Self::CharsetCodon => write!(f, "charset-codon"),
            Self::Nexus => write!(f, "nexus"),
            Self::NexusCodon => write!(f, "nexus-codon"),
            Self::Raxml => write!(f, "raxml"),
            Self::RaxmlCodon => write!(f, "raxml-codon"),
        }
    }
}

impl std::str::FromStr for PartitionFmt {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "charset" => Ok(Self::Charset),
            "charset-codon" => Ok(Self::CharsetCodon),
            "nexus" => Ok(Self::Nexus),
            "nexus-codon" => Ok(Self::NexusCodon),
            "raxml" => Ok(Self::Raxml),
            "raxml-codon" => Ok(Self::RaxmlCodon),
            _ => Err(format!("{} is not a valid format", s)),
        }
    }
}

#[derive(Debug, Clone, Copy, Default, Eq, PartialEq)]
pub enum DnaStrand {
    Forward,
    Reverse,
    #[default]
    Missing,
}

impl std::fmt::Display for DnaStrand {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Forward => write!(f, "+"),
            Self::Reverse => write!(f, "-"),
            Self::Missing => write!(f, "."),
        }
    }
}

impl std::str::FromStr for DnaStrand {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Self::Forward),
            "-" => Ok(Self::Reverse),
            "." => Ok(Self::Missing),
            _ => Err(format!("{} is not a valid DNA strand", s)),
        }
    }
}

impl DnaStrand {
    pub fn from_char(c: char) -> Self {
        match c {
            '+' => DnaStrand::Forward,
            '-' => DnaStrand::Reverse,
            '.' => DnaStrand::Missing,
            _ => panic!("Invalid DNA strand"),
        }
    }

    pub fn to_char(&self) -> char {
        match self {
            DnaStrand::Forward => '+',
            DnaStrand::Reverse => '-',
            DnaStrand::Missing => '.',
        }
    }
}

impl Serialize for DnaStrand {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        match *self {
            DnaStrand::Forward => serializer.serialize_str("+"),
            DnaStrand::Reverse => serializer.serialize_str("-"),
            DnaStrand::Missing => serializer.serialize_str("."),
        }
    }
}

impl<'de> Deserialize<'de> for DnaStrand {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        match s.as_str() {
            "+" => Ok(DnaStrand::Forward),
            "-" => Ok(DnaStrand::Reverse),
            _ => Err(serde::de::Error::custom("Invalid DNA strand")),
        }
    }
}

/// Data types for sequence data
#[derive(PartialEq)]
pub enum DataType {
    /// DNA sequences
    Dna,
    /// Amino acid sequences
    Aa,
    /// Ignore type when parsing sequences
    Ignore,
}

impl std::fmt::Display for DataType {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Dna => write!(f, "DNA"),
            Self::Aa => write!(f, "Amino acid"),
            Self::Ignore => write!(f, "Ignore Data Type"),
        }
    }
}

impl std::str::FromStr for DataType {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "dna" => Ok(Self::Dna),
            "aa" => Ok(Self::Aa),
            "ignore" => Ok(Self::Ignore),
            _ => Err(format!("{} is not a valid data type", s)),
        }
    }
}

/// Data types for sequence data
pub type SeqMatrix = IndexMap<String, String>;

/// A 1-based inclusive site range with a positive stride (`step >= 1`).
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct SiteRange {
    /// 1-based inclusive start position.
    pub start: usize,
    /// 1-based inclusive end position.
    pub end: usize,
    /// Stride between sites (1 for contiguous ranges, 3 for codon positions).
    pub step: usize,
}

impl SiteRange {
    /// Create a contiguous 1-based inclusive range (`step = 1`).
    pub fn new(start: usize, end: usize) -> Self {
        Self {
            start,
            end,
            step: 1,
        }
    }

    /// Create a strided 1-based inclusive range.
    pub fn with_step(start: usize, end: usize, step: usize) -> Self {
        Self { start, end, step }
    }

    /// Create a single-site range (`start = end = pos`, `step = 1`).
    pub fn single(pos: usize) -> Self {
        Self {
            start: pos,
            end: pos,
            step: 1,
        }
    }

    /// Return true when `start >= 1`, `start <= end`, and `step >= 1`.
    pub fn is_valid(&self) -> bool {
        self.start >= 1 && self.start <= self.end && self.step >= 1
    }

    /// Number of sites selected by this range.
    pub fn site_count(&self) -> usize {
        if !self.is_valid() {
            0
        } else {
            (self.end - self.start) / self.step + 1
        }
    }

    /// Iterate over all 1-based site positions in this range.
    pub fn iter_sites(&self) -> std::iter::StepBy<std::ops::RangeInclusive<usize>> {
        (self.start..=self.end).step_by(self.step.max(1))
    }
}

impl std::fmt::Display for SiteRange {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.step > 1 {
            write!(f, "{}-{}\\{}", self.start, self.end, self.step)
        } else if self.start == self.end {
            write!(f, "{}", self.start)
        } else {
            write!(f, "{}-{}", self.start, self.end)
        }
    }
}

/// Data types for partition data
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Partition {
    /// Gene, locus, or subset name
    pub gene: String,
    /// Minimum 1-based start position across all ranges
    pub start: usize,
    /// Maximum 1-based end position across all ranges
    pub end: usize,
    /// One or more 1-based inclusive site ranges belonging to this partition
    pub ranges: Vec<SiteRange>,
    /// Optional substitution model label (e.g., `GTR+I+G`)
    pub model: Option<String>,
}

impl Default for Partition {
    fn default() -> Self {
        Self::new()
    }
}

impl Partition {
    pub fn new() -> Self {
        Self {
            gene: String::new(),
            start: 0,
            end: 0,
            ranges: Vec::new(),
            model: None,
        }
    }

    /// Create a partition covering a single contiguous range.
    pub fn from_range(gene: impl Into<String>, start: usize, end: usize) -> Self {
        Self {
            gene: gene.into(),
            start,
            end,
            ranges: vec![SiteRange::new(start, end)],
            model: None,
        }
    }

    /// Create a partition from one or more site ranges and an optional model.
    pub fn from_ranges(
        gene: impl Into<String>,
        ranges: Vec<SiteRange>,
        model: Option<String>,
    ) -> Self {
        let start = ranges.iter().map(|r| r.start).min().unwrap_or(0);
        let end = ranges.iter().map(|r| r.end).max().unwrap_or(0);
        Self {
            gene: gene.into(),
            start,
            end,
            ranges,
            model,
        }
    }

    /// Return `ranges` when present, or fall back to `[start..=end]` for legacy structs.
    pub fn effective_ranges(&self) -> std::borrow::Cow<'_, [SiteRange]> {
        if !self.ranges.is_empty() {
            std::borrow::Cow::Borrowed(&self.ranges)
        } else if self.start > 0 && self.end >= self.start {
            std::borrow::Cow::Owned(vec![SiteRange::new(self.start, self.end)])
        } else {
            std::borrow::Cow::Borrowed(&[])
        }
    }

    /// Return true if any range in this partition uses a stride greater than 1.
    pub fn has_stride(&self) -> bool {
        self.effective_ranges().iter().any(|r| r.step > 1)
    }

    /// Return the maximum 1-based site position referenced by this partition.
    pub fn max_site(&self) -> usize {
        self.effective_ranges()
            .iter()
            .map(|r| r.end)
            .max()
            .unwrap_or(self.end)
    }

    /// Expand all ranges into sorted, deduplicated 1-based site positions.
    pub fn sorted_site_positions(&self) -> Vec<usize> {
        let ranges = self.effective_ranges();
        let total: usize = ranges.iter().map(SiteRange::site_count).sum();
        let mut sites = Vec::with_capacity(total);
        for range in ranges.iter() {
            sites.extend(range.iter_sites());
        }
        sites.sort_unstable();
        sites.dedup();
        sites
    }

    /// Group sorted unique sites into maximal 0-based half-open `[start, end)` intervals.
    pub fn contiguous_site_runs(&self) -> Vec<(usize, usize)> {
        let sites = self.sorted_site_positions();
        let mut runs = Vec::new();
        let mut iter = sites.into_iter();
        let Some(first) = iter.next() else {
            return runs;
        };
        let mut run_start = first - 1;
        let mut run_end = first;
        for site in iter {
            let zero_idx = site - 1;
            if zero_idx == run_end {
                run_end += 1;
            } else {
                runs.push((run_start, run_end));
                run_start = zero_idx;
                run_end = site;
            }
        }
        runs.push((run_start, run_end));
        runs
    }
}

/// Data types for header data
/// ntax: Number of taxa
/// nchar: Number of characters
/// datatype: Data type
/// missing: Missing data character
/// gap: Gap character
/// aligned: Aligned or unaligned
#[derive(Clone, Debug)]
pub struct Header {
    /// Number of taxa
    pub ntax: usize,
    /// Number of characters
    pub nchar: usize,
    /// Data type
    pub datatype: String,
    /// Missing data character: ?
    pub missing: char,
    /// Gap character: -
    pub gap: char,
    /// Aligned or unaligned
    pub aligned: bool,
}

impl Default for Header {
    fn default() -> Self {
        Self::new()
    }
}

impl Header {
    pub fn new() -> Self {
        Self {
            ntax: 0,
            nchar: 0,
            datatype: String::from("dna"),
            missing: '?',
            gap: '-',
            aligned: false,
        }
    }

    pub fn update(&mut self, matrix: &SeqMatrix) {
        self.ntax = matrix.len();
        self.nchar = matrix.values().next().unwrap().len();
    }

    pub fn from_seq_matrix(&mut self, matrix: &SeqMatrix, aligned: bool) {
        self.ntax = matrix.len();
        self.nchar = matrix.values().next().unwrap().len();
        self.aligned = aligned;
    }
}

/// Data types for taxon data
pub struct TaxonRecords {
    /// Character counts
    pub char_counts: HashMap<char, usize>,
    /// Number of loci
    pub locus_counts: usize,
    /// Total number of characters
    pub total_chars: usize,
    /// Number of G/C characters
    pub gc_count: usize,
    /// Number of A/T characters
    pub at_count: usize,
    /// Number of nucleotides
    pub nucleotides: usize,
    /// Number of missing data
    pub missing_data: usize,
}

impl Default for TaxonRecords {
    fn default() -> Self {
        Self::new()
    }
}

impl TaxonRecords {
    pub fn new() -> Self {
        Self {
            char_counts: HashMap::new(),
            locus_counts: 0,
            total_chars: 0,
            gc_count: 0,
            at_count: 0,
            nucleotides: 0,
            missing_data: 0,
        }
    }
}

#[derive(Debug, Clone, Copy, Default, Eq, PartialEq)]
pub enum SummaryMode {
    /// Only write read counts.
    Minimal,
    /// Write file-level and per-read summary statistics to a single CSV.
    Basic,
    /// Write all summary statistics, per-position profiles, distributions,
    /// duplication, adapters, overrepresented sequences, and JSON report.
    #[default]
    Complete,
}

impl std::fmt::Display for SummaryMode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Minimal => write!(f, "minimal"),
            Self::Basic => write!(f, "basic"),
            Self::Complete => write!(f, "complete"),
        }
    }
}

impl std::str::FromStr for SummaryMode {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "minimal" => Ok(Self::Minimal),
            "basic" | "default" => Ok(Self::Basic),
            "complete" => Ok(Self::Complete),
            _ => Err(format!("{} is not a valid summary mode", s)),
        }
    }
}

pub enum SummaryOutput {
    /// Print summary to stdout
    Stdout,
    /// Print summary comma-delimited file
    Csv,
    /// Print summary to tab-delimited file
    Tsv,
    /// Print to json file
    Json,
}

/// Data types for genetic codes
/// Based on NCBI genetic code table
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum GeneticCodes {
    /// Ncbi Table 1
    StandardCode,
    /// Ncbi Table 2
    VertMtDna,
    /// Ncbi Table 3
    YeastMtDna,
    /// Ncbi Table 4
    MoldProtCoelMtDna,
    /// Ncbi Table 5
    InvertMtDna,
    /// Ncbi Table 6
    CilDasHexNu,
    /// Ncbi Table 9
    EchiFlatwormMtDna,
    /// Ncbi Table 10
    EuplotidNu,
    /// Ncbi Table 11
    BacArchPlantPlast,
    /// Ncbi Table 12
    AltYeastNu,
    /// Ncbi Table 13
    AsciMtDna,
    /// Ncbi Table 14
    AltFlatwormMtDna,
    /// Ncbi Table 16
    ChlorMtDna,
    /// Ncbi Table 21
    TrematodeMtDna,
    /// Ncbi Table 22
    ScenedesmusMtDna,
    /// Ncbi Table 23
    ThrausMtDna,
    /// Ncbi Table 24
    RhabdopMtDna,
    /// Ncbi Table 25
    CaDivSR1GraciBac,
    /// Ncbi Table 26
    PachyNu,
    /// Ncbi Table 29
    MesodiniumNu,
    /// Ncbi Table 30
    PeritrichNu,
    /// Ncbi Table 33
    CephalodiscidaeMtDna,
}

impl std::fmt::Display for GeneticCodes {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::StandardCode => write!(f, "Standard Code"),
            Self::VertMtDna => write!(f, "Vertebrate Mitochondrial DNA"),
            Self::YeastMtDna => write!(f, "Yeast Mitochondrial DNA"),
            Self::MoldProtCoelMtDna => write!(
                f,
                "Mold, Protozoan, and Coelenterate Mitochondrial DNA and the Mycoplasma/Spiroplasma Code"
            ),
            Self::InvertMtDna => write!(f, "Invertebrate Mitochondrial DNA"),
            Self::CilDasHexNu => write!(f, "Ciliate, Dasycladacean and Hexamita Nuclear Code"),
            Self::EchiFlatwormMtDna => write!(f, "Echinoderm and Flatworm Mitochondrial DNA"),
            Self::EuplotidNu => write!(f, "Euplotid Nuclear Code"),
            Self::BacArchPlantPlast => write!(f, "Bacterial, Archaeal and Plant Plastid Code"),
            Self::AltYeastNu => write!(f, "Alternative Yeast Nuclear Code"),
            Self::AsciMtDna => write!(f, "Ascidian Mitochondrial DNA"),
            Self::AltFlatwormMtDna => write!(f, "Alternative Flatworm Mitochondrial DNA"),
            Self::ChlorMtDna => write!(f, "Chlorophycean Mitochondrial DNA"),
            Self::TrematodeMtDna => write!(f, "Trematode Mitochondrial DNA"),
            Self::ScenedesmusMtDna => write!(f, "Scenedesmus obliquus Mitochondrial DNA"),
            Self::ThrausMtDna => write!(f, "Thraustochytrium Mitochondrial DNA"),
            Self::RhabdopMtDna => write!(f, "Rhabdopleuridae Mitochondrial DNA"),
            Self::CaDivSR1GraciBac => write!(f, "Candidate Division SR1 and Gracilibacteria"),
            Self::PachyNu => write!(f, "Pachysolen tannophilus Nuclear Code"),
            Self::MesodiniumNu => write!(f, "Mesodinium Nuclear Code"),
            Self::PeritrichNu => write!(f, "Peritrich Nuclear Code"),
            Self::CephalodiscidaeMtDna => write!(f, "Cephalodiscidae Mitochondrial DNA"),
        }
    }
}

impl std::str::FromStr for GeneticCodes {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "1" => Ok(Self::StandardCode),
            "2" => Ok(Self::VertMtDna),
            "3" => Ok(Self::YeastMtDna),
            "4" => Ok(Self::MoldProtCoelMtDna),
            "5" => Ok(Self::InvertMtDna),
            "6" => Ok(Self::CilDasHexNu),
            "9" => Ok(Self::EchiFlatwormMtDna),
            "10" => Ok(Self::EuplotidNu),
            "11" => Ok(Self::BacArchPlantPlast),
            "12" => Ok(Self::AltYeastNu),
            "13" => Ok(Self::AsciMtDna),
            "14" => Ok(Self::AltFlatwormMtDna),
            "16" => Ok(Self::ChlorMtDna),
            "21" => Ok(Self::TrematodeMtDna),
            "22" => Ok(Self::ScenedesmusMtDna),
            "23" => Ok(Self::ThrausMtDna),
            "24" => Ok(Self::RhabdopMtDna),
            "25" => Ok(Self::CaDivSR1GraciBac),
            "26" => Ok(Self::PachyNu),
            "29" => Ok(Self::MesodiniumNu),
            "30" => Ok(Self::PeritrichNu),
            "33" => Ok(Self::CephalodiscidaeMtDna),
            _ => Err(format!("{} is not a valid genetic code", s)),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_parsing_input_fmt() {
        let file = Path::new("tests/files/simple.nex");
        let input_fmt = infer_input_auto(file);
        assert_eq!(InputFmt::Nexus, input_fmt);
    }

    #[test]
    fn test_summary_mode_parse_and_default() {
        assert_eq!(SummaryMode::default(), SummaryMode::Complete);
        assert_eq!(
            "minimal".parse::<SummaryMode>().unwrap(),
            SummaryMode::Minimal
        );
        assert_eq!("basic".parse::<SummaryMode>().unwrap(), SummaryMode::Basic);
        // Legacy alias for Basic mode.
        assert_eq!(
            "default".parse::<SummaryMode>().unwrap(),
            SummaryMode::Basic
        );
        assert_eq!(
            "complete".parse::<SummaryMode>().unwrap(),
            SummaryMode::Complete
        );
        assert_eq!(SummaryMode::Basic.to_string(), "basic");
        assert!("unknown".parse::<SummaryMode>().is_err());
    }
}
