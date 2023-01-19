//! Data types for supported formats
use ahash::AHashMap as HashMap;
use indexmap::IndexMap;

/// Data types for raw read sequences
#[derive(Debug, PartialEq)]
pub enum RawReadFmt {
    /// Infer format from file extension
    Auto,
    /// Fastq format
    Fastq,
    /// Gzip compressed fastq format
    Gzip,
}

/// Data types for input sequences
#[derive(Debug, PartialEq)]
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

/// Data types for output sequences
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

/// Data types for alignment partitions
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

/// Data types for sequence data
pub type SeqMatrix = IndexMap<String, String>;

/// Data types for partition data
pub struct Partition {
    /// Gene/locus name
    pub gene: String,
    /// Start position
    pub start: usize,
    /// End position
    pub end: usize,
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
        }
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

pub enum SummaryMode {
    /// Only write/print minimal summary information
    Minimal,
    /// Print write/print all essential summary information
    Default,
    /// Print all summary information
    Complete,
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
pub enum GeneticCodes {
    StandardCode,         // Ncbi Table 1
    VertMtDna,            // Ncbi Table 2
    YeastMtDna,           // Ncbi Table 3
    MoldProtCoelMtDna,    // Ncbi Table 4
    InvertMtDna,          // Ncbi Table 5
    CilDasHexNu,          // Ncbi Table 6
    EchiFlatwormMtDna,    // Ncbi Table 9
    EuplotidNu,           // Ncbi Table 10
    BacArchPlantPlast,    // Ncbi Table 11
    AltYeastNu,           // Ncbi Table 12
    AsciMtDna,            // Ncbi Table 13
    AltFlatwormMtDna,     // Ncbi Table 14
    ChlorMtDna,           // Ncbi Table 16
    TrematodeMtDna,       // Ncbi Table 21
    ScenedesmusMtDna,     // Ncbi Table 22
    ThrausMtDna,          // Ncbi Table 23
    RhabdopMtDna,         // Ncbi Table 24
    CaDivSR1GraciBac,     // Ncbi Table 25
    PachyNu,              // Ncbi Table 26
    MesodiniumNu,         // Ncbi Table 29
    PeritrichNu,          // Ncbi Table 30
    CephalodiscidaeMtDna, // Ncbi Table 33
}
