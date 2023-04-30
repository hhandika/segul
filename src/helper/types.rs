//! Data types for supported formats
use ahash::AHashMap as HashMap;
use indexmap::IndexMap;

/// Data types for raw read sequences
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum RawReadFmt {
    /// Infer format from file extension
    Auto,
    /// Fastq format
    Fastq,
    /// Gzip compressed fastq format
    Gzip,
}

impl std::fmt::Display for RawReadFmt {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Auto => write!(f, "auto"),
            Self::Fastq => write!(f, "fastq"),
            Self::Gzip => write!(f, "gzip"),
        }
    }
}

impl std::str::FromStr for RawReadFmt {
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

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum SummaryMode {
    /// Only write/print minimal summary information
    Minimal,
    /// Print write/print all essential summary information
    Default,
    /// Print all summary information
    Complete,
}

impl std::fmt::Display for SummaryMode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Minimal => write!(f, "minimal"),
            Self::Default => write!(f, "default"),
            Self::Complete => write!(f, "complete"),
        }
    }
}

impl std::str::FromStr for SummaryMode {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "minimal" => Ok(Self::Minimal),
            "default" => Ok(Self::Default),
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
