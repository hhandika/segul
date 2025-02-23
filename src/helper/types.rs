//! Data types for supported formats
use std::ffi::OsStr;
use std::path::Path;

use ahash::AHashMap as HashMap;
use indexmap::IndexMap;

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
/// use segul::helper::types::InputFmt;
/// use segul::helper::types::infer_input_auto;
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

pub enum DnaStrand {
    Forward,
    Reverse,
}

impl DnaStrand {
    pub fn from_char(c: char) -> Self {
        match c {
            '+' => DnaStrand::Forward,
            '-' => DnaStrand::Reverse,
            _ => panic!("Invalid DNA strand"),
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

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum SummaryMode {
    /// Only write/print minimal summary information
    Minimal,
    /// Write all essential summary information
    Default,
    /// Write all summary information
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
            Self::MoldProtCoelMtDna => write!(f, "Mold, Protozoan, and Coelenterate Mitochondrial DNA and the Mycoplasma/Spiroplasma Code"),
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
}
