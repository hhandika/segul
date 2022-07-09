use ahash::AHashMap as HashMap;
use indexmap::IndexMap;

#[derive(Debug, PartialEq)]
pub enum InputFmt {
    Auto,
    Fasta,
    Nexus,
    Phylip,
}

pub enum OutputFmt {
    Fasta,
    Nexus,
    Phylip,
    FastaInt,
    NexusInt,
    PhylipInt,
}

pub enum PartitionFmt {
    Charset,
    CharsetCodon,
    Nexus,
    NexusCodon,
    Raxml,
    RaxmlCodon,
}

#[derive(PartialEq)]
pub enum DataType {
    Dna,
    Aa,
    Ignore,
}

pub type SeqMatrix = IndexMap<String, String>;

pub struct Partition {
    pub gene: String,
    pub start: usize,
    pub end: usize,
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

#[derive(Clone, Debug)]
pub struct Header {
    pub ntax: usize,
    pub nchar: usize,
    pub datatype: String,
    pub missing: char,
    pub gap: char,
    pub aligned: bool,
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

pub struct TaxonRecords {
    pub char_counts: HashMap<char, usize>,
    pub locus_counts: usize,
    pub total_chars: usize,
    pub gc_count: usize,
    pub at_count: usize,
    pub nucleotides: usize,
    pub missing_data: usize,
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
