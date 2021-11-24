use indexmap::IndexMap;

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
    None,
}

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

pub enum NCBITable {
    StandardCode,
    MtDna,
}
