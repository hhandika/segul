use std::ffi::OsStr;
use std::path::Path;

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

#[derive(Clone)]
pub enum DataType {
    Dna,
    Aa,
    Ignore,
}

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
}

impl Header {
    pub fn new() -> Self {
        Self {
            ntax: 0,
            nchar: 0,
            datatype: String::from("dna"),
            missing: '?',
            gap: '-',
        }
    }
}

pub fn infer_input_auto(input: &Path) -> InputFmt {
    let ext: &str = input
        .extension()
        .and_then(OsStr::to_str)
        .expect("Failed parsing extension");
    match ext {
        "fas" | "fa" | "fasta" => InputFmt::Fasta,
        "nex" | "nexus" => InputFmt::Nexus,
        "phy" | "phylip" => InputFmt::Phylip,
        _ => panic!(
            "Ups... The program cannot recognize the file extension. \
        Maybe try specify the input format using -f --format option."
        ),
    }
}

pub struct SeqCheck {
    pub shortest: usize,
    pub longest: usize,
    pub is_alignment: bool,
}

impl SeqCheck {
    pub fn new() -> Self {
        Self {
            shortest: 0,
            longest: 0,
            is_alignment: false,
        }
    }

    pub fn get_sequence_info(&mut self, matrix: &IndexMap<String, String>) {
        self.get_shortest_seq_len(matrix);
        self.get_longest_seq_len(matrix);
        self.check_is_alignment();
    }

    fn check_is_alignment(&mut self) {
        self.is_alignment = self.shortest == self.longest;
    }

    fn get_shortest_seq_len(&mut self, matrix: &IndexMap<String, String>) {
        self.shortest = matrix
            .values()
            .map(|s| s.len())
            .min()
            .expect("Failed getting the shortest failed length");
    }

    fn get_longest_seq_len(&mut self, matrix: &IndexMap<String, String>) {
        self.longest = matrix
            .values()
            .map(|s| s.len())
            .max()
            .expect("Failed getting the longest failed length");
    }
}

pub fn check_valid_dna(input: &Path, id: &str, dna: &str) {
    if !is_valid_dna(dna) {
        panic!(
            "Ups... The {} is not a dna sequence. Found {} in file {}",
            id,
            dna,
            input.display()
        );
    }
}

pub fn check_valid_seq(input: &Path, datatype: &DataType, id: &str, seq: &str) {
    match datatype {
        DataType::Dna => check_valid_dna(input, id, seq),
        DataType::Aa => check_valid_aa(input, id, seq),
        DataType::Ignore => (),
    }
}

fn check_valid_aa(input: &Path, id: &str, aa: &str) {
    if !is_valid_aa(aa) {
        panic!(
            "Ups... The {} is not an amino acid sequence. Found {} in file {}",
            id,
            aa,
            input.display()
        );
    }
}

// Alphabeth for dna.
// Include IUPAC characters plus missing, and gaps characters (?, -, *, etc.)
// Source: http://www.iqtree.org/doc/Frequently-Asked-Questions
fn is_valid_dna(dna: &str) -> bool {
    let valid_dna = b"ACGTRYSWKMBDHVNacgtryswkmbdhvn.-?";
    dna.bytes().all(|char| valid_dna.contains(&char))
}

// Alphabeth for amino acid.
// Include 20 IUPAC characters,
// ambiguous, missing, and gaps characters (X,?,-,.,~, *, etc.)
// Source: http://www.iqtree.org/doc/Frequently-Asked-Questions
fn is_valid_aa(aa: &str) -> bool {
    let valid_aa = b"ARNDCQEGHILKMFPSTWYVYXBZJU?-.~*";
    aa.bytes().all(|char| valid_aa.contains(&char))
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn check_valid_dna_test() {
        let dna = String::from("agtc?-");
        assert_eq!(true, is_valid_dna(&dna));
    }

    #[test]
    fn check_invalid_dna_test() {
        let dna = String::from("agtc?-z");
        assert_eq!(false, is_valid_dna(&dna));
    }

    #[test]
    #[should_panic]
    fn check_invalid_dna_panic_test() {
        let sample = Path::new(".");
        let id = "ABCD";
        let dna = String::from("agta?-z");
        check_valid_dna(sample, id, &dna);
    }
}
