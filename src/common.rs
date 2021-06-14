pub enum SeqFormat {
    Fasta,
    Nexus,
    Phylip,
}

pub enum SeqPartition {
    Nexus,
    NexusSeparate,
    Phylip,
    None,
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

pub struct Header {
    pub ntax: Option<usize>,
    pub nchar: Option<usize>,
    pub datatype: Option<String>,
    pub missing: Option<char>,
    pub gap: Option<char>,
}

impl Header {
    pub fn new() -> Self {
        Self {
            ntax: None,
            nchar: None,
            datatype: None,
            missing: None,
            gap: None,
        }
    }
}

// Alphabeth for dna.
// Include IUPAC characters plus missing symbol (?)
pub fn is_valid_dna(dna: &str) -> bool {
    let valid_chars = String::from("ACGTRYSWKMBDHVNacgtryswkmbdhvn.-?");
    dna.chars().all(|char| valid_chars.contains(char))
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn check_valid_dna_test() {
        let dna = String::from("AGTC?-");
        assert_eq!(true, is_valid_dna(&dna));
    }

    #[test]
    fn check_invalid_dna_test() {
        let dna = String::from("AGTC?-Z");
        assert_eq!(false, is_valid_dna(&dna));
    }
}
