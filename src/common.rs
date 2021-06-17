use std::path::Path;

use indexmap::IndexMap;

pub enum InputFormat {
    Nexus,
    Phylip,
    Fasta,
}

pub enum OutputFormat {
    Fasta,
    Nexus,
    Phylip,
}

pub enum PartitionFormat {
    Nexus,
    NexusSeparate,
    Raxml,
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

pub trait SeqCheck {
    fn check_is_alignment(&self, shortest: &usize, longest: &usize) -> bool {
        shortest == longest
    }

    fn get_sequence_len(&self, matrix: &IndexMap<String, String>) -> (usize, usize) {
        let shortest = self.get_shortest_seq_len(matrix);
        let longest = self.get_longest_seq_len(matrix);
        (shortest, longest)
    }

    fn get_shortest_seq_len(&self, matrix: &IndexMap<String, String>) -> usize {
        matrix
            .values()
            .min_by_key(|seq| seq.len())
            .expect("CANNOT GET LONGEST ALIGNMENT LEN")
            .len()
    }

    fn get_longest_seq_len(&self, matrix: &IndexMap<String, String>) -> usize {
        matrix
            .values()
            .max_by_key(|seq| seq.len())
            .expect("CANNOT GET SHORTEST ALIGNMENT LEN")
            .len()
    }
}

pub fn check_valid_dna(input: &Path, id: &str, dna: &str) {
    if !is_valid_dna(dna) {
        panic!(
            "INVALID DNA SEQUENCE FOUND FOR {} IN FILE {}",
            id,
            input.display()
        );
    }
}

// Alphabeth for dna.
// Include IUPAC characters plus missing symbol (?)
fn is_valid_dna(dna: &str) -> bool {
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

    #[test]
    #[should_panic]
    fn check_invalid_dna_panic_test() {
        let sample = Path::new(".");
        let id = "ABCD";
        let dna = String::from("AGTC?-Z");
        check_valid_dna(sample, id, &dna);
    }
}
