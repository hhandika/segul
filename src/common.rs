use std::path::Path;

use indexmap::IndexMap;

pub enum SeqFormat {
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
// We only check lowercase since all the parser
// change the case to lowercase.
fn is_valid_dna(dna: &str) -> bool {
    let valid_chars = String::from("acgtryswkmbdhvn.-?");
    dna.chars().all(|char| valid_chars.contains(char))
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
