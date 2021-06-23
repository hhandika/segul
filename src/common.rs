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
            .expect("CANNOT GET THE SHORTEST SEQUENCE LENGTH");
    }

    fn get_longest_seq_len(&mut self, matrix: &IndexMap<String, String>) {
        self.longest = matrix
            .values()
            .map(|s| s.len())
            .max()
            .expect("CANNOT GET THE LONGEST SEQUENCE LENGTH");
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
