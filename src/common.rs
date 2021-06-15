use indexmap::IndexMap;

pub enum OutputFormat {
    Fasta,
    Nexus,
    Phylip,
}

pub enum PartitionFormat {
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

pub trait SeqCheck {
    fn check_is_alignment(&self, matrix: &IndexMap<String, String>) -> bool {
        let shortest = self.get_shortest_seq_len(matrix);
        let longest = self.get_longest_seq_len(matrix);
        longest == shortest
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
