pub enum SeqFormat {
    Fasta,
    Nexus,
    Phylip,
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
