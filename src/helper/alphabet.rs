use std::path::Path;

use crate::helper::common::DataType;

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
// Include IUPAC characters plus ambiguous, missing, and gap characters (?, -, *, etc.)
// Source: http://www.iqtree.org/doc/Frequently-Asked-Questions
fn is_valid_dna(dna: &str) -> bool {
    let valid_dna = b"ACGTRYSWKMBDHVNacgtryswkmbdhvn.-?";
    dna.bytes().all(|char| valid_dna.contains(&char))
}

// Alphabeth for amino acid.
// Include 20 IUPAC characters,
// ambiguous, missing, and gap characters (X,?,-,.,~, *, etc.)
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
