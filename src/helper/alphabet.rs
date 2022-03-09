use std::path::Path;

use crate::helper::types::DataType;

pub const DNA_STR_UPPERCASE: &str = "?-ACGTNRYSWKMBDHV.";
pub const AA_STR_UPPERCASE: &str = "?-ARNDCQEGHILKMFPSTWYVYXBZJU*.~";

// Alphabeth for dna.
// Include IUPAC characters plus ambiguous, missing, and gap characters (?, -, *, etc.)
// Source: http://www.iqtree.org/doc/Frequently-Asked-Questions
const DNA: &[u8] = b"?-ACGTRYSWKMBDHVNacgtryswkmbdhvn.";

// Alphabeth for amino acid.
// Include 20 IUPAC characters,
// ambiguous, missing, and gap characters (X,?,-,.,~, *, etc.)
// Source: http://www.iqtree.org/doc/Frequently-Asked-Questions
const AA: &[u8] = b"?-ARNDCQEGHILKMFPSTWYVYXBZJU*.~";

pub fn check_valid_dna(input: &Path, id: &str, dna: &str) {
    if !is_valid_dna(dna) {
        panic!(
            "Ups... The sequence {} in file {} is not a dna sequence. \
            Check whether the sequence is amino acid",
            id,
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
            "Ups... The sequence {} in file {} is not an amino acid sequence.",
            id,
            input.display()
        );
    }
}

fn is_valid_dna(dna: &str) -> bool {
    dna.bytes().all(|char| DNA.contains(&char))
}

fn is_valid_aa(aa: &str) -> bool {
    aa.bytes().all(|char| AA.contains(&char))
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_check_valid_dna() {
        let dna = String::from("agtc?-");
        assert_eq!(true, is_valid_dna(&dna));
    }

    #[test]
    fn test_check_invalid_dna() {
        let dna = String::from("agtc?-z");
        assert_eq!(false, is_valid_dna(&dna));
    }

    #[test]
    #[should_panic]
    fn test_check_invalid_dna_panic() {
        let sample = Path::new(".");
        let id = "ABCD";
        let dna = String::from("agta?-z");
        check_valid_dna(sample, id, &dna);
    }
}
