//! Provide utilities to generate and checking DNA and Amino Acid alphabet.
use std::path::Path;

use crate::helper::types::DataType;

/// Alphabeth for dna. In uppercase only.
pub const DNA_STR_UPPERCASE: &str = "?-ACGTNRYSWKMBDHV.";
/// Alphabeth for amino acid. In uppercase only.
pub const AA_STR_UPPERCASE: &str = "?-ARNDCQEGHILKMFPSTWYVYXBZJU*.~";

// Alphabeth for dna.

/// Alphabeth for dna. All cases stored as bytes.
/// Include IUPAC characters plus ambiguous, missing, and gap characters (?, -, *, etc.)
/// Source: http://www.iqtree.org/doc/Frequently-Asked-Questions
const DNA: &[u8] = b"?-ACGTRYSWKMBDHVNacgtryswkmbdhvn.";

// Alphabeth for amino acid.
// Include 20 IUPAC characters,
// ambiguous, missing, and gap characters (X,?,-,.,~, *, etc.)
// Source: http://www.iqtree.org/doc/Frequently-Asked-Questions
/// Alphabeth for amino acid. All cases stored as bytes.
/// Include 20 IUPAC characters,
/// ambiguous, missing, and gap characters (X,?,-,.,~, *, etc.)
/// Source: http://www.iqtree.org/doc/Frequently-Asked-Questions
const AA: &[u8] = b"?-ARNDCQEGHILKMFPSTWYVYXBZJU*.~";

/// Check for valid DNA or amino acid sequence.
/// Panic if the sequence is invalid.
/// # Arguments
/// * `input` - A Path object that holds the path to the input file.
/// * `datatype` - A DataType object that holds the type of the sequence.
/// * `id` - A string slice that holds the sequence id.
///
/// # Example
/// We provide invalid DNA sequence to check_valid_seq function.
/// The function should panic.
/// The panic message should be:
/// "Ups... The sequence seq_1 in file path/to/file is not a dna sequence.
/// Check whether the sequence is amino acid"
/// ```should_panic
/// use std::path::Path;
/// use segul::helper::alphabet;
/// use segul::helper::types::DataType;
///
/// let sample = Path::new("path/to/file");
/// let datatype = DataType::Dna;
/// let id = "seq_1";
/// let seq = String::from("agtc?)-"); // invalid dna sequence
/// alphabet::check_valid_seq(sample, &datatype, id, &seq);
///```
pub fn check_valid_seq(input: &Path, datatype: &DataType, id: &str, seq: &str) {
    match datatype {
        DataType::Dna => check_valid_dna(input, id, seq),
        DataType::Aa => check_valid_aa(input, id, seq),
        DataType::Ignore => (),
    }
}

/// Check for valid DNA sequence only.
/// Similar behavior as check_valid_seq().
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

fn check_valid_aa(input: &Path, id: &str, aa: &str) {
    if !is_valid_aa(aa) {
        panic!(
            "Ups... The sequence {} in file {} is not an amino acid sequence.",
            id,
            input.display()
        );
    }
}

/// Check for valid DNA sequence.
/// Return true if the sequence is valid.
/// Return false if the sequence is not valid.
/// # Arguments
/// * `dna` - A string slice that holds the DNA sequence.
/// # Example
/// ```
/// use std::path::Path;
/// use segul::helper::alphabet;
///
///
/// let dna = String::from("agtc?-");
/// assert_eq!(true, alphabet::is_valid_dna(&dna));
/// ```
// Public for library only. Do not use directly by the Cli app.
pub fn is_valid_dna(dna: &str) -> bool {
    dna.bytes().all(|char| DNA.contains(&char))
}

/// Check for valid amino acid sequence.
/// Return true if the sequence is valid.
/// Return false if the sequence is not valid.
/// # Arguments
/// * `aa` - A string slice that holds the amino acid sequence.
/// # Example
/// ```
/// use std::path::Path;
/// use segul::helper::alphabet;
///
/// let aa = String::from("ARNDCQEGH");
/// assert_eq!(true, alphabet::is_valid_aa(&aa));
/// ```
// Public for library only. Do not use directly by the Cli app.
pub fn is_valid_aa(aa: &str) -> bool {
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
