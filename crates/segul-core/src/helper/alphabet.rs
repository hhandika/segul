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
/// use segul_core::helper::alphabet;
/// use segul_core::helper::types::DataType;
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
/// use segul_core::helper::alphabet;
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
/// use segul_core::helper::alphabet;
///
/// let aa = String::from("ARNDCQEGH");
/// assert_eq!(true, alphabet::is_valid_aa(&aa));
/// ```
// Public for library only. Do not use directly by the Cli app.
pub fn is_valid_aa(aa: &str) -> bool {
    aa.bytes().all(|char| AA.contains(&char))
}

/// Complement table for DNA, including IUPAC ambiguity codes.
/// Case is preserved. Any other byte, such as gaps and missing data, maps to itself.
const DNA_COMPLEMENT: [u8; 256] = build_dna_complement();

const fn build_dna_complement() -> [u8; 256] {
    let mut table = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = i as u8;
        i += 1;
    }
    let pairs: [(u8, u8); 11] = [
        (b'A', b'T'),
        (b'T', b'A'),
        (b'U', b'A'),
        (b'C', b'G'),
        (b'G', b'C'),
        (b'R', b'Y'),
        (b'Y', b'R'),
        (b'K', b'M'),
        (b'M', b'K'),
        (b'B', b'V'),
        (b'V', b'B'),
    ];
    let mut j = 0;
    while j < pairs.len() {
        let (from, to) = pairs[j];
        table[from as usize] = to;
        table[from.to_ascii_lowercase() as usize] = to.to_ascii_lowercase();
        j += 1;
    }
    // D and H complement each other. S, W, and N are their own complement.
    table[b'D' as usize] = b'H';
    table[b'H' as usize] = b'D';
    table[b'd' as usize] = b'h';
    table[b'h' as usize] = b'd';
    table
}

/// Complement a single DNA character. Case is preserved.
/// Gaps, missing data, and unknown characters are returned unchanged.
pub fn complement_dna(base: u8) -> u8 {
    DNA_COMPLEMENT[base as usize]
}

/// Reverse complement a DNA sequence, including IUPAC ambiguity codes.
/// Case is preserved. Gaps and missing data are kept as they are.
/// # Example
/// ```
/// use segul_core::helper::alphabet;
///
/// assert_eq!(alphabet::reverse_complement_dna(b"ACgt-N"), b"N-acGT".to_vec());
/// ```
pub fn reverse_complement_dna(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&base| complement_dna(base)).collect()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_reverse_complement_dna() {
        assert_eq!(reverse_complement_dna(b"ACGTU"), b"AACGT".to_vec());
        assert_eq!(reverse_complement_dna(b"acgt"), b"acgt".to_vec());
        assert_eq!(reverse_complement_dna(b"RYKMBVDH"), b"DHBVKMRY".to_vec());
        assert_eq!(reverse_complement_dna(b"SWN"), b"NWS".to_vec());
        assert_eq!(reverse_complement_dna(b"A-?.*"), b"*.?-T".to_vec());
    }

    #[test]
    fn test_check_valid_dna() {
        let dna = String::from("agtc?-");
        assert!(is_valid_dna(&dna));
    }

    #[test]
    fn test_check_invalid_dna() {
        let dna = String::from("agtc?-z");
        assert!(!is_valid_dna(&dna));
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
