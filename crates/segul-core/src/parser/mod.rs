//! All the parsers for the supported formats.
#[cfg(feature = "maf")]
pub mod bed;
pub mod delimited;
pub mod fasta;
#[cfg(feature = "maf")]
pub mod maf;
pub mod nexus;
pub mod partition;
pub mod phylip;
#[cfg(feature = "read")]
pub mod qscores;
pub mod txt;

#[cfg(target_os = "windows")]
const CAR_RETURN: u8 = b'\r';

#[cfg(feature = "maf")]
const END_OF_LINE: u8 = b'\n';
#[cfg(feature = "maf")]
const EOF: usize = 0;

/// Upper bound for pre-allocating interleaved sequences,
/// guarding against invalid NCHAR values in the header.
pub(crate) const MAX_RESERVE: usize = 1 << 30;

/// Split a matrix line into its ID and sequence.
/// Returns None unless the line has exactly two
/// whitespace-separated tokens.
/// Works on bytes because `split_whitespace` decodes every
/// character of the sequence, which is slow for long sequences.
pub(crate) fn split_id_seq(line: &str) -> Option<(&str, &str)> {
    let line = line.trim();
    let split = line.bytes().position(|b| b.is_ascii_whitespace())?;
    let (id, seq) = (&line[..split], line[split..].trim_start());
    if has_whitespace(seq.as_bytes()) {
        return None;
    }
    Some((id, seq))
}

fn has_whitespace(bytes: &[u8]) -> bool {
    // ASCII whitespace and control characters are all <= b' '.
    // Branch-free per chunk so the compiler can vectorize the scan.
    // Matches are then confirmed with the exact check.
    bytes
        .chunks(64)
        .any(|c| c.iter().fold(false, |acc, &b| acc | (b <= b' ')))
        && bytes.iter().any(|b| b.is_ascii_whitespace())
}

macro_rules! insert_matrix {
    ($self: ident, $id: ident, $seq: ident ) => {
        match $self.matrix.get($id) {
            Some(original_seq) => panic!(
                "Found duplicate IDs for file {}. \
                First duplicate ID found: {}. \
                Both sequences are the same: {}.",
                $self.input.display(),
                $id,
                original_seq == $seq
            ),
            None => {
                $self.matrix.insert($id.to_string(), $seq.to_string());
            }
        }
    };
}

macro_rules! warn_duplicate_ids {
    ($self: ident, $ids: ident) => {
        if $ids.len() != $self.header.ntax {
            log::warn!(
                "{} Found problematic alignment: {}.\n\
                Number of taxa in the matrix ({} taxa) \
                does not match the number of taxa in the header ({} taxa).\n\
                If this is `id` subcommand, try to use the same file as an input \
                for other subcommands.\n\
                If the cause is duplicate IDs, \
                the other subcommands will show the first found duplicate.\n\n",
                "WARNING!".red(),
                $self.input.display(),
                $ids.len(),
                $self.header.ntax
            );
        }
    };
}

pub(crate) use insert_matrix;
pub(crate) use warn_duplicate_ids;

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_split_id_seq() {
        assert_eq!(Some(("ABCD", "ACGT")), split_id_seq("ABCD ACGT"));
        assert_eq!(Some(("ABCD", "ACGT")), split_id_seq(" ABCD \t  ACGT\r\n"));
        assert_eq!(None, split_id_seq("AB CD ACGT"));
        assert_eq!(None, split_id_seq("ABCD ACG T"));
        assert_eq!(None, split_id_seq("ABCD"));
        assert_eq!(None, split_id_seq(""));
    }
}
