//! Check input and output sequences
use std::path::Path;

use crate::helper::types;
use crate::helper::types::{DataType, Header, InputFmt, SeqMatrix};
use crate::parser::fasta::Fasta;
use crate::parser::nexus::Nexus;
use crate::parser::phylip::Phylip;

macro_rules! parse_sequence {
    ($self:ident, $format:ident) => {{
        let mut seq = $format::new($self.file, $self.datatype);
        seq.parse();
        (seq.matrix, seq.header)
    }};
}

/// Parse sequence files.
pub struct SeqParser<'a> {
    /// Path to the sequence file.
    file: &'a Path,
    /// The data type of the sequences.
    datatype: &'a DataType,
}

impl<'a> SeqParser<'a> {
    /// Create a new `SeqParser` instance.
    pub fn new(file: &'a Path, datatype: &'a DataType) -> Self {
        Self { file, datatype }
    }

    /// Parse sequence based on the input format and check if the sequences are aligned.
    /// Return a tuple of the sequence matrix and the header.
    ///
    /// # Example
    /// ```
    /// use std::path::Path;
    /// use segul_core::helper::types::{DataType, InputFmt};
    /// use segul_core::helper::sequence::SeqParser;
    ///
    /// let file = Path::new("tests/files/simple.fas");
    /// let datatype = &DataType::Dna;
    /// let input_fmt = &InputFmt::Fasta;
    ///
    /// let seq = SeqParser::new(&file, datatype);
    /// let (matrix, header) = seq.get_alignment(input_fmt);
    /// assert_eq!(matrix.len(), 2);
    /// assert_eq!(header.aligned, true);
    /// ```
    pub fn get_alignment(&self, input_fmt: &'a InputFmt) -> (SeqMatrix, Header) {
        let (matrix, header) = self.parse(input_fmt);

        assert!(
            header.aligned,
            "Found an invalid alignment file. \
            {} is not an alignment. \
            SEGUL assumes the sequences are aligned if they are the same length.",
            self.file.display()
        );

        (matrix, header)
    }

    /// Parse sequence based on the input format.
    /// Similar to `get_alignment` but does not check if the sequences are aligned.
    pub fn parse(&self, input_fmt: &'a InputFmt) -> (SeqMatrix, Header) {
        match input_fmt {
            InputFmt::Fasta => parse_sequence!(self, Fasta),
            InputFmt::Nexus => parse_sequence!(self, Nexus),
            InputFmt::Phylip => parse_sequence!(self, Phylip),
            InputFmt::Auto => {
                let input_fmt = types::infer_input_auto(self.file);
                self.parse(&input_fmt)
            }
        }
    }
}

/// Get a sequence shortest and longest length and check if it is aligned.
/// # Example
/// ```
/// use std::path::Path;
/// use segul_core::helper::types::{DataType, InputFmt};
/// use segul_core::helper::sequence::{SeqParser, SeqCheck};
///
/// let file = Path::new("tests/files/simple.fas");
/// let datatype = &DataType::Dna;
/// let input_fmt = &InputFmt::Fasta;
///
/// let seq = SeqParser::new(&file, datatype);
/// let (matrix, _) = seq.parse(input_fmt);
///
/// let mut seq_check = SeqCheck::new();
/// seq_check.check(&matrix);
/// assert_eq!(seq_check.shortest, 6);
/// assert_eq!(seq_check.longest, 6);
/// assert_eq!(seq_check.is_alignment, true);
/// ```
pub struct SeqCheck {
    /// The shortest sequence length.
    pub shortest: usize,
    /// The longest sequence length.
    pub longest: usize,
    /// A boolean indicating if the sequences are aligned.
    pub is_alignment: bool,
}

impl Default for SeqCheck {
    fn default() -> Self {
        Self::new()
    }
}

impl SeqCheck {
    /// Create a new `SeqCheck` instance.
    pub fn new() -> Self {
        Self {
            shortest: 0,
            longest: 0,
            is_alignment: false,
        }
    }

    /// Check if the sequences are aligned.
    pub fn check(&mut self, matrix: &SeqMatrix) {
        let mut iter = matrix.values();
        let first = iter
            .next()
            .expect("Failed getting the shortest sequence length")
            .len();
        let mut min_len = first;
        let mut max_len = first;
        for seq in iter {
            let len = seq.len();
            min_len = min_len.min(len);
            max_len = max_len.max(len);
        }
        self.shortest = min_len;
        self.longest = max_len;
        self.is_alignment = min_len == max_len;
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_alignment_simple() {
        let file = Path::new("tests/files/simple.nex");
        let datatype = DataType::Dna;
        let input_fmt = InputFmt::Nexus;
        let aln = SeqParser::new(file, &datatype);
        let (matrix, header) = aln.get_alignment(&input_fmt);
        assert_eq!(1, header.ntax);
        assert_eq!(6, header.nchar);
        assert_eq!(1, matrix.len());
    }
}
