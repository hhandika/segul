//! Check input and output sequences

use std::ffi::OsStr;
use std::path::Path;

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

/// Infer input format automatically based on the file extension.
/// Return the input format.
/// # Example
/// ```
/// use std::path::Path;
/// use segul::helper::types::InputFmt;
/// use segul::helper::sequence::infer_input_auto;
///
/// let file = Path::new("tests/files/simple.fas");
/// let input_fmt = infer_input_auto(&file);
/// assert_eq!(input_fmt, InputFmt::Fasta);
/// ```
pub fn infer_input_auto(input: &Path) -> InputFmt {
    let ext: &str = input
        .extension()
        .and_then(OsStr::to_str)
        .expect("Failed parsing extension");
    match ext {
        "fas" | "fa" | "fasta" => InputFmt::Fasta,
        "nex" | "nexus" => InputFmt::Nexus,
        "phy" | "phylip" => InputFmt::Phylip,
        _ => panic!(
            "The program cannot recognize the file extension. \
        Maybe try to specify the input format using the -f or --input-format option."
        ),
    }
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
    /// use segul::helper::types::{DataType, InputFmt};
    /// use segul::helper::sequence::SeqParser;
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
                let input_fmt = infer_input_auto(self.file);
                self.parse(&input_fmt)
            }
        }
    }
}

pub struct SeqCheck {
    pub shortest: usize,
    pub longest: usize,
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
        self.shortest_seq_len(matrix);
        self.longest_seq_len(matrix);
        self.check_is_alignment();
    }

    fn check_is_alignment(&mut self) {
        self.is_alignment = self.shortest == self.longest;
    }

    fn shortest_seq_len(&mut self, matrix: &SeqMatrix) {
        self.shortest = matrix
            .values()
            .map(|s| s.len())
            .min_by(|a, b| a.cmp(b))
            .expect("Failed getting the shortest sequence length");
    }

    fn longest_seq_len(&mut self, matrix: &SeqMatrix) {
        self.longest = matrix
            .values()
            .map(|s| s.len())
            .max_by(|a, b| a.cmp(b))
            .expect("Failed getting the longest sequence length");
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

    #[test]
    fn test_parsing_input_fmt() {
        let file = Path::new("tests/files/simple.nex");
        let input_fmt = infer_input_auto(file);
        assert_eq!(InputFmt::Nexus, input_fmt);
    }
}
