//! Concatenate multiple alignments into a single alignment.
use std::ffi::OsStr;
use std::path::{Path, PathBuf};

use indexmap::{IndexMap, IndexSet};
use indicatif::ProgressBar;

use crate::helper::finder::IDs;
use crate::helper::sequence::SeqParser;
use crate::helper::types::{DataType, Header, InputFmt, Partition, SeqMatrix};

/// Concatenate multiple alignments into a single alignment.
/// # Arguments
/// * `files` - A mutable slice of PathBuf that holds the alignment files.
/// * `input_fmt` - The alignment input format.
/// * `datatype` - The input data type.
/// # Example
/// ```
/// use std::path::PathBuf;
/// use segul::helper::concat::Concat;
/// use segul::helper::types::{DataType, InputFmt};
///
/// let mut files = vec![
///     PathBuf::from("tests/files/concat/gene_1.nex"),
///    PathBuf::from("tests/files/concat/gene_2.nex"),
/// ];
/// let spinner = segul::helper::utils::set_spinner();
/// let input_fmt = InputFmt::Nexus;
/// let datatype = DataType::Dna;
/// let mut concat = Concat::new(&mut files, &input_fmt, &datatype);
/// concat.concat_alignment(&spinner);
/// ```
pub struct Concat<'a> {
    /// The concatenated alignment.
    pub alignment: SeqMatrix,
    /// The header of the concatenated alignment.
    pub header: Header,
    /// The partitions of the concatenated alignment.
    pub partition: Vec<Partition>,
    /// The input data type.
    datatype: &'a DataType,
    /// The alignment input format.
    input_fmt: &'a InputFmt,
    /// The alignment files.
    files: &'a mut [PathBuf],
}

impl<'a> Concat<'a> {
    pub fn new(files: &'a mut [PathBuf], input_fmt: &'a InputFmt, datatype: &'a DataType) -> Self {
        Self {
            input_fmt,
            datatype,
            alignment: IndexMap::new(),
            header: Header::new(),
            partition: Vec::new(),
            files,
        }
    }

    /// Concatenate alignments, spinner included for CLI.
    pub fn concat_alignment(&mut self, spin: &ProgressBar) {
        alphanumeric_sort::sort_path_slice(self.files);
        spin.set_message("Indexing alignments...");
        let id = IDs::new(self.files, self.input_fmt, self.datatype).id_unique();
        spin.set_message("Concatenating alignments...");
        self.concat(&id);
        self.header.ntax = self.alignment.len();
        self.match_header_datatype();
    }

    /// Concatenate alignments, no spinner for GUI.
    pub fn concat_alignment_no_spinner(&mut self) {
        alphanumeric_sort::sort_path_slice(self.files);
        let id = IDs::new(self.files, self.input_fmt, self.datatype).id_unique();
        self.concat(&id);
        self.header.ntax = self.alignment.len();
        self.match_header_datatype();
    }

    fn concat(&mut self, id: &IndexSet<String>) {
        let mut alignment = IndexMap::with_capacity(id.len());
        let mut nchar = 0;
        let mut gene_start = 1;
        let mut partition = Vec::new();
        self.files.iter().for_each(|file| {
            let (matrix, header) = self.get_alignment(file);
            nchar += header.nchar; // increment sequence length using the value from parser
            let gene_name = self.parse_aln_name(file);
            let part = self.get_partition(&gene_name, gene_start, nchar);
            partition.push(part);
            gene_start = nchar + 1;
            id.iter().for_each(|id| match matrix.get(id) {
                Some(seq) => {
                    self.insert_alignment(&mut alignment, id, seq);
                }
                None => {
                    let seq = self.get_missings(header.nchar);
                    self.insert_alignment(&mut alignment, id, &seq);
                }
            });
        });

        self.alignment = alignment;
        self.header.nchar = nchar;
        self.partition = partition;
    }

    fn get_alignment(&self, file: &Path) -> (SeqMatrix, Header) {
        let aln = SeqParser::new(file, self.datatype);
        let (matrix, header) = aln.get_alignment(self.input_fmt);
        assert!(
            header.ntax != 0,
            "Found an empty alignment {}",
            file.display()
        );
        (matrix, header)
    }

    fn get_partition(&self, gene_name: &str, start: usize, end: usize) -> Partition {
        let mut part = Partition::new();
        part.gene = gene_name.to_string();
        part.start = start;
        part.end = end;
        part
    }

    fn insert_alignment(&self, alignment: &mut SeqMatrix, id: &str, seq: &str) {
        match alignment.get_mut(id) {
            Some(seqs) => seqs.push_str(seq),
            None => {
                alignment.insert(id.to_string(), seq.to_string());
            }
        }
    }

    fn parse_aln_name(&self, file: &Path) -> String {
        file.file_stem()
            .and_then(OsStr::to_str)
            .expect("Failed getting alignment name from the file")
            .to_string()
    }

    #[inline]
    fn get_missings(&self, len: usize) -> String {
        "?".repeat(len)
    }

    #[inline]
    fn match_header_datatype(&mut self) {
        if let DataType::Aa = self.datatype {
            self.header.datatype = String::from("protein")
        };
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::helper::finder::Files;
    use crate::helper::utils;

    const DNA: DataType = DataType::Dna;

    #[test]
    fn test_concat_nexus() {
        let path = Path::new("tests/files/concat/");
        let mut files = Files::new(path, &InputFmt::Nexus).find();
        let mut concat = Concat::new(&mut files, &InputFmt::Nexus, &DNA);
        let spin = utils::set_spinner();
        concat.concat_alignment(&spin);
        assert_eq!(3, concat.alignment.len());
    }

    #[test]
    #[should_panic]
    fn test_get_alignment_panic() {
        let path = Path::new("tests/files/concat/");
        let mut files = Files::new(path, &InputFmt::Nexus).find();
        let concat = Concat::new(&mut files, &InputFmt::Nexus, &DNA);
        concat.get_alignment(Path::new("."));
    }

    #[test]
    fn test_concat_check_result() {
        let path = Path::new("tests/files/concat/");
        let mut files = Files::new(path, &InputFmt::Nexus).find();
        let mut concat = Concat::new(&mut files, &InputFmt::Nexus, &DNA);
        let spin = utils::set_spinner();
        concat.concat_alignment(&spin);
        let abce = concat.alignment.get("ABCE").unwrap();
        let res = "??????????????gatattagtata";
        assert_eq!(res, abce);
    }

    #[test]
    fn test_concat_partition() {
        let path = Path::new("tests/files/concat/");
        let mut files = Files::new(path, &InputFmt::Nexus).find();
        let mut concat = Concat::new(&mut files, &InputFmt::Nexus, &DNA);
        let spin = utils::set_spinner();
        concat.concat_alignment(&spin);
        assert_eq!(1, concat.partition[0].start);
        assert_eq!(6, concat.partition[0].end);
        assert_eq!(7, concat.partition[1].start);
        assert_eq!(14, concat.partition[1].end);
        assert_eq!(15, concat.partition[2].start);
        assert_eq!(20, concat.partition[2].end);
        assert_eq!(21, concat.partition[3].start);
        assert_eq!(26, concat.partition[3].end);
    }

    #[test]
    fn test_get_gaps() {
        let len = 5;
        let gaps = "?????";
        assert_eq!(
            gaps,
            Concat::new(&mut [PathBuf::from(".")], &InputFmt::Fasta, &DNA).get_missings(len)
        )
    }

    #[test]
    fn test_header_datatype() {
        let path = Path::new("tests/files/concat/");
        let mut files = Files::new(path, &InputFmt::Nexus).find();
        let mut concat = Concat::new(&mut files, &InputFmt::Nexus, &DataType::Aa);
        concat.match_header_datatype();
        assert_eq!(concat.header.datatype, String::from("protein"));
    }
}
