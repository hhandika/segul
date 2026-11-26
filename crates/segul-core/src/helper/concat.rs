//! Concatenate multiple alignments into a single alignment.
use std::ffi::OsStr;
use std::path::{Path, PathBuf};

use indexmap::IndexSet;
use indicatif::ProgressBar;
use rayon::prelude::*;

use crate::helper::finder::IDs;
use crate::helper::sequence::SeqParser;
use crate::helper::types::{DataType, Header, InputFmt, Partition, SeqMatrix};

use super::types::{OutputFmt, PartitionFmt};

/// Concatenation parameters.
/// Used by features that support concatenation,
/// such as `AlignmentFiltering` and `SequenceFiltering`.
pub struct ConcatParams {
    pub output_fmt: OutputFmt,
    pub part_fmt: PartitionFmt,
    pub prefix: PathBuf,
}

impl ConcatParams {
    pub fn new(output_fmt: &OutputFmt, part_fmt: &PartitionFmt, prefix: &Path) -> Self {
        Self {
            prefix: prefix.to_path_buf(),
            part_fmt: *part_fmt,
            output_fmt: *output_fmt,
        }
    }
}

/// Concatenate multiple alignments into a single alignment.
/// # Arguments
/// * `files` - A mutable slice of PathBuf that holds the alignment files.
/// * `input_fmt` - The alignment input format.
/// * `datatype` - The input data type.
/// # Example
/// ```
/// use std::path::PathBuf;
/// use segul_core::helper::concat::Concat;
/// use segul_core::helper::types::{DataType, InputFmt};
///
/// let mut files = vec![
///     PathBuf::from("tests/files/alignments/gene_1.nex"),
///    PathBuf::from("tests/files/alignments/gene_2.nex"),
/// ];
/// let spinner = segul_core::helper::utils::set_spinner();
/// let input_fmt = InputFmt::Nexus;
/// let datatype = DataType::Dna;
/// let mut concat = Concat::new(&mut files, &input_fmt, &datatype);
/// concat.concat_alignment(&spinner);
///
/// // With no spinner
/// concat.concat_alignment_no_spinner();
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
    /// Create a new Concat instance.
    pub fn new(files: &'a mut [PathBuf], input_fmt: &'a InputFmt, datatype: &'a DataType) -> Self {
        Self {
            input_fmt,
            datatype,
            alignment: SeqMatrix::new(),
            header: Header::new(),
            partition: Vec::new(),
            files,
        }
    }

    /// Concatenate alignments, required a spinner for stdout.
    pub fn concat_alignment(&mut self, spin: &ProgressBar) {
        alphanumeric_sort::sort_path_slice(self.files);
        spin.set_message("Indexing alignments...");
        let (id, file_nchars) =
            IDs::new(self.files, self.input_fmt, self.datatype).id_unique_with_nchars();
        spin.set_message("Concatenating alignments...");
        self.concat(id, &file_nchars);
        self.header.ntax = self.alignment.len();
        self.match_header_datatype();
    }

    /// Concatenate alignments, without a spinner.
    pub fn concat_alignment_no_spinner(&mut self) {
        alphanumeric_sort::sort_path_slice(self.files);
        let (id, file_nchars) =
            IDs::new(self.files, self.input_fmt, self.datatype).id_unique_with_nchars();
        self.concat(id, &file_nchars);
        self.header.ntax = self.alignment.len();
        self.match_header_datatype();
    }

    fn concat(&mut self, id: IndexSet<String>, file_nchars: &[usize]) {
        let ntax = id.len();
        let total_nchar: usize = file_nchars.iter().sum();
        let total_bytes = ntax
            .checked_mul(total_nchar)
            .expect("Concatenated alignment size overflow");
        // Pre-fill with `?` once so missing taxa require no per-locus filling.
        let mut data = vec![b'?'; total_bytes];
        let mut col_offset = 0usize;
        let mut gene_start = 1usize;
        let mut partition = Vec::with_capacity(self.files.len());
        // Parse a small batch of files in parallel, then copy them in file order.
        // Batching bounds how many parsed alignments are held in memory at once.
        let batch_size = rayon::current_num_threads().max(1) * 2;
        for (files, expected_nchars) in self
            .files
            .chunks(batch_size)
            .zip(file_nchars.chunks(batch_size))
        {
            let parsed: Vec<(SeqMatrix, Header)> = files
                .par_iter()
                .map(|file| self.get_alignment(file))
                .collect();
            for ((file, &expected_nchar), (matrix, header)) in
                files.iter().zip(expected_nchars).zip(parsed)
            {
                assert_eq!(
                    header.nchar,
                    expected_nchar,
                    "Sequence length mismatch in {}",
                    file.display()
                );
                let locus_end = col_offset + header.nchar;
                let gene_name = self.parse_aln_name(file);
                let part = self.get_partition(&gene_name, gene_start, locus_end);
                partition.push(part);
                gene_start = locus_end + 1;
                for (taxon_id, seq) in matrix.iter() {
                    if let Some(tax_idx) = id.get_index_of(taxon_id) {
                        let dst_start = tax_idx * total_nchar + col_offset;
                        data[dst_start..dst_start + header.nchar].copy_from_slice(seq);
                    }
                }
                col_offset = locus_end;
            }
        }

        self.alignment = SeqMatrix::from_uniform(id, data, total_nchar);
        self.header.nchar = total_nchar;
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
        Partition::from_range(gene_name, start, end)
    }

    fn parse_aln_name(&self, file: &Path) -> String {
        file.file_stem()
            .and_then(OsStr::to_str)
            .expect("Failed getting alignment name from the file")
            .to_string()
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
    use crate::helper::finder::SeqFileFinder;
    use crate::helper::utils;

    const DNA: DataType = DataType::Dna;

    #[test]
    fn test_concat_nexus() {
        let path = Path::new("tests/files/alignments/");
        let mut files = SeqFileFinder::new(path).find(&InputFmt::Nexus);
        let mut concat = Concat::new(&mut files, &InputFmt::Nexus, &DNA);
        let spin = utils::set_spinner();
        concat.concat_alignment(&spin);
        assert_eq!(3, concat.alignment.len());
    }

    #[test]
    #[should_panic]
    fn test_get_alignment_panic() {
        let path = Path::new("tests/files/alignments/");
        let mut files = SeqFileFinder::new(path).find(&InputFmt::Nexus);
        let concat = Concat::new(&mut files, &InputFmt::Nexus, &DNA);
        concat.get_alignment(Path::new("."));
    }

    #[test]
    fn test_concat_check_result() {
        let path = Path::new("tests/files/alignments/");
        let mut files = SeqFileFinder::new(path).find(&InputFmt::Nexus);
        let mut concat = Concat::new(&mut files, &InputFmt::Nexus, &DNA);
        let spin = utils::set_spinner();
        concat.concat_alignment(&spin);
        let abce = concat.alignment.get("ABCE").unwrap();
        let res = b"??????????????gatattagtata";
        assert_eq!(res.as_slice(), abce);
    }

    #[test]
    fn test_concat_partition() {
        let path = Path::new("tests/files/alignments/");
        let mut files = SeqFileFinder::new(path).find(&InputFmt::Nexus);
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
    fn test_header_datatype() {
        let path = Path::new("tests/files/alignments/");
        let mut files = SeqFileFinder::new(path).find(&InputFmt::Nexus);
        let mut concat = Concat::new(&mut files, &InputFmt::Nexus, &DataType::Aa);
        concat.match_header_datatype();
        assert_eq!(concat.header.datatype, String::from("protein"));
    }
}
