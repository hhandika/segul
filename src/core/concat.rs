//! Multi-Sequence Alignment Module.
//! Contains methods for working with multi-sequence alignments.

use std::ffi::OsStr;
use std::path::{Path, PathBuf};

use ansi_term::Colour::Yellow;
use indexmap::{IndexMap, IndexSet};
use indicatif::ProgressBar;

use crate::helper::finder::IDs;
use crate::helper::sequence::Sequence;
use crate::helper::types::{
    DataType, Header, InputFmt, OutputFmt, Partition, PartitionFmt, SeqMatrix,
};
use crate::helper::utils;
use crate::writer::sequences::SeqWriter;

pub struct ConcatHandler<'a> {
    input_fmt: &'a InputFmt,
    output: &'a Path,
    output_fmt: &'a OutputFmt,
    part_fmt: &'a PartitionFmt,
}

impl<'a> ConcatHandler<'a> {
    pub fn new(
        input_fmt: &'a InputFmt,
        output: &'a Path,
        output_fmt: &'a OutputFmt,
        part_fmt: &'a PartitionFmt,
    ) -> Self {
        Self {
            input_fmt,
            output,
            output_fmt,
            part_fmt,
        }
    }

    pub fn concat_alignment(&mut self, files: &mut [PathBuf], datatype: &DataType) {
        let mut concat = Concat::new(files, self.input_fmt, datatype);
        let spin = utils::set_spinner();
        self.write_alignment(&mut concat, &spin);
    }

    fn write_alignment(&mut self, concat: &mut Concat, spin: &ProgressBar) {
        concat.concat_alignment(spin);
        let mut save = SeqWriter::new(
            self.output,
            &concat.alignment,
            &concat.header,
            Some(&concat.partition),
            self.part_fmt,
        );
        let part_fname = self.get_partition_path();
        save.set_partition_name(&part_fname);
        spin.set_message("Writing output files...");
        save.write_sequence(self.output_fmt)
            .expect("Failed writing the output file");
        spin.finish_with_message("Finished concatenating alignments!\n");
        self.print_output_info(concat.partition.len(), &concat.header, &part_fname);
    }

    fn get_partition_path(&mut self) -> PathBuf {
        match self.part_fmt {
            PartitionFmt::Charset | PartitionFmt::CharsetCodon => {
                PathBuf::from("charset (in-file)")
            }
            PartitionFmt::Nexus | PartitionFmt::NexusCodon => self.get_part_fname("nex"),
            PartitionFmt::Raxml | PartitionFmt::RaxmlCodon => self.get_part_fname("txt"),
            _ => unreachable!("Please, define a valid partition format!"),
        }
    }

    fn get_part_fname(&self, ext: &str) -> PathBuf {
        let fname = format!(
            "{}_partition.{}",
            self.output
                .file_stem()
                .expect("Failed getting file name for partition")
                .to_string_lossy(),
            ext
        );

        self.output
            .parent()
            .expect("Failed getting output parent directory")
            .join(Path::new(&fname))
    }

    fn print_output_info(&self, count: usize, header: &Header, part_file: &Path) {
        log::info!("{}", Yellow.paint("Output"));
        log::info!("{:18}: {}", "Taxa", utils::fmt_num(&header.ntax));
        log::info!("{:18}: {}", "Loci", utils::fmt_num(&count));
        log::info!(
            "{:18}: {}",
            "Alignment length",
            utils::fmt_num(&header.nchar)
        );
        log::info!("{:18}: {}", "Alignment file", self.output.display());
        log::info!("{:18}: {}", "Partition file", &part_file.display(),);
    }
}

struct Concat<'a> {
    input_fmt: &'a InputFmt,
    alignment: SeqMatrix,
    datatype: &'a DataType,
    header: Header,
    partition: Vec<Partition>,
    files: &'a mut [PathBuf],
}

impl<'a> Concat<'a> {
    fn new(files: &'a mut [PathBuf], input_fmt: &'a InputFmt, datatype: &'a DataType) -> Self {
        Self {
            input_fmt,
            datatype,
            alignment: IndexMap::new(),
            header: Header::new(),
            partition: Vec::new(),
            files,
        }
    }

    fn concat_alignment(&mut self, spin: &ProgressBar) {
        alphanumeric_sort::sort_path_slice(self.files);
        spin.set_message("Indexing alignments...");
        let id = IDs::new(self.files, self.input_fmt, self.datatype).get_id_unique();
        spin.set_message("Concatenating alignments...");
        self.concat(&id);
        self.header.ntax = self.alignment.len();
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
        let aln = Sequence::new(file, self.datatype);
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
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::helper::finder::Files;

    const DNA: DataType = DataType::Dna;

    #[test]
    fn test_concat_nexus() {
        let path = Path::new("test_files/concat/");
        let mut files = Files::new(path, &InputFmt::Nexus).get_files();
        let mut concat = Concat::new(&mut files, &InputFmt::Nexus, &DNA);
        let spin = utils::set_spinner();
        concat.concat_alignment(&spin);
        assert_eq!(3, concat.alignment.len());
    }

    #[test]
    #[should_panic]
    fn test_get_alignment_panic() {
        let path = Path::new("test_files/concat/");
        let mut files = Files::new(path, &InputFmt::Nexus).get_files();
        let concat = Concat::new(&mut files, &InputFmt::Nexus, &DNA);
        concat.get_alignment(Path::new("."));
    }

    #[test]
    fn test_concat_check_result() {
        let path = Path::new("test_files/concat/");
        let mut files = Files::new(path, &InputFmt::Nexus).get_files();
        let mut concat = Concat::new(&mut files, &InputFmt::Nexus, &DNA);
        let spin = utils::set_spinner();
        concat.concat_alignment(&spin);
        let abce = concat.alignment.get("ABCE").unwrap();
        let res = "??????????????gatattagtata";
        assert_eq!(res, abce);
    }

    #[test]
    fn test_concat_partition() {
        let path = Path::new("test_files/concat/");
        let mut files = Files::new(path, &InputFmt::Nexus).get_files();
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
}