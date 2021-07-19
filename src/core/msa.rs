//! Multi-Sequence Alignment Module.
//! Contains methods for working with multi-sequence alignments.

use std::io::{self, Result, Write};
use std::iter;
use std::path::{Path, PathBuf};

use indexmap::{IndexMap, IndexSet};
use indicatif::ProgressBar;

use crate::helper::common::{DataType, Header, InputFmt, OutputFmt, Partition, PartitionFmt};
use crate::helper::finder::IDs;
use crate::helper::sequence::Sequence;
use crate::helper::utils;
use crate::writer::seqwriter::SeqWriter;

pub struct MSAlignment<'a> {
    input_fmt: &'a InputFmt,
    output: &'a str,
    output_fmt: &'a OutputFmt,
    part_fmt: &'a PartitionFmt,
}

impl<'a> MSAlignment<'a> {
    pub fn new(
        input_fmt: &'a InputFmt,
        output: &'a str,
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

    pub fn concat_alignment(&self, files: &mut [PathBuf], datatype: &DataType) {
        let mut concat = Concat::new(files, &self.input_fmt, datatype);
        let spin = utils::set_spinner();
        self.write_alignment(&mut concat, &spin);
    }

    fn write_alignment(&self, concat: &mut Concat, spin: &ProgressBar) {
        concat.concat_alignment(&spin);
        let output = Path::new(self.output);
        let mut save = SeqWriter::new(
            output,
            &concat.alignment,
            concat.header.clone(),
            Some(&concat.partition),
            &self.part_fmt,
        );
        spin.set_message("Writing output files...");
        save.write_sequence(&self.output_fmt)
            .expect("Failed writing the output file");
        spin.finish_with_message("DONE!\n");
        self.print_alignment_stats(concat.partition.len(), &concat.header)
            .unwrap();
        save.print_save_path()
            .expect("Cannot write save path to stdout");
        save.print_partition_path()
            .expect("Cannot write partition path to stdout");
    }

    fn print_alignment_stats(&self, count: usize, header: &Header) -> Result<()> {
        let io = io::stdout();
        let mut writer = io::BufWriter::new(io);
        writeln!(writer, "\x1b[0;33mAlignment\x1b[0m")?;
        writeln!(writer, "Taxa\t\t: {}", utils::fmt_num(&header.ntax))?;
        writeln!(writer, "Loci\t\t: {}", utils::fmt_num(&count))?;
        writeln!(writer, "Length\t\t: {}", utils::fmt_num(&header.nchar))?;

        Ok(())
    }
}

struct Concat<'a> {
    input_fmt: &'a InputFmt,
    alignment: IndexMap<String, String>,
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
        let id = IDs::new(&self.files, &self.input_fmt, self.datatype).get_id_all();
        spin.set_message("Concatenating alignments...");
        let (alignment, nchar, partition) = self.concat(&id);
        self.alignment = alignment;
        self.header.ntax = self.alignment.len();
        self.header.nchar = nchar;
        self.partition = partition;
    }

    fn get_alignment(&self, file: &Path) -> Sequence {
        let mut aln = Sequence::new();
        aln.get_alignment(file, &self.input_fmt, self.datatype);
        assert!(
            aln.header.ntax != 0,
            "Found an empty alignment {}",
            file.display()
        );
        aln
    }

    fn concat(
        &mut self,
        id: &IndexSet<String>,
    ) -> (IndexMap<String, String>, usize, Vec<Partition>) {
        let mut alignment = IndexMap::new();
        let mut nchar = 0;
        let mut gene_start = 1;
        let mut partition = Vec::new();
        self.files.iter().for_each(|file| {
            let aln = self.get_alignment(file);
            nchar += aln.header.nchar; // increment sequence length using the value from parser
            self.get_partition(&mut partition, &aln.name, gene_start, nchar);
            gene_start = nchar + 1;
            id.iter().for_each(|id| match aln.alignment.get(id) {
                Some(seq) => {
                    self.insert_alignment(&mut alignment, id, seq);
                }
                None => {
                    let seq = self.get_missings(aln.header.nchar);
                    self.insert_alignment(&mut alignment, id, &seq);
                }
            });
        });
        (alignment, nchar, partition)
    }

    fn get_partition(
        &self,
        partition: &mut Vec<Partition>,
        gene_name: &str,
        start: usize,
        end: usize,
    ) {
        let mut part = Partition::new();
        part.gene = gene_name.to_string();
        part.start = start;
        part.end = end;
        partition.push(part);
    }

    fn insert_alignment(&self, alignment: &mut IndexMap<String, String>, id: &str, seq: &str) {
        match alignment.get_mut(id) {
            Some(seqs) => seqs.push_str(seq),
            None => {
                alignment.insert(id.to_string(), seq.to_string());
            }
        }
    }

    fn get_missings(&self, len: usize) -> String {
        iter::repeat('?').take(len).collect()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::helper::finder::Files;

    const DNA: DataType = DataType::Dna;

    #[test]
    fn concat_nexus_test() {
        let path = "test_files/concat/";
        let mut files = Files::new(path, &InputFmt::Nexus).get_files();
        let mut concat = Concat::new(&mut files, &InputFmt::Nexus, &DNA);
        let spin = utils::set_spinner();
        concat.concat_alignment(&spin);
        assert_eq!(3, concat.alignment.len());
    }

    #[test]
    #[should_panic]
    fn get_alignment_panic_test() {
        let path = "test_files/concat/";
        let mut files = Files::new(path, &InputFmt::Nexus).get_files();
        let concat = Concat::new(&mut files, &InputFmt::Nexus, &DNA);
        concat.get_alignment(Path::new("."));
    }

    #[test]
    fn concat_check_result_test() {
        let path = "test_files/concat/";
        let mut files = Files::new(path, &InputFmt::Nexus).get_files();
        let mut concat = Concat::new(&mut files, &InputFmt::Nexus, &DNA);
        let spin = utils::set_spinner();
        concat.concat_alignment(&spin);
        let abce = concat.alignment.get("ABCE").unwrap();
        let res = "??????????????gatattagtata";
        assert_eq!(res, abce);
    }

    #[test]
    fn concat_partition_test() {
        let path = "test_files/concat/";
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
    fn get_gaps_test() {
        let len = 5;
        let gaps = "?????";
        assert_eq!(
            gaps,
            Concat::new(&mut [PathBuf::from(".")], &InputFmt::Fasta, &DNA).get_missings(len)
        )
    }
}
