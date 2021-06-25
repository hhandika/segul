//! Multi-Sequence Alignment Module.
//! Contains methods for working with multi-sequence alignments.

use std::io::{self, Result, Write};
use std::iter;
use std::path::{Path, PathBuf};

use alphanumeric_sort;
use indexmap::{IndexMap, IndexSet};
use indicatif::ProgressBar;

use crate::alignment::Alignment;
use crate::common::{Header, Partition, PartitionFormat, SeqFormat};
use crate::finder::IDs;
use crate::utils;
use crate::writer::SeqWriter;

pub struct MSAlignment<'a> {
    input_format: &'a SeqFormat,
    output: &'a str,
    output_format: SeqFormat,
    part_format: PartitionFormat,
}

impl<'a> MSAlignment<'a> {
    pub fn new(
        input_format: &'a SeqFormat,
        output: &'a str,
        output_format: SeqFormat,
        part_format: PartitionFormat,
    ) -> Self {
        Self {
            input_format,
            output,
            output_format,
            part_format,
        }
    }

    pub fn concat_alignment(&self, files: &mut [PathBuf]) {
        let mut concat = Concat::new(files, &self.input_format);
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
            &self.part_format,
        );
        spin.set_message("Writing output files...");
        save.write_sequence(&self.output_format);
        spin.finish_with_message("DONE!\n");
        self.display_alignment_stats(concat.partition.len(), &concat.header)
            .unwrap();
        save.display_save_path();
        save.display_partition_path();
    }

    fn display_alignment_stats(&self, count: usize, header: &Header) -> Result<()> {
        let io = io::stdout();
        let mut writer = io::BufWriter::new(io);
        writeln!(writer, "\x1b[0;33mAlignment\x1b[0m")?;
        writeln!(writer, "#Loci\t\t: {}", utils::format_thousand_sep(&count))?;
        writeln!(
            writer,
            "#Taxa\t\t: {}",
            utils::format_thousand_sep(&header.ntax)
        )?;
        writeln!(
            writer,
            "#Chars\t\t: {} bp",
            utils::format_thousand_sep(&header.nchar)
        )?;

        Ok(())
    }
}

struct Concat<'a> {
    input_format: &'a SeqFormat,
    alignment: IndexMap<String, String>,
    header: Header,
    partition: Vec<Partition>,
    files: &'a mut [PathBuf],
}

impl<'a> Concat<'a> {
    fn new(files: &'a mut [PathBuf], input_format: &'a SeqFormat) -> Self {
        Self {
            input_format,
            alignment: IndexMap::new(),
            header: Header::new(),
            partition: Vec::new(),
            files,
        }
    }

    fn concat_alignment(&mut self, spin: &ProgressBar) {
        alphanumeric_sort::sort_path_slice(self.files);
        spin.set_message("Indexing alignments...");
        let id = IDs::new(&self.files, &self.input_format).get_id_all();
        spin.set_message("Concatenating alignments...");
        let (alignment, nchar, partition) = self.concat(&id);
        self.alignment = alignment;
        self.header.ntax = self.alignment.len();
        self.header.nchar = nchar;
        self.partition = partition;
    }

    fn get_alignment(&self, file: &Path) -> Alignment {
        let mut aln = Alignment::new();
        aln.get_aln_any(file, &self.input_format);
        assert!(
            aln.header.ntax != 0,
            "ZERO TAXON FOUND IN ALIGNMENT FILES {}",
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
    use crate::finder::Files;

    #[test]
    fn concat_nexus_test() {
        let path = "test_files/concat/";
        let mut files = Files::new(path, &SeqFormat::Nexus).get_files();
        let mut concat = Concat::new(&mut files, &SeqFormat::Nexus);
        let spin = utils::set_spinner();
        concat.concat_alignment(&spin);
        assert_eq!(3, concat.alignment.len());
    }

    #[test]
    #[should_panic]
    fn get_alignment_panic_test() {
        let path = "test_files/concat/";
        let mut files = Files::new(path, &SeqFormat::Nexus).get_files();
        let concat = Concat::new(&mut files, &SeqFormat::Nexus);
        concat.get_alignment(Path::new("."));
    }

    #[test]
    fn concat_check_result_test() {
        let path = "test_files/concat/";
        let mut files = Files::new(path, &SeqFormat::Nexus).get_files();
        let mut concat = Concat::new(&mut files, &SeqFormat::Nexus);
        let spin = utils::set_spinner();
        concat.concat_alignment(&spin);
        let abce = concat.alignment.get("ABCE").unwrap();
        let res = "??????????????gatattagtata";
        assert_eq!(res, abce);
    }

    #[test]
    fn concat_partition_test() {
        let path = "test_files/concat/";
        let mut files = Files::new(path, &SeqFormat::Nexus).get_files();
        let mut concat = Concat::new(&mut files, &SeqFormat::Nexus);
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
            Concat::new(&mut [PathBuf::from(".")], &SeqFormat::Fasta).get_missings(len)
        )
    }
}
