//! Multi-Sequence Alignment Module.
//! Contains methods for working with multi-sequence alignments.

use std::io::{self, Result, Write};
use std::iter;
use std::path::{Path, PathBuf};

use indexmap::{IndexMap, IndexSet};
use indicatif::ProgressBar;

use crate::alignment::Alignment;
use crate::common::{Header, Partition, PartitionFormat, SeqFormat};
use crate::finder::{Files, IDs};
use crate::utils;
use crate::writer::SeqWriter;

pub struct MSAlignment<'a> {
    dir: &'a str,
    output: &'a str,
    input_format: &'a SeqFormat,
    output_format: SeqFormat,
    part_format: PartitionFormat,
}

impl<'a> MSAlignment<'a> {
    pub fn new(
        dir: &'a str,
        output: &'a str,
        input_format: &'a SeqFormat,
        output_format: SeqFormat,
        part_format: PartitionFormat,
    ) -> Self {
        Self {
            dir,
            output,
            input_format,
            output_format,
            part_format,
        }
    }

    pub fn concat_alignment(&self) {
        let mut concat = self.get_aln_format();
        let spin = concat.concat_alignment(self.dir);
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

    fn get_aln_format(&self) -> Concat {
        match self.input_format {
            SeqFormat::Fasta | SeqFormat::FastaInt => Concat::new(SeqFormat::Fasta, false),
            SeqFormat::Nexus | SeqFormat::NexusInt => Concat::new(SeqFormat::Nexus, false),
            SeqFormat::Phylip => Concat::new(SeqFormat::Phylip, false),
            SeqFormat::PhylipInt => Concat::new(SeqFormat::Phylip, true),
        }
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

struct Concat {
    input_format: SeqFormat,
    interleave: bool,
    alignment: IndexMap<String, String>,
    header: Header,
    partition: Vec<Partition>,
    files: Vec<PathBuf>,
}

impl Concat {
    fn new(input_format: SeqFormat, interleave: bool) -> Self {
        Self {
            input_format,
            interleave,
            alignment: IndexMap::new(),
            header: Header::new(),
            partition: Vec::new(),
            files: Vec::new(),
        }
    }

    fn concat_alignment(&mut self, dir: &str) -> ProgressBar {
        let spin = utils::set_spinner();
        self.files = Files::new(dir, &self.input_format).get_files();
        self.files.sort();
        spin.set_message("Indexing alignments...");
        let id = IDs::new(&self.files, &self.input_format).get_id_all(self.interleave);
        spin.set_message("Concatenating alignments...");
        let (alignment, nchar, partition) = self.concat(&id);
        self.alignment = alignment;
        self.header.ntax = self.alignment.len();
        self.header.nchar = nchar;
        self.partition = partition;
        spin
    }

    fn get_alignment(&self, file: &Path) -> Alignment {
        let mut aln = Alignment::new();
        aln.get_aln_any(file, &self.input_format, self.interleave);
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

    #[test]
    fn concat_nexus_test() {
        let path = "test_files/concat/";
        let mut concat = Concat::new(SeqFormat::Nexus, false);
        concat.concat_alignment(path);
        assert_eq!(3, concat.alignment.len());
    }

    #[test]
    fn concat_check_result_test() {
        let path = "test_files/concat/";
        let mut concat = Concat::new(SeqFormat::Nexus, false);
        concat.concat_alignment(path);
        let abce = concat.alignment.get("ABCE").unwrap();
        let res = "??????????????gatattagtata";
        assert_eq!(res, abce);
    }

    #[test]
    fn concat_partition_test() {
        let path = "test_files/concat/";
        let mut concat = Concat::new(SeqFormat::Nexus, false);
        concat.concat_alignment(path);
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
        assert_eq!(gaps, Concat::new(SeqFormat::Fasta, false).get_missings(len))
    }
}
