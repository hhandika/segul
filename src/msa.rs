//!Multi-Sequence Alignment Module.
//! Contains methods for working with multi-sequence alignments.

use std::io::{self, Result, Write};
use std::iter;
use std::path::{Path, PathBuf};

use indexmap::IndexMap;
use indexmap::IndexSet;

use crate::alignment::Alignment;
use crate::common::{Header, InputFormat, OutputFormat, Partition, PartitionFormat};
use crate::finder::{Files, IDs};
use crate::utils;
use crate::writer::SeqWriter;

pub struct MSAlignment<'a> {
    dir: &'a str,
    output: &'a str,
    output_format: OutputFormat,
    part_format: PartitionFormat,
}

impl<'a> MSAlignment<'a> {
    pub fn new(
        dir: &'a str,
        output: &'a str,
        output_format: OutputFormat,
        part_format: PartitionFormat,
    ) -> Self {
        Self {
            dir,
            output,
            output_format,
            part_format,
        }
    }

    pub fn concat_nexus(&mut self) {
        let mut nex = Concat::new(InputFormat::Nexus);
        nex.concat_alignment(self.dir);
        let header = nex.get_header();
        self.write_alignment(&nex.alignment, &nex.partition, header);
    }

    pub fn concat_phylip(&mut self) {
        let mut phy = Concat::new(InputFormat::Phylip);
        phy.concat_alignment(self.dir);
        let header = phy.get_header();
        self.write_alignment(&phy.alignment, &phy.partition, header);
    }

    pub fn concat_fasta(&mut self) {
        let mut fas = Concat::new(InputFormat::Fasta);
        fas.concat_alignment(self.dir);
        let header = fas.get_header();
        self.write_alignment(&fas.alignment, &fas.partition, header);
    }

    fn write_alignment(&self, aln: &IndexMap<String, String>, part: &[Partition], header: Header) {
        let output = Path::new(self.output);
        self.display_alignment_stats(part.len(), &header).unwrap();
        let mut save = SeqWriter::new(output, aln, header, Some(part), &self.part_format);

        match self.output_format {
            OutputFormat::Nexus => save.write_sequence(&self.output_format),
            OutputFormat::Phylip => save.write_sequence(&self.output_format),
            OutputFormat::Fasta => save.write_fasta(),
        };

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

struct Concat {
    input: InputFormat,
    alignment: IndexMap<String, String>,
    ntax: usize,
    nchar: usize,
    datatype: String,
    missing: char,
    gap: char,
    partition: Vec<Partition>,
    files: Vec<PathBuf>,
}

impl Concat {
    fn new(input: InputFormat) -> Self {
        Self {
            input,
            alignment: IndexMap::new(),
            datatype: String::from("dna"),
            ntax: 0,
            nchar: 0,
            missing: '?',
            gap: '-',
            partition: Vec::new(),
            files: Vec::new(),
        }
    }

    fn concat_alignment(&mut self, dir: &str) {
        self.files = Files::new(dir, &self.input).get_files();
        self.files.sort();
        let id = IDs::new(&self.files, &self.input).get_id_all();
        let (alignment, nchar, partition) = self.concat(&id);
        self.alignment = alignment;
        self.ntax = self.alignment.len();
        self.nchar = nchar;
        self.partition = partition;
    }

    fn get_header(&self) -> Header {
        let mut header = Header::new();
        header.ntax = self.ntax;
        header.nchar = self.nchar;
        header.datatype = Some(self.datatype.clone());
        header.missing = Some(self.missing);
        header.gap = Some(self.gap);
        header
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
            let mut aln = Alignment::new();
            aln.get_aln_any(file, &self.input);
            nchar += aln.nchar; // increment sequence length using the value from parser
            let gene_name = file.file_stem().unwrap().to_string_lossy();
            self.get_partition(&mut partition, &gene_name, gene_start, nchar);
            gene_start = nchar + 1;
            id.iter().for_each(|id| {
                if !aln.alignment.contains_key(id) {
                    let seq = self.get_gaps(aln.nchar);
                    self.insert_alignment(&mut alignment, id, seq)
                } else {
                    let seq = aln.alignment.get(id).unwrap().to_string();
                    self.insert_alignment(&mut alignment, id, seq)
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

    fn insert_alignment(&self, alignment: &mut IndexMap<String, String>, id: &str, values: String) {
        if !alignment.contains_key(id) {
            alignment.insert(id.to_string(), values);
        } else if let Some(value) = alignment.get_mut(id) {
            value.push_str(&values);
        }
    }

    fn get_gaps(&self, len: usize) -> String {
        iter::repeat('-').take(len).collect()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn concat_nexus_test() {
        let path = "test_files/concat/";
        let mut concat = Concat::new(InputFormat::Nexus);
        concat.concat_alignment(path);
        assert_eq!(3, concat.alignment.len());
    }

    #[test]
    fn concat_check_result_test() {
        let path = "test_files/concat/";
        let mut concat = Concat::new(InputFormat::Nexus);
        concat.concat_alignment(path);
        let abce = concat.alignment.get("ABCE").unwrap();
        let res = "--------------gatattagtata";
        assert_eq!(res, abce);
    }

    #[test]
    fn concat_partition_test() {
        let path = "test_files/concat/";
        let mut concat = Concat::new(InputFormat::Nexus);
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
        let gaps = "-----";
        assert_eq!(gaps, Concat::new(InputFormat::Fasta).get_gaps(len))
    }
}
