use std::iter;
use std::path::{Path, PathBuf};

use glob::glob;
use indexmap::IndexMap;
use indexmap::IndexSet;

use crate::common::{Header, OutputFormat, Partition, PartitionFormat};
use crate::fasta::Fasta;
use crate::nexus::Nexus;
use crate::phylip::Phylip;
use crate::writer::SeqWriter;

/// Multi-Sequence Alignment Module.
/// Contains methods for working with multi-sequence alignments.

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
        nex.concat_from_nexus(self.dir);
        let header = nex.get_header();
        self.write_alignment(&nex.alignment, &nex.partition, header);
    }

    pub fn concat_phylip(&mut self) {
        let mut phy = Concat::new(InputFormat::Phylip);
        phy.concat_from_nexus(self.dir);
        let header = phy.get_header();
        self.write_alignment(&phy.alignment, &phy.partition, header);
    }

    pub fn concat_fasta(&mut self) {
        let mut fas = Concat::new(InputFormat::Fasta);
        fas.concat_from_nexus(self.dir);
        let header = fas.get_header();
        self.write_alignment(&fas.alignment, &fas.partition, header);
    }

    fn write_alignment(&self, aln: &IndexMap<String, String>, part: &[Partition], header: Header) {
        let path = Path::new(self.dir).join(self.output);
        let mut save = SeqWriter::new(&path, aln, header, Some(part), &self.part_format);

        match self.output_format {
            OutputFormat::Nexus => save.write_sequence(&self.output_format),
            OutputFormat::Phylip => save.write_sequence(&self.output_format),
            OutputFormat::Fasta => save.write_fasta(),
        };
    }
}

struct Alignment {
    alignment: IndexMap<String, String>,
    nchar: usize,
}

impl Alignment {
    fn new() -> Self {
        Self {
            alignment: IndexMap::new(),
            nchar: 0,
        }
    }

    fn get_aln_from_nexus(&mut self, file: &Path) {
        let mut nex = Nexus::new(file);
        nex.read().expect("CANNOT READ A NEXUS FILE");
        self.check_is_alignment(&file, nex.is_alignment);
        self.get_alignment(nex.matrix, nex.nchar)
    }

    fn get_aln_from_phylip(&mut self, file: &Path) {
        let mut phy = Phylip::new(file);
        phy.read().expect("CANNOT READ A PHYLIP FILE");
        self.check_is_alignment(file, phy.is_alignment);
        self.get_alignment(phy.matrix, phy.nchar);
    }

    fn get_aln_from_fasta(&mut self, file: &Path) {
        let mut fas = Fasta::new(file);
        fas.read();
        self.check_is_alignment(file, fas.is_alignment);
        let nchar = self.get_nchar(&fas.matrix);
        self.get_alignment(fas.matrix, nchar);
    }

    fn get_alignment(&mut self, alignment: IndexMap<String, String>, nchar: usize) {
        self.alignment = alignment;
        self.nchar = nchar;
    }

    // Count char for Fasta. Get the char length from the first sequence.
    // This is fine since Fasta struct check for the char
    // is the same length.
    fn get_nchar(&mut self, alignment: &IndexMap<String, String>) -> usize {
        alignment.values().next().unwrap().len()
    }

    fn check_is_alignment(&self, file: &Path, aligned: bool) {
        if !aligned {
            panic!(
                "INVALID INPUT FILES. {} IS NOT AN ALIGNMENT",
                file.display()
            );
        }
    }
}

enum InputFormat {
    Nexus,
    Phylip,
    Fasta,
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

    fn concat_from_nexus(&mut self, dir: &str) {
        let pattern = self.get_pattern(dir);
        self.get_files(&pattern);
        self.check_glob_results();
        self.files.sort();
        let id = self.get_id_all();
        self.alignment = self.concat_nexus(&id);
        self.ntax = self.alignment.len();
    }

    fn get_pattern(&self, dir: &str) -> String {
        match self.input {
            InputFormat::Nexus => format!("{}/*.nex*", dir),
            InputFormat::Phylip => format!("{}/*.phy*", dir),
            InputFormat::Fasta => format!("{}/*.fa*", dir),
        }
    }

    fn get_header(&self) -> Header {
        let mut header = Header::new();
        header.ntax = Some(self.ntax);
        header.nchar = Some(self.nchar);
        header.datatype = Some(self.datatype.clone());
        header.missing = Some(self.missing);
        header.gap = Some(self.gap);
        header
    }

    fn get_files(&mut self, pattern: &str) {
        self.files = glob(pattern)
            .expect("COULD NOT FIND FILES")
            .filter_map(|ok| ok.ok())
            .collect();
    }

    fn check_glob_results(&self) {
        if self.files.is_empty() {
            panic!("NO VALID ALIGNMENT FILES FOUND");
        }
    }

    fn get_id_all(&self) -> IndexSet<String> {
        let mut id = IndexSet::new();
        match self.input {
            InputFormat::Nexus => self.get_id_from_nexus(&mut id),
            InputFormat::Phylip => self.get_id_from_phylip(&mut id),
            InputFormat::Fasta => self.get_id_from_fasta(&mut id),
        };
        id
    }

    fn get_id_from_phylip(&self, id: &mut IndexSet<String>) {
        self.files.iter().for_each(|file| {
            let mut phy = Phylip::new(file);
            phy.read().expect("CANNOT READ A NEXUS FILE");
            self.get_id(&phy.matrix, id);
        });
    }

    fn get_id_from_nexus(&self, id: &mut IndexSet<String>) {
        self.files.iter().for_each(|file| {
            let mut nex = Nexus::new(file);
            nex.read().expect("CANNOT READ A NEXUS FILE");
            self.get_id(&nex.matrix, id);
        });
    }

    fn get_id_from_fasta(&self, id: &mut IndexSet<String>) {
        self.files.iter().for_each(|file| {
            let mut fas = Fasta::new(file);
            fas.read();
            self.get_id(&fas.matrix, id);
        });
    }

    fn get_id(&self, matrix: &IndexMap<String, String>, id: &mut IndexSet<String>) {
        matrix.keys().for_each(|key| {
            if !id.contains(key) {
                id.insert(key.to_string());
            }
        });
    }

    fn concat_nexus(&mut self, id: &IndexSet<String>) -> IndexMap<String, String> {
        let mut alignment = IndexMap::new();
        let mut nchar = 0;
        let mut gene_start = 1;
        let mut partition = Vec::new();
        self.files.iter().for_each(|file| {
            let mut aln = Alignment::new();
            match self.input {
                InputFormat::Nexus => aln.get_aln_from_nexus(file),
                InputFormat::Phylip => aln.get_aln_from_phylip(file),
                InputFormat::Fasta => aln.get_aln_from_fasta(file),
            }
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
        self.nchar = nchar;
        self.partition = partition;
        alignment
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
    fn get_files_test() {
        let path = "test_files/concat/";
        let pattern = format!("{}/*.nex*", path);
        let mut concat = Concat::new(InputFormat::Nexus);
        concat.get_files(&pattern);
        assert_eq!(4, concat.files.len());
    }

    #[test]
    #[should_panic]
    fn check_empty_files_test() {
        let path = "test_files/empty/";
        let pattern = format!("{}/*.nex*", path);
        let mut concat = Concat::new(InputFormat::Nexus);
        concat.get_files(&pattern);
        concat.check_glob_results();
    }

    #[test]
    fn concat_nexus_test() {
        let path = "test_files/concat/";
        let mut concat = Concat::new(InputFormat::Nexus);
        concat.concat_from_nexus(path);
        assert_eq!(3, concat.alignment.len());
    }

    #[test]
    fn concat_check_result_test() {
        let path = "test_files/concat/";
        let mut concat = Concat::new(InputFormat::Nexus);
        concat.concat_from_nexus(path);
        let abce = concat.alignment.get("ABCE").unwrap();
        let res = "--------------gatattagtata";
        assert_eq!(res, abce);
    }

    #[test]
    fn concat_partition_test() {
        let path = "test_files/concat/";
        let mut concat = Concat::new(InputFormat::Nexus);
        concat.concat_from_nexus(path);
        assert_eq!(1, concat.partition[0].start);
        assert_eq!(6, concat.partition[0].end);
        assert_eq!(7, concat.partition[1].start);
        assert_eq!(14, concat.partition[1].end);
        assert_eq!(15, concat.partition[2].start);
        assert_eq!(20, concat.partition[2].end);
        assert_eq!(21, concat.partition[3].start);
        assert_eq!(26, concat.partition[3].end);
    }
}
