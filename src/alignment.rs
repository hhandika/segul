use std::iter;
use std::path::{Path, PathBuf};

use glob::glob;
use indexmap::IndexMap;
use indexmap::IndexSet;

use crate::common::{Header, Partition, SeqFormat, SeqPartition};
use crate::nexus::Nexus;
use crate::phylip::Phylip;
use crate::writer::SeqWriter;

pub fn concat_nexus(dir: &str, outname: &str, filetype: SeqFormat, partition: SeqPartition) {
    let mut nex = Concat::new();
    let path = Path::new(dir).join(outname);
    nex.concat_from_nexus(dir);
    let header = nex.get_header();
    let mut save = SeqWriter::new(
        &path,
        &nex.alignment,
        header,
        Some(nex.partition),
        partition,
    );

    match filetype {
        SeqFormat::Nexus => save.write_sequence(&filetype),
        SeqFormat::Phylip => save.write_sequence(&filetype),
        SeqFormat::Fasta => save.write_fasta(),
    };
}

struct Concat {
    alignment: IndexMap<String, String>,
    ntax: usize,
    nchar: usize,
    datatype: String,
    missing: char,
    gap: char,
    partition: Vec<Partition>,
    files: Vec<PathBuf>,
}

#[allow(dead_code)]
impl Concat {
    fn new() -> Self {
        Self {
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
        let pattern = format!("{}/*.nex*", dir);
        self.get_files(&pattern);
        self.files.sort();
        let id = self.get_id_from_nexus();
        self.alignment = self.concat_nexus(&id);
        self.ntax = self.alignment.len();
    }

    fn concat_from_phylip(&mut self, dir: &str) {
        let pattern = format!("{}/*.phy*", dir);
        self.get_files(&pattern);
        self.files.sort();
        let id = self.get_id_from_phylip();
        self.alignment = self.concat_phylip(&id);
        self.ntax = self.alignment.len();
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

    fn get_id_from_phylip(&mut self) -> IndexSet<String> {
        let mut id = IndexSet::new();
        self.files.iter().for_each(|file| {
            let mut phy = Phylip::new(file);
            phy.read().expect("CANNOT READ A NEXUS FILE");
            self.get_id(&phy.matrix, &mut id);
        });

        id
    }

    fn get_id_from_nexus(&mut self) -> IndexSet<String> {
        let mut id = IndexSet::new();
        self.files.iter().for_each(|file| {
            let mut nex = Nexus::new(file);
            nex.read().expect("CANNOT READ A NEXUS FILE");
            self.get_id(&nex.matrix, &mut id);
        });

        id
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
            let mut nex = Nexus::new(file);
            nex.read().expect("CANNOT READ A NEXUS FILE");
            self.check_is_alignment(&file, nex.is_alignment);
            nchar += nex.nchar; // increment sequence length using the value from parser
            let gene_name = file.file_stem().unwrap().to_string_lossy();
            self.get_partition(&mut partition, &gene_name, gene_start, nchar);
            gene_start = nchar + 1;
            id.iter().for_each(|id| {
                if !nex.matrix.contains_key(id) {
                    let seq = self.get_gaps(nex.nchar);
                    self.insert_alignment(&mut alignment, id, seq)
                } else {
                    let seq = nex.matrix.get(id).unwrap().to_string();
                    self.insert_alignment(&mut alignment, id, seq)
                }
            });
        });
        self.nchar = nchar;
        self.partition = partition;
        alignment
    }

    fn concat_phylip(&mut self, id: &IndexSet<String>) -> IndexMap<String, String> {
        let mut alignment = IndexMap::new();
        let mut nchar = 0;
        let mut gene_start = 1;
        let mut partition = Vec::new();
        self.files.iter().for_each(|file| {
            let mut phy = Phylip::new(file);
            phy.read().expect("CANNOT READ A PHYLIP FILE");
            self.check_is_alignment(&file, phy.is_alignment);
            nchar += phy.nchar; // increment sequence length using the value from parser
            let gene_name = file.file_stem().unwrap().to_string_lossy();
            self.get_partition(&mut partition, &gene_name, gene_start, nchar);
            gene_start = nchar + 1;
            id.iter().for_each(|id| {
                if !phy.matrix.contains_key(id) {
                    let seq = self.get_gaps(phy.nchar);
                    self.insert_alignment(&mut alignment, id, seq)
                } else {
                    let seq = phy.matrix.get(id).unwrap().to_string();
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

    fn check_is_alignment(&self, file: &Path, aligned: bool) {
        if !aligned {
            panic!(
                "INVALID INPUT FILES. {} IS NOT AN ALIGNMENT",
                file.display()
            );
        }
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
        let mut concat = Concat::new();
        concat.get_files(&pattern);
        assert_eq!(4, concat.files.len());
    }

    #[test]
    fn concat_nexus_test() {
        let path = "test_files/concat/";
        let mut concat = Concat::new();
        concat.concat_from_nexus(path);
        assert_eq!(3, concat.alignment.len());
    }

    #[test]
    fn concat_check_result_test() {
        let path = "test_files/concat/";
        let mut concat = Concat::new();
        concat.concat_from_nexus(path);
        let abce = concat.alignment.get("ABCE").unwrap();
        let res = "--------------gatattagtata";
        assert_eq!(res, abce);
    }

    #[test]
    fn concat_partition_test() {
        let path = "test_files/concat/";
        let mut concat = Concat::new();
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
