use std::iter;
use std::path::{Path, PathBuf};

use glob::glob;
use linked_hash_map::LinkedHashMap;
use linked_hash_set::LinkedHashSet;

use crate::common::SeqFormat;
use crate::nexus::Nexus;
use crate::writer::SeqWriter;

pub fn concat_nexus(dir: &str, outname: &str, filetype: SeqFormat) {
    let mut nex = ConcatNexus::new();
    let path = Path::new(dir).join(outname);
    nex.concat_from_nexus(dir);
    let mut save = SeqWriter::new(
        &path,
        &nex.alignment,
        Some(nex.ntax),
        Some(nex.nchar),
        Some(nex.datatype),
        Some(nex.missing),
        Some(nex.gap),
    );

    match filetype {
        SeqFormat::Nexus => save.write_sequence(&filetype),
        SeqFormat::Phylip => save.write_sequence(&filetype),
        SeqFormat::Fasta => save.write_fasta(),
    };
}

#[allow(dead_code)]
struct ConcatNexus {
    start: usize,
    end: usize,
    genes_pos: LinkedHashMap<usize, usize>,
    alignment: LinkedHashMap<String, String>,
    ntax: usize,
    nchar: usize,
    datatype: String,
    missing: char,
    gap: char,
    id: LinkedHashSet<String>,
    files: Vec<PathBuf>,
}

#[allow(dead_code)]
impl ConcatNexus {
    fn new() -> Self {
        Self {
            start: 0,
            end: 0,
            genes_pos: LinkedHashMap::new(),
            alignment: LinkedHashMap::new(),
            datatype: String::from("dna"),
            ntax: 0,
            nchar: 0,
            missing: '?',
            gap: '-',
            id: LinkedHashSet::new(),
            files: Vec::new(),
        }
    }

    fn concat_from_nexus(&mut self, dir: &str) {
        let pattern = format!("{}/*.nex*", dir);
        self.get_files(&pattern);
        self.files.sort();
        self.id = self.get_id_from_nexus();
        self.alignment = self.concat();
        self.ntax = self.alignment.len();
    }

    fn get_files(&mut self, pattern: &str) {
        self.files = glob(pattern)
            .expect("COULD NOT FIND FILES")
            .filter_map(|ok| ok.ok())
            .collect();
    }

    fn get_id_from_nexus(&mut self) -> LinkedHashSet<String> {
        let mut id = LinkedHashSet::new();
        self.files.iter().for_each(|file| {
            let mut nex = Nexus::new();
            nex.read(file).expect("CANNOT READ A NEXUS FILE");
            nex.matrix.keys().for_each(|key| {
                if !id.contains(key) {
                    id.insert(key.to_string());
                }
            });
        });

        id
    }

    fn concat(&mut self) -> LinkedHashMap<String, String> {
        let mut alignment = LinkedHashMap::new();
        let mut nchar = 0;
        self.files.iter().for_each(|file| {
            let mut nex = Nexus::new();
            nex.read(file).expect("CANNOT READ A NEXUS FILE");
            nchar += nex.nchar;
            self.id.iter().for_each(|id| {
                if !nex.matrix.contains_key(id) {
                    let seq = self.get_gaps(nex.nchar);
                    self.insert_alignment(&mut alignment, id, seq)
                } else {
                    let seq = nex.matrix.get(id).unwrap().to_string();
                    self.insert_alignment(&mut alignment, id, seq)
                }
            })
        });
        self.nchar = nchar;
        alignment
    }

    fn insert_alignment(
        &self,
        alignment: &mut LinkedHashMap<String, String>,
        id: &str,
        values: String,
    ) {
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
        let mut concat = ConcatNexus::new();
        concat.get_files(&pattern);
        assert_eq!(4, concat.files.len());
    }

    #[test]
    fn concat_nexus_test() {
        let path = "test_files/concat/";
        let mut concat = ConcatNexus::new();
        concat.concat_from_nexus(path);
        let alignment = concat.concat();
        assert_eq!(3, alignment.len());
    }

    #[test]
    fn concat_check_result_test() {
        let path = "test_files/concat/";
        let mut concat = ConcatNexus::new();
        concat.concat_from_nexus(path);
        let alignment = concat.concat();
        let abce = alignment.get("ABCE").unwrap();
        let res = "--------------gatattagtata";
        assert_eq!(res, abce);
    }
}
