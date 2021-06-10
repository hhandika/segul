use std::collections::BTreeMap;

pub fn read_phylip(path: &str) {
    let phylip = Phylip::new();
    phylip.read(path);
}

#[allow(dead_code)]
struct Phylip {
    matrix: BTreeMap<String, String>,
    num_tax: usize,
    num_seq: usize,
}

impl Phylip {
    fn new() -> Self {
        Self {
            matrix: BTreeMap::new(),
            num_tax: 0,
            num_seq: 0,
        }
    }

    fn read(&self, path: &str) {
        println!("Phylip {}", path);
    }
}
