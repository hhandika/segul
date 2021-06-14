use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, LineWriter, Lines};
use std::path::Path;

use indexmap::IndexMap;

pub fn parse_fasta_id(path: &str) {
    let ids = get_ids(path);
    let mut records: Vec<IDs> = Vec::with_capacity(ids.len());
    ids.iter().for_each(|r| {
        let mut fasta = IDs::new();
        fasta.parse_id(r);
        records.push(fasta);
    });
    write_records(path, &mut records);
}

pub fn parse_fasta(path: &str) {
    let file = File::open(path).expect("CANNOT OPEN THE FILE");
    let buff = BufReader::new(file);

    let fas = get_records(buff);
    fas.iter().for_each(|(id, seq)| {
        println!("ID: {}", id);
        println!("Seq: {}", seq);
    })
}

#[allow(dead_code)]
pub fn get_stats(path: &str) {
    let file = File::open(path).expect("CANNOT OPEN THE FILE");
    let buff = BufReader::new(file);
    let fasta = FastaReader::new(buff);
    fasta.into_iter().for_each(|fast| {
        fast.get_seq_stats();
    });
}

fn get_records<R: Read>(buff: R) -> IndexMap<String, String> {
    let fasta = FastaReader::new(buff);
    let mut records: IndexMap<String, String> = IndexMap::new();
    fasta.into_iter().for_each(|fas| {
        #[allow(clippy::all)]
        if records.contains_key(&fas.id) {
            panic!("DUPLICATE SAMPLES. FIRST DUPLICATE FOUND: {}", fas.id);
        } else {
            records.insert(fas.id, fas.seq);
        }
    });

    records
}

fn get_ids<P: AsRef<Path>>(path: P) -> Vec<String> {
    let file = File::open(path).expect("CANNOT OPEN THE FILE");
    let buff = BufReader::new(file);
    buff.lines()
        .filter_map(|ok| ok.ok())
        .filter(|line| line.starts_with('>'))
        .map(|line| line.replace('>', ""))
        .collect::<Vec<String>>()
}

fn write_records(path: &str, records: &mut [IDs]) {
    let fstem = Path::new(path)
        .file_stem()
        .expect("CANNOT GET THE FILENAME FROM PATH");
    let fname = format!("{}.csv", fstem.to_string_lossy());
    let output = File::create(&fname).expect("CANNOT WRITE RESULTS");
    let mut line = LineWriter::new(output);
    records.sort_by(|a, b| a.genus.cmp(&b.genus));
    writeln!(line, "genus,species,voucher_id,geograpy").unwrap();
    records.iter().for_each(|r| {
        writeln!(
            line,
            "{},{},{},{}",
            r.genus, r.species, r.voucher_id, r.geography
        )
        .unwrap();
    });

    println!("DONE!");
    println!("The result is saved as {}", &fname);
}

struct IDs {
    genus: String,
    species: String,
    voucher_id: String,
    geography: String,
}

impl IDs {
    fn new() -> Self {
        Self {
            genus: String::new(),
            species: String::new(),
            voucher_id: String::new(),
            geography: String::new(),
        }
    }

    fn parse_id(&mut self, id: &str) {
        let ids: Vec<&str> = id.split('_').collect();
        assert!(
            ids.len() > 2,
            "INVALID ID FORMATS: {}; LEN: {}",
            id,
            ids.len()
        );
        self.genus.push_str(&ids[0]);
        self.species.push_str(&ids[1]);
        self.voucher_id.push_str(&ids[2]);
        if ids.len() > 3 {
            self.geography.push_str(&ids[3]);
        }
    }
}

pub struct Fasta {
    id: String,
    seq: String,
}

impl Fasta {
    fn new() -> Self {
        Self {
            id: String::new(),
            seq: String::new(),
        }
    }

    fn get_seq_stats(&self) {
        let (gc, len) = self.get_gc_content();
        println!("{}", self.id);
        println!("Sequence len\t: {} bp", len);
        println!("GC Content\t: {:.2}\n", gc as f64 / len as f64);
    }

    fn get_gc_content(&self) -> (usize, usize) {
        let mut gc = 0;
        let mut len = 0;

        self.seq.chars().for_each(|base| match base {
            'G' | 'g' | 'C' | 'c' => {
                gc += 1;
                len += 1;
            }
            '-' | '?' => (),
            _ => len += 1,
        });

        (gc, len)
    }
}

pub struct FastaReader<R> {
    reader: Lines<BufReader<R>>,
    pub id: Option<String>,
    pub seq: String,
}

impl<R: Read> FastaReader<R> {
    fn new(file: R) -> Self {
        Self {
            reader: BufReader::new(file).lines(),
            id: None,
            seq: String::new(),
        }
    }
}

impl<R: Read> Iterator for FastaReader<R> {
    type Item = Fasta;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(Ok(line)) = self.reader.next() {
            if let Some(id) = line.strip_prefix('>') {
                if self.id.is_some() {
                    let mut res = Fasta::new();
                    res.id.push_str(&self.id.as_ref().unwrap());
                    res.seq.push_str(&self.seq);
                    self.id = Some(String::from(id));
                    self.seq.clear();
                    return Some(res);
                } else {
                    self.id = Some(String::from(id));
                    self.seq.clear();
                }
                continue;
            }
            self.seq.push_str(line.trim());
        }
        if self.id.is_some() {
            let mut res = Fasta::new();
            res.id.push_str(&self.id.as_ref().unwrap());
            res.seq.push_str(&self.seq);
            self.id = None;
            self.seq.clear();
            self.seq.shrink_to_fit();
            return Some(res);
        }
        None
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn read_fasta_simple_test() {
        let path = "test_files/simple.fas";
        let file = File::open(path).unwrap();
        let buff = BufReader::new(file);
        let res = get_records(buff);

        assert_eq!(2, res.len());
    }
}
