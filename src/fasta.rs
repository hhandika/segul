use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, LineWriter};
use std::path::Path;

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

fn get_ids(path: &str) -> Vec<String> {
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
