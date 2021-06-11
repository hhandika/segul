use std::collections::BTreeMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::LineWriter;
use std::path::Path;

pub fn write_fasta(matrix: &BTreeMap<String, String>, path: &str) {
    let name = Path::new(path).file_stem().unwrap();
    let fname = format!("{}.fas", name.to_string_lossy());
    let file = File::create(fname).expect("CANNOT CREATE FASTA FILE");
    let mut writer = LineWriter::new(file);

    matrix.iter().for_each(|(id, seq)| {
        writeln!(writer, ">{}", id).unwrap();
        writeln!(writer, "{}", seq).unwrap();
    });
}
