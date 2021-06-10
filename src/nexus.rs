use std::collections::BTreeMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader, LineWriter, Lines, Read, Result};
use std::path::Path;

pub fn read_nexus<P: AsRef<Path>>(path: &P) {
    let input = File::open(path).unwrap();
    let buff = BufReader::new(input);
    let mut nex = NexusReader::new();
    nex.read(buff).unwrap();

    let matrix = nex.parse_matrix();

    matrix.iter().for_each(|(id, seq)| {
        println!(">{}", id);
        println!("{}", seq);
    });
}

pub fn convert_to_fasta(path: &str) {
    let input = File::open(path).unwrap();
    let buff = BufReader::new(input);
    let mut nex = NexusReader::new();

    nex.read(buff).expect("CANNOT READ NEXUS FILES");
    let matrix = nex.parse_matrix();
    write_fasta(matrix, path)
}

fn write_fasta(matrix: BTreeMap<String, String>, path: &str) {
    let name = Path::new(path).file_stem().unwrap();
    let fname = format!("{}.fas", name.to_string_lossy());
    let file = File::create(fname).expect("CANNOT CREATE FASTA FILE");
    let mut writer = LineWriter::new(file);

    matrix.iter().for_each(|(id, seq)| {
        writeln!(writer, ">{}", id).unwrap();
        writeln!(writer, "{}", seq).unwrap();
    });
}

struct NexusReader {
    matrix: String,
}

impl NexusReader {
    fn new() -> Self {
        Self {
            matrix: String::new(),
        }
    }

    fn read<R: Read>(&mut self, reader: R) -> Result<()> {
        let reader = Reader::new(reader);
        reader.into_iter().for_each(|r| {
            if r.to_lowercase().contains("matrix") {
                self.matrix = r.trim().to_string();
            }
        });

        Ok(())
    }

    fn parse_matrix(&mut self) -> BTreeMap<String, String> {
        self.matrix.pop(); // remove terminated semicolon.
        let matrix: Vec<&str> = self.matrix.split('\n').collect();
        let mut seqs = BTreeMap::new();
        matrix[1..]
            .iter()
            .filter(|l| !l.is_empty())
            .map(|l| l.trim())
            .for_each(|l| {
                let seq: Vec<&str> = l.split_whitespace().collect();
                if seq.len() == 2 {
                    let id = seq[0].to_string();
                    let dna = seq[1].to_string();
                    if seqs.contains_key(&id) {
                        panic!("DUPLICATE SAMPLES. FIRST DUPLICATE FOUND: {}", id);
                    } else {
                        seqs.insert(id, dna);
                    }
                }
            });
        seqs
    }
}

struct Reader<R> {
    reader: Lines<BufReader<R>>,
    buffer: String,
    content: String,
}

impl<R: Read> Reader<R> {
    fn new(file: R) -> Self {
        Self {
            reader: BufReader::new(file).lines(),
            buffer: String::new(),
            content: String::new(),
        }
    }
}

// Iterate over the file.
// Collect each of the nexus block terminated by semi-colon.
impl<R: Read> Iterator for Reader<R> {
    type Item = String;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(Ok(line)) = self.reader.next() {
            self.buffer.push_str(&line);
            if !line.is_empty() {
                self.buffer.push('\n');
            }
            if line.ends_with(';') {
                self.content.push_str(&self.buffer);
                self.buffer.clear();
            }
            let token = self.content.trim().to_string();
            self.content.clear();
            Some(token)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    #[should_panic]
    fn delimit_matrix_test() {
        let sample = "test_files/duplicates.nex";
        read_nexus(&sample);
    }

    // #[test]
    // fn regex_seq_id_test() {
    //     let text = "agga--";
    //     let re = regex_seq_id();
    //     assert_eq!(true, re.is_match(text))
    // }

    // #[test]
    // fn regex_seq_end_test() {
    //     let text = "agga--\n";
    //     let re = regex_seq_id();
    //     assert_eq!(true, re.is_match(text))
    // }
}
