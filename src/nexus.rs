// use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;

pub fn parse_nexus_id(path: &str) {
    let file = File::open(path).expect("CANNOT OPEN THE FILE");
    let buff = BufReader::new(file);

    let mut matrix = false;
    let mut inblock = false;

    // buff.read_line(buff);
    buff.lines().filter_map(|ok| ok.ok()).for_each(|l| {
        match l.to_lowercase() {
            l if l.starts_with("#nexus") => (),
            l if l.starts_with("begin") => inblock = true,
            l if l.starts_with("matrix") => matrix = true,
            l if l.contains(';') => matrix = false,
            _ => (),
        }

        if matrix && inblock {
            if !l.to_lowercase().starts_with("matrix") && l.len() > 1 {
                let cigar: Vec<&str> = l.split_whitespace().collect();
                println!("ID: {}", cigar[0].trim());
                println!("{}", cigar[1].trim());
            }
        }
    });

    // let fasta = NexusReader::new(buff);
    // fasta.into_iter().for_each(|fast| {
    //     println!("ID: {}", fast.id);
    //     println!("{}", fast.seq);
    // });
}

// struct Nexus {
//     id: String,
//     seq: String,
//     // ntax: usize,
//     // nchar: usize,
// }

// impl Nexus {
//     fn new() -> Self {
//         Self {
//             id: String::new(),
//             seq: String::new(),
//             // ntax: 0,
//             // nchar: 0,
//         }
//     }
// }

// struct NexusReader<R> {
//     inblock: bool,
//     matrix: bool,
//     id: String,
//     seq: String,
//     reader: Lines<BufReader<R>>,
// }

// impl<R: Read> NexusReader<R> {
//     fn new(file: R) -> Self {
//         Self {
//             inblock: false,
//             matrix: false,
//             id: String::new(),
//             seq: String::new(),
//             reader: BufReader::new(file).lines(),
//         }
//     }

//     // fn get_blocks(&self, line: &str) {}
// }

// impl<R: Read> Iterator for NexusReader<R> {
//     type Item = Nexus;

//     fn next(&mut self) -> Option<Self::Item> {
//         while let Some(Ok(line)) = self.reader.next() {
//             match line.to_lowercase() {
//                 l if l.starts_with("#nexus") => continue,
//                 l if l.starts_with("begin") => self.inblock = true,
//                 l if l.starts_with("matrix") => self.inblock = true,
//                 l if l.starts_with("end") => self.inblock = false,
//                 _ => (),
//             }

//             if self.inblock && self.matrix {
//                 let mut res = Nexus::new();
//                 let cigar: Vec<&str> = line.split('\t').collect();
//                 res.id.push_str(&self.id);
//                 res.seq.push_str(&self.seq);
//                 self.seq = Some(&cigar[1]);
//                 return Some(res);
//             } else {
//                  None
//             }
//         }
//     }
// }
