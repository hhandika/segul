//! Write sequences to file
use std::collections::BTreeMap;
use std::io::prelude::*;
use std::path::Path;

use anyhow::Result;
use indexmap::IndexMap;

use crate::helper::types::{Header, OutputFmt};
use crate::writer::FileWriter;

impl FileWriter for SeqWriter<'_> {}

pub struct SeqWriter<'a> {
    output: &'a Path,
    matrix: &'a IndexMap<String, String>,
    id_len: usize,
    header: &'a Header,
}

impl<'a> SeqWriter<'a> {
    pub fn new(output: &'a Path, matrix: &'a IndexMap<String, String>, header: &'a Header) -> Self {
        Self {
            output,
            id_len: 0,
            matrix,
            header,
        }
    }

    pub fn write_sequence(&mut self, output_fmt: &OutputFmt) -> Result<()> {
        match output_fmt {
            OutputFmt::Nexus => self.write_nexus(false)?,
            OutputFmt::NexusInt => self.write_nexus(true)?,
            OutputFmt::Phylip => self.write_phylip(false)?,
            OutputFmt::PhylipInt => self.write_phylip(true)?,
            OutputFmt::Fasta => self.write_fasta(false)?,
            OutputFmt::FastaInt => self.write_fasta(true)?,
        }

        Ok(())
    }

    fn write_fasta(&mut self, interleave: bool) -> Result<()> {
        let mut writer = self
            .create_output_file(self.output)
            .expect("Failed writing a fasta formatted file");
        let n = self.get_interleave_len();
        self.matrix.iter().for_each(|(id, seq)| {
            writeln!(writer, ">{}", id).unwrap();
            if !interleave {
                writeln!(writer, "{}", seq).unwrap();
            } else {
                let chunks = self.chunk_seq(seq, n);
                chunks.iter().for_each(|chunk| {
                    writeln!(writer, "{}", chunk).unwrap();
                })
            }
        });

        writer.flush()?;
        Ok(())
    }

    fn write_nexus(&mut self, interleave: bool) -> Result<()> {
        let mut writer = self
            .create_output_file(self.output)
            .expect("Failed writing a NEXUS formatted file");
        self.write_nex_header(&mut writer, interleave)?;

        // We write only instead of write line.
        // This allow for no whitespace
        // before semicolon before the end of matrix.
        write!(writer, "matrix")?;

        if !interleave {
            self.write_matrix(&mut writer)?;
        } else {
            self.write_matrix_nex_int(&mut writer);
        }

        writeln!(writer, ";")?;
        writeln!(writer, "end;")?;

        writer.flush()?;
        Ok(())
    }

    fn write_phylip(&mut self, interleave: bool) -> Result<()> {
        let mut writer = self
            .create_output_file(self.output)
            .expect("Failed writing a philip formatted file");
        write!(writer, "{} {}", self.header.ntax, self.header.nchar)?;

        if !interleave {
            self.write_matrix(&mut writer)?;
        } else {
            self.write_matrix_phy_int(&mut writer);
        }

        writer.flush()?;
        Ok(())
    }

    fn write_nex_header<W: Write>(&self, writer: &mut W, interleave: bool) -> Result<()> {
        writeln!(writer, "#NEXUS")?;
        writeln!(writer, "begin data;")?;
        writeln!(
            writer,
            "dimensions ntax={} nchar={};",
            self.header.ntax, self.header.nchar
        )?;

        write!(
            writer,
            "format datatype={} missing={} gap={}",
            self.header.datatype, self.header.missing, self.header.gap,
        )?;

        if interleave {
            write!(writer, " interleave")?;
        }

        writeln!(writer, ";")?;

        Ok(())
    }

    fn write_matrix<W: Write>(&mut self, writer: &mut W) -> Result<()> {
        // we insert newline for
        // the non-terminated new line matrix commmand.
        writeln!(writer)?;
        self.matrix.iter().for_each(|(taxa, seq)| {
            self.write_padded_seq(writer, taxa, seq)
                .expect("Failed writing nexus data matrix")
        });

        Ok(())
    }

    fn write_matrix_nex_int<W: Write>(&mut self, writer: &mut W) {
        let mat_int = self.get_matrix_int();
        mat_int.values().for_each(|seq| {
            writeln!(writer).unwrap(); // insert newline before each group.
            seq.iter().for_each(|s| {
                self.write_padded_seq(writer, &s.id, &s.seq)
                    .expect("Failed writing nexus interleave data matrix");
            });
        });
    }

    fn write_matrix_phy_int<W: Write>(&mut self, writer: &mut W) {
        let mat_int = self.get_matrix_int();
        mat_int.iter().for_each(|(idx, seq)| {
            writeln!(writer).unwrap(); // insert newline before each group.
            seq.iter().for_each(|s| match idx {
                0 => self
                    .write_padded_seq(writer, &s.id, &s.seq)
                    .expect("Failed writing phylip data matrix"),
                _ => writeln!(writer, "{}", s.seq).unwrap(),
            });
        });
    }

    fn get_matrix_int(&self) -> BTreeMap<usize, Vec<Sequence>> {
        let mut vec: BTreeMap<usize, Vec<Sequence>> = BTreeMap::new();
        let n = self.get_interleave_len();
        self.matrix.iter().for_each(|(id, seq)| {
            let chunks = self.chunk_seq(seq, n);
            chunks.iter().enumerate().for_each(|(idx, seqs)| {
                let mat = Sequence::new(id, seqs);
                match vec.get_mut(&idx) {
                    Some(value) => value.push(mat),
                    None => {
                        vec.insert(idx, vec![mat]);
                    }
                }
            })
        });

        vec
    }

    fn chunk_seq(&self, seq: &str, n: usize) -> Vec<String> {
        seq.as_bytes()
            .chunks(n)
            .map(|chunk| {
                std::str::from_utf8(chunk)
                    .expect("Failed chunking sequence")
                    .to_string()
            })
            .collect()
    }

    fn write_padded_seq<W: Write>(&mut self, writer: &mut W, taxa: &str, seq: &str) -> Result<()> {
        self.get_max_id_len();
        write!(writer, "{}", taxa)?;
        write!(writer, "{}", self.insert_whitespaces(taxa, self.id_len))?;
        writeln!(writer, "{}", seq)?;
        Ok(())
    }

    fn get_interleave_len(&self) -> usize {
        if self.header.nchar < 2000 {
            80
        } else {
            500
        }
    }

    fn get_max_id_len(&mut self) {
        self.id_len = self.matrix.keys().map(|id| id.len()).max().unwrap();
    }

    fn insert_whitespaces(&self, id: &str, max_len: usize) -> String {
        let len = id.len();
        let spaces = 1;
        if len < max_len {
            let inserts = (max_len - len) + spaces;
            " ".repeat(inserts)
        } else {
            " ".repeat(spaces)
        }
    }
}

struct Sequence {
    id: String,
    seq: String,
}

impl Sequence {
    fn new(id: &str, seq: &str) -> Self {
        Self {
            id: String::from(id),
            seq: String::from(seq),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_insert_whitespaces() {
        let max_len = 10;
        let id = "ABCDE";
        let matrix = IndexMap::new();
        let header = Header::new();
        let convert = SeqWriter::new(Path::new("."), &matrix, &header);
        assert_eq!(6, convert.insert_whitespaces(id, max_len).len())
    }

    #[test]
    fn test_chunk_seq() {
        let path = Path::new(".");
        let matrix = IndexMap::new();
        let header = Header::new();
        let convert = SeqWriter::new(path, &matrix, &header);
        let seq = "AGTCAGTC";
        let chunk = String::from("AGTC");
        let chunk2 = String::from("AGTC");
        let res = vec![chunk, chunk2];
        assert_eq!(res, convert.chunk_seq(seq, 4));
    }

    #[test]
    fn test_matrix_int() {
        let path = Path::new(".");
        let mut matrix = IndexMap::new();
        let header = Header::new();

        let id = String::from("ABC");

        let seq = String::from(
            "ATGTGTGTGTGTGTGTAAAA\
        ATGTGTGTGTGTGTGTAAAA\
        ATGTGTGTGTGTGTGTAAAA\
        ATGTGTGTGTGTGTGTAAAA\
        ATGTGTGTGTGTGTGTAAAA",
        );

        // Expected result first chunk
        let res0 = String::from(
            "ATGTGTGTGTGTGTGTAAAA\
        ATGTGTGTGTGTGTGTAAAA\
        ATGTGTGTGTGTGTGTAAAA\
        ATGTGTGTGTGTGTGTAAAA",
        );

        // Expected result second chunk
        let res1 = String::from("ATGTGTGTGTGTGTGTAAAA");

        matrix.insert(id.clone(), seq);
        let convert = SeqWriter::new(path, &matrix, &header);
        let int = convert.get_matrix_int();
        let mat_int = int.get(&0).unwrap();
        let mat_int1 = int.get(&1).unwrap();
        assert_eq!(id, mat_int[0].id.to_string());
        assert_eq!(id, mat_int1[0].id.to_string());
        assert_eq!(res0, mat_int[0].seq.to_string());
        assert_eq!(res1, mat_int1[0].seq.to_string());
    }
}
