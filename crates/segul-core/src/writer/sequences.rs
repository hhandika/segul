//! Write sequences to file
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
                seq.as_bytes().chunks(n).for_each(|chunk| {
                    writer.write_all(chunk).unwrap();
                    writeln!(writer).unwrap();
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
        self.get_max_id_len();
        self.matrix.iter().for_each(|(taxa, seq)| {
            self.write_padded_seq(writer, taxa, seq.as_bytes())
                .expect("Failed writing nexus data matrix")
        });

        Ok(())
    }

    fn write_matrix_nex_int<W: Write>(&mut self, writer: &mut W) {
        self.get_max_id_len();
        self.write_matrix_int(writer, true);
    }

    fn write_matrix_phy_int<W: Write>(&mut self, writer: &mut W) {
        self.get_max_id_len();
        self.write_matrix_int(writer, false);
    }

    /// Write interleaved blocks by slicing each sequence in place,
    /// avoiding copying the matrix into chunks.
    /// For phylip, the IDs are only written in the first block.
    fn write_matrix_int<W: Write>(&self, writer: &mut W, id_every_block: bool) {
        let n = self.get_interleave_len();
        let longest = self.matrix.values().map(|s| s.len()).max().unwrap_or(0);
        for (idx, start) in (0..longest).step_by(n).enumerate() {
            writeln!(writer).unwrap(); // insert newline before each group.
            self.matrix
                .iter()
                .filter(|(_, seq)| start < seq.len())
                .for_each(|(id, seq)| {
                    let end = (start + n).min(seq.len());
                    let chunk = &seq.as_bytes()[start..end];
                    if idx == 0 || id_every_block {
                        self.write_padded_seq(writer, id, chunk)
                            .expect("Failed writing interleaved data matrix");
                    } else {
                        writer.write_all(chunk).unwrap();
                        writeln!(writer).unwrap();
                    }
                });
        }
    }

    fn write_padded_seq<W: Write>(&self, writer: &mut W, taxa: &str, seq: &[u8]) -> Result<()> {
        // Pad IDs to the longest ID plus one space.
        write!(writer, "{:<width$} ", taxa, width = self.id_len)?;
        writer.write_all(seq)?;
        writeln!(writer)?;
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
}

#[cfg(test)]
mod test {
    use super::*;

    fn write_int(matrix: &IndexMap<String, String>, id_every_block: bool) -> String {
        let header = Header::new();
        let mut convert = SeqWriter::new(Path::new("."), matrix, &header);
        convert.get_max_id_len();
        let mut out = Vec::new();
        convert.write_matrix_int(&mut out, id_every_block);
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn test_padded_seq() {
        let mut matrix = IndexMap::new();
        matrix.insert(String::from("ABCDEFGHIJ"), String::from("A"));
        let header = Header::new();
        let mut convert = SeqWriter::new(Path::new("."), &matrix, &header);
        convert.get_max_id_len();
        let mut out = Vec::new();
        convert
            .write_padded_seq(&mut out, "ABCDE", b"AGTC")
            .unwrap();
        assert_eq!("ABCDE      AGTC\n", String::from_utf8(out).unwrap());
    }

    #[test]
    fn test_matrix_int() {
        let mut matrix = IndexMap::new();
        let seq = "ATGTGTGTGTGTGTGTAAAA".repeat(5);
        matrix.insert(String::from("ABC"), seq);
        matrix.insert(String::from("ABCD"), String::from("GG"));

        // Header nchar is 0, so the interleave length is 80.
        let res0 = "ATGTGTGTGTGTGTGTAAAA".repeat(4);
        let res1 = "ATGTGTGTGTGTGTGTAAAA";
        let nex = format!("\nABC  {}\nABCD GG\n\nABC  {}\n", res0, res1);
        let phy = format!("\nABC  {}\nABCD GG\n\n{}\n", res0, res1);
        assert_eq!(nex, write_int(&matrix, true));
        assert_eq!(phy, write_int(&matrix, false));
    }
}
