//! Write sequences to file
use std::io::prelude::*;
use std::path::Path;

use anyhow::Result;

use crate::helper::types::{Header, OutputFmt, SeqMatrix};
use crate::writer::FileWriter;

const SPACES: [u8; 64] = [b' '; 64];

impl FileWriter for SeqWriter<'_> {}

pub struct SeqWriter<'a> {
    output: &'a Path,
    matrix: &'a SeqMatrix,
    id_len: usize,
    header: &'a Header,
}

impl<'a> SeqWriter<'a> {
    pub fn new(output: &'a Path, matrix: &'a SeqMatrix, header: &'a Header) -> Self {
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

    /// Write selected `rows` and contiguous column `runs` directly without
    /// allocating an intermediate per-partition `SeqMatrix`.
    #[cfg(feature = "align")]
    pub(crate) fn write_runs(
        &mut self,
        output_fmt: &OutputFmt,
        rows: &[usize],
        runs: &[(usize, usize)],
    ) -> Result<()> {
        match output_fmt {
            OutputFmt::Fasta => self.write_fasta_runs(rows, runs, false),
            OutputFmt::FastaInt => self.write_fasta_runs(rows, runs, true),
            OutputFmt::Nexus => self.write_nexus_runs(rows, runs, false),
            OutputFmt::NexusInt => self.write_nexus_runs(rows, runs, true),
            OutputFmt::Phylip => self.write_phylip_runs(rows, runs, false),
            OutputFmt::PhylipInt => self.write_phylip_runs(rows, runs, true),
        }
    }

    fn write_fasta(&mut self, interleave: bool) -> Result<()> {
        let mut writer = self
            .create_output_file(self.output)
            .expect("Failed writing a fasta formatted file");
        let n = self.get_interleave_len();
        for (id, seq) in self.matrix.iter() {
            writer.write_all(b">")?;
            writer.write_all(id.as_bytes())?;
            writer.write_all(b"\n")?;
            if !interleave {
                writer.write_all(seq)?;
                writer.write_all(b"\n")?;
            } else {
                for chunk in seq.chunks(n) {
                    writer.write_all(chunk)?;
                    writer.write_all(b"\n")?;
                }
            }
        }

        writer.flush()?;
        Ok(())
    }

    #[cfg(feature = "align")]
    fn write_fasta_runs(
        &mut self,
        rows: &[usize],
        runs: &[(usize, usize)],
        interleave: bool,
    ) -> Result<()> {
        let mut writer = self
            .create_output_file(self.output)
            .expect("Failed writing a fasta formatted file");
        let n = self.get_interleave_len();
        for &row in rows {
            let (id, seq) = self.matrix.get_index(row).expect("Row index out of bounds");
            writer.write_all(b">")?;
            writer.write_all(id.as_bytes())?;
            writer.write_all(b"\n")?;
            if !interleave {
                for &(start, end) in runs {
                    writer.write_all(&seq[start..end])?;
                }
                writer.write_all(b"\n")?;
            } else {
                write_wrapped_runs(&mut writer, seq, runs, n)?;
            }
        }

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
        writer.write_all(b"matrix")?;

        if !interleave {
            self.write_matrix(&mut writer)?;
        } else {
            self.write_matrix_nex_int(&mut writer);
        }

        writer.write_all(b";\nend;\n")?;

        writer.flush()?;
        Ok(())
    }

    #[cfg(feature = "align")]
    fn write_nexus_runs(
        &mut self,
        rows: &[usize],
        runs: &[(usize, usize)],
        interleave: bool,
    ) -> Result<()> {
        let mut writer = self
            .create_output_file(self.output)
            .expect("Failed writing a NEXUS formatted file");
        self.write_nex_header(&mut writer, interleave)?;
        writer.write_all(b"matrix")?;
        self.id_len = self.max_id_len_rows(rows);

        if !interleave {
            self.write_matrix_runs(&mut writer, rows, runs)?;
        } else {
            self.write_matrix_runs_int(&mut writer, rows, runs, true)?;
        }

        writer.write_all(b";\nend;\n")?;
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

    #[cfg(feature = "align")]
    fn write_phylip_runs(
        &mut self,
        rows: &[usize],
        runs: &[(usize, usize)],
        interleave: bool,
    ) -> Result<()> {
        let mut writer = self
            .create_output_file(self.output)
            .expect("Failed writing a philip formatted file");
        write!(writer, "{} {}", self.header.ntax, self.header.nchar)?;
        self.id_len = self.max_id_len_rows(rows);

        if !interleave {
            self.write_matrix_runs(&mut writer, rows, runs)?;
        } else {
            self.write_matrix_runs_int(&mut writer, rows, runs, false)?;
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
            writer.write_all(b" interleave")?;
        }

        writer.write_all(b";\n")?;

        Ok(())
    }

    fn write_matrix<W: Write>(&mut self, writer: &mut W) -> Result<()> {
        // we insert newline for
        // the non-terminated new line matrix commmand.
        writer.write_all(b"\n")?;
        self.get_max_id_len();
        for (taxa, seq) in self.matrix.iter() {
            self.write_padded_seq(writer, taxa, seq)
                .expect("Failed writing nexus data matrix");
        }

        Ok(())
    }

    #[cfg(feature = "align")]
    fn write_matrix_runs<W: Write>(
        &self,
        writer: &mut W,
        rows: &[usize],
        runs: &[(usize, usize)],
    ) -> Result<()> {
        writer.write_all(b"\n")?;
        for &row in rows {
            let (taxa, seq) = self.matrix.get_index(row).expect("Row index out of bounds");
            self.write_padded_prefix(writer, taxa)?;
            for &(start, end) in runs {
                writer.write_all(&seq[start..end])?;
            }
            writer.write_all(b"\n")?;
        }
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
            writer.write_all(b"\n").unwrap(); // insert newline before each group.
            self.matrix
                .iter()
                .filter(|(_, seq)| start < seq.len())
                .for_each(|(id, seq)| {
                    let end = (start + n).min(seq.len());
                    let chunk = &seq[start..end];
                    if idx == 0 || id_every_block {
                        self.write_padded_seq(writer, id, chunk)
                            .expect("Failed writing interleaved data matrix");
                    } else {
                        writer.write_all(chunk).unwrap();
                        writer.write_all(b"\n").unwrap();
                    }
                });
        }
    }

    #[cfg(feature = "align")]
    fn write_matrix_runs_int<W: Write>(
        &self,
        writer: &mut W,
        rows: &[usize],
        runs: &[(usize, usize)],
        id_every_block: bool,
    ) -> Result<()> {
        let n = self.get_interleave_len();
        let total_len = self.header.nchar;
        for (idx, block_start) in (0..total_len).step_by(n).enumerate() {
            writer.write_all(b"\n")?;
            let block_end = (block_start + n).min(total_len);
            for &row in rows {
                let (id, seq) = self.matrix.get_index(row).expect("Row index out of bounds");
                if idx == 0 || id_every_block {
                    self.write_padded_prefix(writer, id)?;
                }
                write_logical_window(writer, seq, runs, block_start, block_end)?;
                writer.write_all(b"\n")?;
            }
        }
        Ok(())
    }

    fn write_padded_prefix<W: Write>(&self, writer: &mut W, taxa: &str) -> Result<()> {
        writer.write_all(taxa.as_bytes())?;
        // Pad IDs to the longest ID plus one space.
        let mut pad = self.id_len.saturating_sub(taxa.len()) + 1;
        while pad > 0 {
            let take = pad.min(SPACES.len());
            writer.write_all(&SPACES[..take])?;
            pad -= take;
        }
        Ok(())
    }

    fn write_padded_seq<W: Write>(&self, writer: &mut W, taxa: &str, seq: &[u8]) -> Result<()> {
        self.write_padded_prefix(writer, taxa)?;
        writer.write_all(seq)?;
        writer.write_all(b"\n")?;
        Ok(())
    }

    fn get_interleave_len(&self) -> usize {
        if self.header.nchar < 2000 { 80 } else { 500 }
    }

    fn get_max_id_len(&mut self) {
        self.id_len = self.matrix.keys().map(|id| id.len()).max().unwrap_or(0);
    }

    #[cfg(feature = "align")]
    fn max_id_len_rows(&self, rows: &[usize]) -> usize {
        rows.iter()
            .map(|&row| {
                self.matrix
                    .get_index(row)
                    .expect("Row index out of bounds")
                    .0
                    .len()
            })
            .max()
            .unwrap_or(0)
    }
}

#[cfg(feature = "align")]
fn write_wrapped_runs<W: Write>(
    writer: &mut W,
    seq: &[u8],
    runs: &[(usize, usize)],
    width: usize,
) -> Result<()> {
    let mut col = 0usize;
    for &(start, end) in runs {
        let mut slice = &seq[start..end];
        while !slice.is_empty() {
            let take = (width - col).min(slice.len());
            writer.write_all(&slice[..take])?;
            slice = &slice[take..];
            col += take;
            if col == width {
                writer.write_all(b"\n")?;
                col = 0;
            }
        }
    }
    if col > 0 {
        writer.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(feature = "align")]
fn write_logical_window<W: Write>(
    writer: &mut W,
    seq: &[u8],
    runs: &[(usize, usize)],
    block_start: usize,
    block_end: usize,
) -> Result<()> {
    let mut pos = 0usize;
    for &(rs, re) in runs {
        let run_len = re - rs;
        let run_end = pos + run_len;
        if run_end > block_start && pos < block_end {
            let lo = block_start.max(pos);
            let hi = block_end.min(run_end);
            writer.write_all(&seq[rs + (lo - pos)..rs + (hi - pos)])?;
        }
        pos = run_end;
        if pos >= block_end {
            break;
        }
    }
    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;

    fn write_int(matrix: &SeqMatrix, id_every_block: bool) -> String {
        let header = Header::new();
        let mut convert = SeqWriter::new(Path::new("."), matrix, &header);
        convert.get_max_id_len();
        let mut out = Vec::new();
        convert.write_matrix_int(&mut out, id_every_block);
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn test_padded_seq() {
        let mut matrix = SeqMatrix::new();
        matrix.insert(String::from("ABCDEFGHIJ"), b"A");
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
        let mut matrix = SeqMatrix::new();
        let seq = b"ATGTGTGTGTGTGTGTAAAA".repeat(5);
        matrix.insert(String::from("ABC"), &seq);
        matrix.insert(String::from("ABCD"), b"GG");

        // Header nchar is 0, so the interleave length is 80.
        let res0 = "ATGTGTGTGTGTGTGTAAAA".repeat(4);
        let res1 = "ATGTGTGTGTGTGTGTAAAA";
        let nex = format!("\nABC  {}\nABCD GG\n\nABC  {}\n", res0, res1);
        let phy = format!("\nABC  {}\nABCD GG\n\n{}\n", res0, res1);
        assert_eq!(nex, write_int(&matrix, true));
        assert_eq!(phy, write_int(&matrix, false));
    }
}
