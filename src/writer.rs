use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::prelude::*;
use std::io::{self, BufWriter, Result};
use std::iter;
use std::path::{Path, PathBuf};

use crate::common::{Header, Partition, PartitionFormat, SeqFormat};
use indexmap::IndexMap;

pub struct SeqWriter<'a> {
    path: &'a Path,
    output: PathBuf,
    matrix: &'a IndexMap<String, String>,
    id_len: usize,
    header: Header,
    partition: Option<&'a [Partition]>,
    part_format: &'a PartitionFormat,
    part_file: PathBuf,
}

impl<'a> SeqWriter<'a> {
    pub fn new(
        path: &'a Path,
        matrix: &'a IndexMap<String, String>,
        header: Header,
        partition: Option<&'a [Partition]>,
        part_format: &'a PartitionFormat,
    ) -> Self {
        Self {
            path,
            output: PathBuf::new(),
            id_len: 0,
            matrix,
            header,
            partition,
            part_format,
            part_file: PathBuf::new(),
        }
    }

    pub fn write_sequence(&mut self, output_format: &SeqFormat) -> Result<()> {
        self.get_output_name(output_format);

        if self.partition.is_some() {
            self.get_partition_path();
        }

        match output_format {
            SeqFormat::Nexus => self.write_nexus(false)?,
            SeqFormat::NexusInt => self.write_nexus(true)?,
            SeqFormat::Phylip => self.write_phylip(false)?,
            SeqFormat::PhylipInt => self.write_phylip(true)?,
            SeqFormat::Fasta => self.write_fasta(false)?,
            SeqFormat::FastaInt => self.write_fasta(true)?,
        }

        Ok(())
    }

    pub fn display_save_path(&self) {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "Output\t\t: {}", self.output.display()).unwrap();
    }

    pub fn display_partition_path(&self) {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "Partition\t: {}", &self.part_file.display()).unwrap();
    }

    fn write_fasta(&mut self, interleave: bool) -> Result<()> {
        let mut writer = self.create_output_file(&self.output);
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

        if self.partition.is_some() {
            self.write_partition_sep();
        }

        writer.flush()?;
        Ok(())
    }

    fn write_nexus(&mut self, interleave: bool) -> Result<()> {
        let mut writer = self.create_output_file(&self.output);
        self.write_nex_header(&mut writer)?;
        self.write_nex_format(&mut writer, interleave)?;

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

        if self.partition.is_some() {
            match self.part_format {
                PartitionFormat::Nexus => self
                    .write_part_nexus(&mut writer, false)
                    .expect("CANNOT WRITER NEXUS PARTITION"),
                PartitionFormat::NexusCodon => self
                    .write_part_nexus(&mut writer, true)
                    .expect("CANNOT WRITER NEXUS PARTITION"),
                _ => self.write_partition_sep(),
            }
        }

        writer.flush()?;
        Ok(())
    }

    fn write_phylip(&mut self, interleave: bool) -> Result<()> {
        let mut writer = self.create_output_file(&self.output);
        write!(writer, "{} {}", self.header.ntax, self.header.nchar)?;

        if !interleave {
            self.write_matrix(&mut writer)?;
        } else {
            self.write_matrix_phy_int(&mut writer);
        }

        if self.partition.is_some() {
            self.write_partition_sep();
        }

        writer.flush()?;
        Ok(())
    }

    fn write_nex_header<W: Write>(&self, writer: &mut W) -> Result<()> {
        writeln!(writer, "#NEXUS")?;
        writeln!(writer, "begin data;")?;
        writeln!(
            writer,
            "dimensions ntax={} nchar={};",
            self.header.ntax, self.header.nchar
        )?;

        Ok(())
    }

    fn write_nex_format<W: Write>(&self, writer: &mut W, interleave: bool) -> Result<()> {
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
                .expect("CANNOT WRITE SEQ MATRIX")
        });

        Ok(())
    }

    fn write_matrix_nex_int<W: Write>(&mut self, writer: &mut W) {
        let mat_int = self.get_matrix_int();
        mat_int.values().for_each(|seq| {
            writeln!(writer).unwrap(); // insert newline before each group.
            seq.iter().for_each(|s| {
                self.write_padded_seq(writer, &s.id, &s.seq)
                    .expect("CANNOT WRITE PADDED SEQ MATRIX");
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
                    .expect("CANNOT WRITE PADDED SEQ MATRIX"),
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

    fn write_padded_seq<W: Write>(&mut self, writer: &mut W, taxa: &str, seq: &str) -> Result<()> {
        self.get_max_id_len();
        write!(writer, "{}", taxa)?;
        write!(writer, "{}", self.insert_whitespaces(taxa, self.id_len))?;
        writeln!(writer, "{}", seq)?;
        Ok(())
    }

    fn write_partition_sep(&self) {
        match self.part_format {
            PartitionFormat::Charset => self.write_part_nexus_sep(false),
            PartitionFormat::CharsetCodon => self.write_part_nexus_sep(true),
            PartitionFormat::Raxml => self.write_part_raxml(false),
            PartitionFormat::RaxmlCodon => self.write_part_raxml(true),
            _ => eprintln!("UNKNOWN PARTITION FORMAT"),
        }
    }

    fn write_part_raxml(&self, codon: bool) {
        let mut writer = self.create_output_file(Path::new(&self.part_file));
        match &self.partition {
            Some(partition) => partition.iter().for_each(|part| {
                if codon {
                    self.write_raxml_codon(&mut writer, part).unwrap();
                } else {
                    writeln!(writer, "DNA, {} = {}-{}", part.gene, part.start, part.end).unwrap();
                }
            }),
            None => eprintln!("CANNOT FIND PARTITION DATA"),
        }
    }

    fn write_part_nexus_sep(&self, codon: bool) {
        let mut writer = self.create_output_file(&self.part_file);
        writeln!(writer, "#nexus").unwrap();
        self.write_part_nexus(&mut writer, codon)
            .expect("CANNOT WRITE NEXUS PARTITION");
    }

    fn write_part_nexus<W: Write>(&self, writer: &mut W, codon: bool) -> Result<()> {
        writeln!(writer, "begin sets;")?;
        match &self.partition {
            Some(partition) => partition.iter().for_each(|part| {
                if codon {
                    self.write_nex_codon(writer, &part).unwrap();
                } else {
                    writeln!(
                        writer,
                        "charset {} = {}-{};",
                        part.gene, part.start, part.end
                    )
                    .unwrap();
                }
            }),
            None => panic!("CANNOT READ PARTITION DATA"),
        }
        writeln!(writer, "end;")?;
        Ok(())
    }

    fn write_raxml_codon<W: Write>(&self, writer: &mut W, part: &Partition) -> Result<()> {
        writeln!(
            writer,
            "DNA, {}-Subset1 = {}-{}\\3",
            part.gene, part.start, part.end
        )?;
        writeln!(
            writer,
            "DNA, {}-Subset2 = {}-{}\\3",
            part.gene,
            part.start + 1,
            part.end
        )?;
        writeln!(
            writer,
            "DNA, {}-Subset3 = {}-{}\\3",
            part.gene,
            part.start + 2,
            part.end
        )?;

        Ok(())
    }

    fn write_nex_codon<W: Write>(&self, writer: &mut W, part: &Partition) -> Result<()> {
        writeln!(
            writer,
            "charset {}-Subset1 = {}-{}\\3;",
            part.gene, part.start, part.end
        )?;
        writeln!(
            writer,
            "charset {}-Subset2 = {}-{}\\3;",
            part.gene,
            part.start + 1,
            part.end
        )?;
        writeln!(
            writer,
            "charset {}-Subset3 = {}-{}\\3;",
            part.gene,
            part.start + 2,
            part.end
        )?;

        Ok(())
    }

    fn get_partition_path(&mut self) {
        match self.part_format {
            PartitionFormat::Charset | PartitionFormat::CharsetCodon => {
                self.part_file = self
                    .output
                    .parent()
                    .unwrap()
                    .join(&self.get_partition_name("nex"));
            }
            PartitionFormat::Raxml | PartitionFormat::RaxmlCodon => {
                self.part_file = self
                    .output
                    .parent()
                    .unwrap()
                    .join(&self.get_partition_name("txt"));
            }
            PartitionFormat::Nexus | PartitionFormat::NexusCodon => {
                self.part_file = PathBuf::from("in-file")
            }
            _ => (),
        }
    }

    fn get_partition_name(&self, ext: &str) -> PathBuf {
        let fname = format!(
            "{}_partition.{}",
            self.output.file_stem().unwrap().to_string_lossy(),
            ext
        );

        PathBuf::from(fname)
    }

    fn get_output_name(&mut self, ext: &SeqFormat) {
        self.output = match ext {
            SeqFormat::Fasta | SeqFormat::FastaInt => self.path.with_extension("fas"),
            SeqFormat::Nexus | SeqFormat::NexusInt => self.path.with_extension("nex"),
            SeqFormat::Phylip | SeqFormat::PhylipInt => self.path.with_extension("phy"),
        };
    }

    fn create_output_file(&self, fname: &Path) -> BufWriter<File> {
        fs::create_dir_all(fname.parent().unwrap()).expect("CANNOT CREATE A TARGET DIRECTORY");
        let file = File::create(&fname).expect("CANNOT CREATE OUTPUT FILE");
        BufWriter::new(file)
    }

    fn chunk_seq(&self, seq: &str, n: usize) -> Vec<String> {
        seq.as_bytes()
            .chunks(n)
            .map(|chunk| {
                std::str::from_utf8(chunk)
                    .expect("FAILED CHUNKING THE SEQ OUTPUT")
                    .to_string()
            })
            .collect()
    }

    fn get_interleave_len(&self) -> usize {
        if self.header.nchar < 1000 {
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
            iter::repeat(' ').take(inserts).collect()
        } else {
            iter::repeat(' ').take(spaces).collect()
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
    fn insert_whitespaces_test() {
        let max_len = 10;
        let id = "ABCDE";
        let matrix = IndexMap::new();
        let header = Header::new();
        let convert = SeqWriter::new(
            Path::new("."),
            &matrix,
            header,
            None,
            &PartitionFormat::None,
        );
        assert_eq!(6, convert.insert_whitespaces(id, max_len).len())
    }

    #[test]
    fn get_output_name_test() {
        let path = Path::new("sanger/cytb");
        let matrix = IndexMap::new();
        let header = Header::new();
        let mut convert = SeqWriter::new(path, &matrix, header, None, &PartitionFormat::None);
        let output = PathBuf::from("sanger/cytb.fas");
        convert.get_output_name(&SeqFormat::Fasta);
        assert_eq!(output, convert.output);
    }
}
