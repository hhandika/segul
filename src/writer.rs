use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::prelude::*;
use std::io::{self, BufWriter, LineWriter, Result};
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

    pub fn write_sequence(&mut self, output_format: &SeqFormat) {
        self.get_output_name(output_format);
        self.get_max_id_len();

        if self.partition.is_some() {
            self.get_partition_path();
        }

        match output_format {
            SeqFormat::Nexus => self.write_nexus(false).expect("CANNOT WRITE A NEXUS FILE."),
            SeqFormat::NexusInt => self
                .write_nexus(true)
                .expect("CANNOT WRITE A NEXUS INTERLEAVE FILE."),
            SeqFormat::Phylip => self
                .write_phylip(false)
                .expect("CANNOT WRITE A PHYLIP FILE."),
            SeqFormat::PhylipInt => self
                .write_phylip(true)
                .expect("CANNOT WRITE A PHYLIP FILE."),
            SeqFormat::Fasta => self.write_fasta(false),
            SeqFormat::FastaInt => self.write_fasta(true),
        }
    }

    fn write_fasta(&mut self, interleave: bool) {
        let mut writer = self.create_output_file(&self.output);
        let n = 500;
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

    fn write_phylip(&mut self, interleave: bool) -> Result<()> {
        let mut writer = self.create_output_file(&self.output);
        writeln!(writer, "{} {}", self.header.ntax, self.header.nchar)?;

        if !interleave {
            self.write_matrix(&mut writer);
        } else {
            self.write_matrix_phy_int(&mut writer);
        }

        if self.partition.is_some() {
            self.write_partition_sep();
        }

        Ok(())
    }

    fn write_nexus(&mut self, interleave: bool) -> Result<()> {
        let mut writer = self.create_output_file(&self.output);
        self.write_nex_header(&mut writer)?;
        self.write_nex_format(&mut writer, interleave)?;
        writeln!(writer, "matrix")?;

        if !interleave {
            self.write_matrix(&mut writer);
        } else {
            self.write_matrix_nex_int(&mut writer);
        }

        writeln!(writer, ";")?;
        writeln!(writer, "end;")?;

        if self.partition.is_some() {
            match self.part_format {
                PartitionFormat::Nexus => self
                    .write_part_nexus(&mut writer)
                    .expect("CANNOT WRITER NEXUS PARTITION"),
                PartitionFormat::Raxml => self.write_part_phylip(),
                _ => self.write_part_nexus_sep(),
            }
        }

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

    fn write_matrix<W: Write>(&self, writer: &mut W) {
        self.matrix.iter().for_each(|(taxa, seq)| {
            self.write_padded_seq(writer, taxa, seq)
                .expect("CANNOT WRITE SEQ MATRIX")
        });
    }

    fn write_matrix_phy_int<W: Write>(&self, writer: &mut W) {
        let mat_int = self.get_matrix_int();
        mat_int.iter().for_each(|(idx, seq)| {
            seq.iter().for_each(|s| match idx {
                0 => self
                    .write_padded_seq(writer, &s.id, &s.seq)
                    .expect("CANNOT WRITE PADDED SEQ MATRIX"),
                _ => writeln!(writer, "{}", s.seq).unwrap(),
            });

            writeln!(writer).unwrap();
        });
    }

    fn write_matrix_nex_int<W: Write>(&self, writer: &mut W) {
        let mat_int = self.get_matrix_int();
        mat_int.values().for_each(|seq| {
            seq.iter().for_each(|s| {
                self.write_padded_seq(writer, &s.id, &s.seq)
                    .expect("CANNOT WRITE PADDED SEQ MATRIX");
            });
            writeln!(writer).unwrap();
        });
    }

    fn write_padded_seq<W: Write>(&self, writer: &mut W, taxa: &str, seq: &str) -> Result<()> {
        write!(writer, "{}", taxa)?;
        write!(writer, "{}", self.insert_whitespaces(taxa, self.id_len))?;
        writeln!(writer, "{}", seq)?;
        Ok(())
    }

    fn write_partition_sep(&self) {
        match self.part_format {
            PartitionFormat::Charset => self.write_part_nexus_sep(),
            PartitionFormat::Raxml => self.write_part_phylip(),
            _ => eprintln!("UNKNOWN PARTITION FORMAT"),
        }
    }

    fn write_part_phylip(&self) {
        let mut writer = self.create_output_file(Path::new(&self.part_file));
        match &self.partition {
            Some(partition) => partition.iter().for_each(|part| {
                writeln!(writer, "DNA, {} = {}-{}", part.gene, part.start, part.end).unwrap();
            }),
            None => eprintln!("CANNOT READ PARTITION DATA"),
        }
    }

    fn write_part_nexus_sep(&self) {
        let mut writer = self.create_output_file(&self.part_file);
        writeln!(writer, "#nexus").unwrap();
        self.write_part_nexus(&mut writer)
            .expect("CANNOT WRITE NEXUS PARTITION");
    }

    fn write_part_nexus<W: Write>(&self, writer: &mut W) -> Result<()> {
        writeln!(writer, "begin sets;")?;
        match &self.partition {
            Some(partition) => partition.iter().for_each(|part| {
                writeln!(
                    writer,
                    "charset {} = {}-{};",
                    part.gene, part.start, part.end
                )
                .unwrap();
            }),
            None => panic!("CANNOT READ PARTITION DATA"),
        }
        writeln!(writer, "end;")?;
        Ok(())
    }

    fn get_partition_path(&mut self) {
        match self.part_format {
            PartitionFormat::Charset => {
                self.part_file = self
                    .output
                    .parent()
                    .unwrap()
                    .join(&self.get_partition_name("nex"));
            }
            PartitionFormat::Raxml => {
                self.part_file = self
                    .output
                    .parent()
                    .unwrap()
                    .join(&self.get_partition_name("txt"));
            }
            PartitionFormat::Nexus => self.part_file = PathBuf::from("in-file"),
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

    fn create_output_file(&self, fname: &Path) -> LineWriter<File> {
        fs::create_dir_all(fname.parent().unwrap()).expect("CANNOT CREATE A TARGET DIRECTORY");
        let file = File::create(&fname).expect("CANNOT CREATE OUTPUT FILE");
        LineWriter::new(file)
    }

    fn get_matrix_int(&self) -> BTreeMap<usize, Vec<Sequence>> {
        let mut vec: BTreeMap<usize, Vec<Sequence>> = BTreeMap::new();
        let n = 500;
        self.matrix.iter().for_each(|(id, seq)| {
            let chunks = self.chunk_seq(seq, n);
            chunks.iter().enumerate().for_each(|(idx, seqs)| {
                let mat = Sequence::new(id, seqs);
                if vec.contains_key(&idx) {
                    if let Some(value) = vec.get_mut(&idx) {
                        value.push(mat);
                    }
                } else {
                    vec.insert(idx, vec![mat]);
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
                    .expect("FAILED CHUNKING THE SEQ OUTPUT")
                    .to_string()
            })
            .collect()
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
