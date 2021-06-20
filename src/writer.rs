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
            SeqFormat::Nexus => self.write_nexus().expect("CANNOT WRITE A NEXUS FILE."),
            SeqFormat::Phylip => self.write_phylip().expect("CANNOT WRITE A PHYLIP FILE."),
            _ => self.write_fasta(),
        }
    }

    pub fn write_fasta(&mut self) {
        self.get_output_name(&SeqFormat::Fasta);
        let mut writer = self.create_output_file(&self.output);
        self.matrix.iter().for_each(|(id, seq)| {
            writeln!(writer, ">{}", id).unwrap();
            writeln!(writer, "{}", seq).unwrap();
        });
    }

    pub fn display_save_path(&self) {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "Output\t\t: {}", self.output.display()).unwrap();
    }

    pub fn display_partition_path(&self) {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "Partition\t: {}\n", &self.part_file.display()).unwrap();
    }

    fn write_phylip(&mut self) -> Result<()> {
        let mut writer = self.create_output_file(&self.output);
        writeln!(writer, "{} {}", self.header.ntax, self.header.nchar)?;
        self.write_matrix(&mut writer);
        if self.partition.is_some() {
            self.write_partition_sep();
        }
        Ok(())
    }

    fn write_nexus(&mut self) -> Result<()> {
        let mut writer = self.create_output_file(&self.output);
        writeln!(writer, "#NEXUS")?;
        writeln!(writer, "begin data;")?;
        writeln!(
            writer,
            "dimensions ntax={} nchar={};",
            self.header.ntax, self.header.nchar
        )?;
        writeln!(
            writer,
            "format datatype={} missing={} gap={};",
            self.header.datatype, self.header.missing, self.header.gap,
        )?;
        writeln!(writer, "matrix")?;
        self.write_matrix(&mut writer);
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

    fn write_matrix<W: Write>(&self, writer: &mut W) {
        self.matrix.iter().for_each(|(taxa, seq)| {
            write!(writer, "{}", taxa).unwrap();
            write!(writer, "{}", self.insert_whitespaces(taxa, self.id_len)).unwrap();
            writeln!(writer, "{}", seq).unwrap();
        });
    }

    fn write_partition_sep(&self) {
        match self.part_format {
            PartitionFormat::NexusSeparate => self.write_part_nexus_sep(),
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
            PartitionFormat::NexusSeparate => {
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
            SeqFormat::Fasta => self.path.with_extension("fas"),
            SeqFormat::Nexus => self.path.with_extension("nex"),
            SeqFormat::Phylip => self.path.with_extension("phy"),
        };
    }

    fn create_output_file(&self, fname: &Path) -> LineWriter<File> {
        fs::create_dir_all(fname.parent().unwrap()).expect("CANNOT CREATE A TARGET DIRECTORY");
        let file = File::create(&fname).expect("CANNOT CREATE OUTPUT FILE");
        LineWriter::new(file)
    }

    fn get_max_id_len(&mut self) {
        let longest = self.matrix.keys().max_by_key(|id| id.len()).unwrap();
        self.id_len = longest.len();
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
