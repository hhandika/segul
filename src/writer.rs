use std::fs::{self, File};
use std::io::prelude::*;
use std::io::{self, BufWriter, LineWriter, Result};
use std::iter;
use std::path::{Path, PathBuf};

use crate::common::{Header, OutputFormat, Partition, PartitionFormat};
use indexmap::IndexMap;

pub struct SeqWriter<'a> {
    path: &'a Path,
    output: PathBuf,
    matrix: &'a IndexMap<String, String>,
    id_len: usize,
    header: Header,
    partition: Option<&'a [Partition]>,
    part_format: &'a PartitionFormat,
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
        }
    }

    pub fn write_sequence(&mut self, filetype: &OutputFormat) {
        self.get_output_name(filetype);
        self.get_ntax();
        self.get_nchar();
        self.get_max_id_len();

        match filetype {
            OutputFormat::Nexus => self.write_nexus().expect("CANNOT WRITE A NEXUS FILE."),
            OutputFormat::Phylip => self.write_phylip().expect("CANNOT WRITE A PHYLIP FILE."),
            _ => self.write_fasta(),
        }
    }

    pub fn write_fasta(&mut self) {
        self.get_output_name(&OutputFormat::Fasta);
        let mut writer = self.create_output_file(&self.output);
        self.matrix.iter().for_each(|(id, seq)| {
            writeln!(writer, ">{}", id).unwrap();
            writeln!(writer, "{}", seq).unwrap();
        });
    }

    pub fn display_save_path(&self) {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "Output\t\t: {}\n", self.output.display()).unwrap();
    }

    fn write_phylip(&mut self) -> Result<()> {
        let mut writer = self.create_output_file(&self.output);
        writeln!(
            writer,
            "{} {}",
            self.header.ntax.as_ref().unwrap(),
            self.header.nchar.as_ref().unwrap()
        )?;
        self.write_matrix(&mut writer);
        if self.partition.is_some() {
            self.write_partition_sep();
        }
        Ok(())
    }

    fn write_nexus(&mut self) -> Result<()> {
        self.get_datatype();
        self.get_missing();
        self.get_gap();
        let mut writer = self.create_output_file(&self.output);
        writeln!(writer, "#NEXUS")?;
        writeln!(writer, "begin data;")?;
        writeln!(
            writer,
            "dimensions ntax={} nchar={};",
            self.header.ntax.as_ref().unwrap(),
            self.header.nchar.as_ref().unwrap()
        )?;
        writeln!(
            writer,
            "format datatype={} missing={} gap={};",
            self.header.datatype.as_ref().unwrap(),
            self.header.missing.as_ref().unwrap(),
            self.header.gap.as_ref().unwrap()
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
            self.check_sequence_len(seq.len(), taxa);
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
        let fname = format!(
            "{}_partition.txt",
            self.output.file_stem().unwrap().to_string_lossy()
        );
        let mut writer = self.create_output_file(Path::new(&fname));
        match &self.partition {
            Some(partition) => partition.iter().for_each(|part| {
                writeln!(writer, "DNA, {} = {}-{}", part.gene, part.start, part.end).unwrap();
            }),
            None => eprintln!("CANNOT READ PARTITION DATA"),
        }
    }

    fn write_part_nexus_sep(&self) {
        let fname = format!(
            "{}_partition.nex",
            self.output.file_stem().unwrap().to_string_lossy()
        );
        let mut writer = self.create_output_file(Path::new(&fname));
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

    fn check_sequence_len(&self, len: usize, taxa: &str) {
        if len != *self.header.nchar.as_ref().unwrap() {
            panic!(
                "DIFFERENT SEQUENCE LENGTH FOUND AT {}. \
                MAKE SURE THE INPUT IS AN ALIGNMENT",
                taxa
            )
        }
    }

    fn get_output_name(&mut self, ext: &OutputFormat) {
        self.output = match ext {
            OutputFormat::Fasta => self.path.with_extension("fas"),
            OutputFormat::Nexus => self.path.with_extension("nex"),
            OutputFormat::Phylip => self.path.with_extension("phy"),
        };
    }

    fn create_output_file(&self, fname: &Path) -> LineWriter<File> {
        fs::create_dir_all(fname.parent().unwrap()).expect("CANNOT CREATE A TARGET DIRECTORY");
        let file = File::create(&fname).expect("CANNOT CREATE OUTPUT FILE");
        LineWriter::new(file)
    }

    fn get_ntax(&mut self) {
        if self.header.ntax.is_none() {
            self.header.ntax = Some(self.matrix.len());
        }
    }

    fn get_nchar(&mut self) {
        if self.header.nchar.is_none() {
            let (_, chars) = self.matrix.iter().next().unwrap();
            self.header.nchar = Some(chars.len())
        }
    }

    fn get_datatype(&mut self) {
        if self.header.datatype.is_none() {
            self.header.datatype = Some(String::from("dna"));
        }
    }

    fn get_missing(&mut self) {
        if self.header.missing.is_none() {
            self.header.missing = Some('?');
        }
    }

    fn get_gap(&mut self) {
        if self.header.gap.is_none() {
            self.header.gap = Some('-');
        }
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
        convert.get_output_name(&OutputFormat::Fasta);
        assert_eq!(output, convert.output);
    }
}
