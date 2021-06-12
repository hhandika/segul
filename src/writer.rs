use std::collections::BTreeMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::{LineWriter, Result};
use std::iter;
use std::path::{Path, PathBuf};

use crate::common::SeqFormat;

pub struct SeqWriter<'m> {
    path: PathBuf,
    outname: PathBuf,
    matrix: &'m BTreeMap<String, String>,
    id_len: usize,
    ntax: Option<usize>,
    nchar: Option<usize>,
    datatype: Option<String>,
    missing: Option<char>,
    gap: Option<char>,
}

impl<'m> SeqWriter<'m> {
    pub fn new(
        path: &str,
        matrix: &'m BTreeMap<String, String>,
        ntax: Option<usize>,
        nchar: Option<usize>,
        datatype: Option<String>,
        missing: Option<char>,
        gap: Option<char>,
    ) -> Self {
        Self {
            path: PathBuf::from(path),
            outname: PathBuf::new(),
            id_len: 0,
            matrix,
            ntax,
            nchar,
            datatype,
            missing,
            gap,
        }
    }

    pub fn write_sequence(&mut self, filetype: &SeqFormat) {
        self.get_output_name(filetype);
        self.get_ntax();
        self.get_nchar();
        self.get_max_id_len();

        match filetype {
            SeqFormat::Nexus => self.write_nexus().expect("CANNOT WRITE A NEXUS FILE."),
            SeqFormat::Phylip => self.write_phylip().expect("CANNOT WRITE A PHYLIP FILE."),
            _ => self.write_fasta(),
        }
    }

    pub fn write_fasta(&mut self) {
        self.get_output_name(&SeqFormat::Fasta);
        let mut writer = self.create_output_file(&self.outname);
        self.matrix.iter().for_each(|(id, seq)| {
            writeln!(writer, ">{}", id).unwrap();
            writeln!(writer, "{}", seq).unwrap();
        });
        self.display_save_path();
    }

    fn write_phylip(&mut self) -> Result<()> {
        let mut writer = self.create_output_file(&self.outname);
        writeln!(
            writer,
            "{} {}",
            self.ntax.as_ref().unwrap(),
            self.nchar.as_ref().unwrap()
        )?;
        self.write_matrix(&mut writer);
        self.display_save_path();
        Ok(())
    }

    fn write_nexus(&mut self) -> Result<()> {
        self.get_datatype();
        self.get_missing();
        self.get_gap();
        let mut writer = self.create_output_file(&self.outname);
        writeln!(writer, "#NEXUS")?;
        writeln!(writer, "begin data;")?;
        writeln!(
            writer,
            "dimensions ntax={} nchar={};",
            self.ntax.as_ref().unwrap(),
            self.nchar.as_ref().unwrap()
        )?;
        writeln!(
            writer,
            "format datatype={} missing={} gap={};",
            self.datatype.as_ref().unwrap(),
            self.missing.as_ref().unwrap(),
            self.gap.as_ref().unwrap()
        )?;
        writeln!(writer, "matrix")?;
        self.write_matrix(&mut writer);
        writeln!(writer, ";")?;
        writeln!(writer, "end;")?;
        self.display_save_path();
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

    fn check_sequence_len(&self, len: usize, taxa: &str) {
        if len != *self.nchar.as_ref().unwrap() {
            panic!(
                "DIFFERENT SEQUENCE LENGTH FOUND AT {}. \
                MAKE SURE THE INPUT IS AN ALIGNMENT",
                taxa
            )
        }
    }
    fn display_save_path(&self) {
        println!("Save as {}", self.outname.display());
    }

    fn get_output_name(&mut self, ext: &SeqFormat) {
        let name = Path::new(self.path.file_name().unwrap());
        match ext {
            SeqFormat::Fasta => self.outname = name.with_extension("fas"),
            SeqFormat::Nexus => self.outname = name.with_extension("nex"),
            SeqFormat::Phylip => self.outname = name.with_extension("phy"),
        };
    }

    fn create_output_file(&self, fname: &Path) -> LineWriter<File> {
        let file = File::create(&fname).expect("CANNOT CREATE OUTPUT FILE");
        LineWriter::new(file)
    }

    fn get_ntax(&mut self) {
        if self.ntax.is_none() {
            self.ntax = Some(self.matrix.len());
        }
    }

    fn get_nchar(&mut self) {
        if self.nchar.is_none() {
            let (_, chars) = self.matrix.iter().next().unwrap();
            self.nchar = Some(chars.len())
        }
    }

    fn get_datatype(&mut self) {
        if self.datatype.is_none() {
            self.datatype = Some(String::from("dna"));
        }
    }

    fn get_missing(&mut self) {
        if self.missing.is_none() {
            self.missing = Some('?');
        }
    }

    fn get_gap(&mut self) {
        if self.gap.is_none() {
            self.gap = Some('-');
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
            iter::repeat(' ').take(inserts).collect::<String>()
        } else {
            iter::repeat(' ').take(spaces).collect::<String>()
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
        let ntax = Some(2);
        let nchar = Some(5);
        let datatype = Some(String::from("dna"));
        let missing = Some('?');
        let gap = Some('-');
        let matrix = BTreeMap::new();
        let convert = SeqWriter::new(".", &matrix, ntax, nchar, datatype, missing, gap);
        assert_eq!(6, convert.insert_whitespaces(id, max_len).len())
    }
}
