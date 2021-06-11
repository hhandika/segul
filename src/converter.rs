use std::collections::BTreeMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::LineWriter;
use std::iter;
use std::path::{Path, PathBuf};

pub struct Converter<'m> {
    path: PathBuf,
    outname: PathBuf,
    matrix: &'m BTreeMap<String, String>,
}

impl<'m> Converter<'m> {
    pub fn new(path: &str, matrix: &'m BTreeMap<String, String>) -> Self {
        Self {
            path: PathBuf::from(path),
            outname: PathBuf::new(),
            matrix,
        }
    }

    pub fn write_fasta(&mut self) {
        self.get_output_name("fas");
        let mut writer = self.create_output_file(&self.outname);
        self.matrix.iter().for_each(|(id, seq)| {
            writeln!(writer, ">{}", id).unwrap();
            writeln!(writer, "{}", seq).unwrap();
        });
        println!("Save as {}", self.outname.display());
    }

    pub fn write_phylip(&mut self) {
        self.get_output_name("phy");
        let (num_taxa, chars) = self.get_num_taxa_chars();
        let max_id_len = self.get_max_id_len();
        let mut writer = self.create_output_file(&self.outname);
        writeln!(writer, "{} {}", num_taxa, chars).unwrap();

        self.matrix.iter().for_each(|(taxa, seq)| {
            if seq.len() != chars {
                panic!(
                    "DIFFERENT SEQUENCE LENGTH FOUND AT {}. \
                MAKE SURE THE INPUT IS AN ALIGNMENT",
                    taxa
                )
            }
            write!(writer, "{}", taxa).unwrap();
            write!(writer, "{}", self.insert_whitespaces(taxa, max_id_len)).unwrap();
            writeln!(writer, "{}", seq).unwrap();
        });

        println!("Save as {}", self.outname.display());
    }

    fn get_output_name(&mut self, ext: &str) {
        let name = Path::new(self.path.file_name().unwrap());
        self.outname = name.with_extension(ext);
    }

    fn create_output_file(&self, fname: &Path) -> LineWriter<File> {
        let file = File::create(&fname).expect("CANNOT CREATE OUTPUT FILE");
        LineWriter::new(file)
    }

    fn get_num_taxa_chars(&self) -> (usize, usize) {
        let taxa = self.matrix.len();
        let (_, seq) = self.matrix.iter().next().unwrap();
        let chars: usize = seq.len();
        (taxa, chars)
    }

    fn get_max_id_len(&self) -> usize {
        let longest = self.matrix.keys().max_by_key(|id| id.len()).unwrap();
        longest.len()
    }

    fn insert_whitespaces(&self, id: &str, max_len: usize) -> String {
        let len = id.len();
        let spaces = 5;
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
        let matrix = BTreeMap::new();
        let convert = Converter::new(".", &matrix);
        assert_eq!(10, convert.insert_whitespaces(id, max_len).len())
    }
}
