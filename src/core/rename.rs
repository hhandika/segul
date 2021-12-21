// use std::fs;
use std::path::{Path, PathBuf};

#[allow(unused_imports)]
use crate::helper::sequence::{SeqCheck, Sequence};
#[allow(unused_imports)]
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, PartitionFmt, SeqMatrix};
use crate::parser::delimited;
use crate::writer::sequences::SeqWriter;

#[allow(dead_code)]
pub struct Rename<'a> {
    input_fmt: &'a InputFmt,
    datatype: &'a DataType,
    ids: &'a Path,
}

#[allow(dead_code)]
impl<'a> Rename<'a> {
    pub fn new(input_fmt: &'a InputFmt, datatype: &'a DataType, ids: &'a Path) -> Self {
        Self {
            input_fmt,
            datatype,
            ids,
        }
    }

    pub fn dry_run(&self) {
        let names = self.get_names();
        names
            .iter()
            .for_each(|(origin, destination)| log::info!("{} -> {}", origin, destination));
    }

    #[allow(unused_variables)]
    pub fn rename(&self, files: &[PathBuf], outdir: &Path, output_fmt: &OutputFmt) {
        let names = self.get_names();
        files.iter().for_each(|file| {
            let (mut seq, header) = Sequence::new(file, self.datatype).get(self.input_fmt);
            let original_size = seq.len();
            names.iter().for_each(|(origin, destination)| {
                let values = seq.remove(origin);
                match values {
                    Some(value) => {
                        seq.insert(destination.to_string(), value);
                    }
                    None => (),
                }
            });

            assert_eq!(original_size, seq.len());
            let mut writer = SeqWriter::new(outdir, &seq, &header, None, &PartitionFmt::None);
            writer
                .write_sequence(output_fmt)
                .expect("Failed writing output sequence");
        });
    }

    fn get_names(&self) -> Vec<(String, String)> {
        delimited::parse_delimited_text(self.ids)
    }
}
