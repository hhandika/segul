//! Write a list of string to a file

use std::{
    io::{Result, Write},
    path::{Path, PathBuf},
};

use indexmap::IndexSet;

use crate::core::sequence::id::IdRecords;

use super::FileWriter;

const ID_EXTENSION: &str = "txt";
const MAP_EXTENSION: &str = "csv";
const ID_SUFFIX: &str = "id";
const MAP_SUFFIX: &str = "map";

pub struct IdWriter<'a> {
    output: &'a Path,
    id_list: &'a IndexSet<String>,
    prefix: Option<&'a str>,
}

impl FileWriter for IdWriter<'_> {}

impl<'a> IdWriter<'a> {
    pub fn new(output: &'a Path, id_list: &'a IndexSet<String>, prefix: Option<&'a str>) -> Self {
        Self {
            output,
            id_list,
            prefix,
        }
    }

    pub fn write_unique_id(&self) -> Result<()> {
        let output_path = self.create_final_output_path(ID_SUFFIX, ID_EXTENSION);
        let mut writer = self
            .create_output_file(&output_path)
            .expect("Failed creating output file");
        self.id_list.iter().for_each(|id| {
            writeln!(writer, "{}", id).unwrap();
        });
        writer.flush()?;
        Ok(())
    }

    pub fn write_mapped_id(&self, mapped_ids: &[IdRecords]) -> Result<()> {
        let output_path = self.create_final_output_path(MAP_SUFFIX, MAP_EXTENSION);
        let mut writer = self
            .create_output_file(&output_path)
            .expect("Failed creating output file");
        write!(writer, "locus")?;
        self.id_list.iter().for_each(|id| {
            write!(writer, ",{}", id).expect("Failed writing a csv header");
        });
        writeln!(writer)?;
        mapped_ids.iter().for_each(|rec| {
            write!(writer, "{}", rec.name).expect("Failed writing a csv header");
            rec.records.iter().for_each(|is_id| {
                write!(writer, ",{}", is_id).expect("Failed writing id map");
            });
            writeln!(writer).expect("Failed writing id map");
        });
        writer.flush()?;
        Ok(())
    }

    fn create_final_output_path(&self, suffix: &str, extension: &str) -> PathBuf {
        match self.prefix {
            Some(prefix) => {
                let file_name = format!("{}_{}", prefix, suffix);
                self.output.join(file_name).with_extension(extension)
            }
            None => self.output.join(suffix).with_extension(extension),
        }
    }
}
