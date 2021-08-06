use std::fs::OpenOptions;
use std::io::prelude::*;
use std::io::{BufWriter, Result};
use std::path::{Path, PathBuf};

use ansi_term::Colour::Yellow;
use indexmap::IndexSet;

use crate::helper::finder::IDs;
use crate::helper::types::{DataType, InputFmt};
use crate::helper::utils;

pub struct Id<'a> {
    pub output: &'a Path,
    pub input_fmt: &'a InputFmt,
    pub datatype: &'a DataType,
}

impl<'a> Id<'a> {
    pub fn new(output: &'a Path, input_fmt: &'a InputFmt, datatype: &'a DataType) -> Self {
        Self {
            output,
            input_fmt,
            datatype,
        }
    }

    pub fn generate_id(&self, files: &[PathBuf]) {
        let spin = utils::set_spinner();
        spin.set_message("Indexing IDs..");
        let ids = IDs::new(files, self.input_fmt, self.datatype).get_id_unique();
        spin.finish_with_message("DONE!\n");
        self.write_results(&ids).expect("Failed writing results");
        self.print_output(ids.len());
    }

    #[allow(dead_code)]
    pub fn map_id(&self, files: &[PathBuf]) {
        let spin = utils::set_spinner();
        spin.set_message("Indexing IDs..");
        let ids = IDs::new(files, self.input_fmt, self.datatype);
        let taxon_id = ids.get_id_unique();
        println!("{}", taxon_id.len());
    }

    fn write_results(&self, ids: &IndexSet<String>) -> Result<()> {
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&self.output)
            .expect("Failed writing id results");
        let mut writer = BufWriter::new(file);
        ids.iter().for_each(|id| {
            writeln!(writer, "{}", id).unwrap();
        });
        writer.flush()?;
        Ok(())
    }

    fn print_output(&self, ids: usize) {
        log::info!("\n{}", Yellow.paint("Output"));
        log::info!("{:18}: {}", "Total unique IDs", ids);
        log::info!("{:18}: {}", "File output", self.output.display());
    }
}

// #[cfg(test)]
// mod test {
//     use super::*;
//     use crate::helper::finder::Files;
