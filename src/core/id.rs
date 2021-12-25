use std::ffi::OsStr;
use std::fs::{File, OpenOptions};
use std::io::prelude::*;
use std::io::{BufWriter, Result};
use std::path::{Path, PathBuf};
// use std::sync::{Arc, Mutex};

use ansi_term::Colour::Yellow;
use indexmap::IndexSet;
// use rayon::prelude::*;

use crate::helper::finder::IDs;
use crate::helper::sequence::Sequence;
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
        let ids = self.get_unique_id(files);
        spin.finish_with_message("DONE!\n");
        self.write_unique_id(&ids).expect("Failed writing results");
        self.print_output(ids.len());
    }

    pub fn map_id(&self, files: &[PathBuf], output_id: &Path) {
        let spin = utils::set_spinner();
        spin.set_message("Mapping IDs..");
        let ids = self.get_unique_id(files);
        let mapped_ids = self.map_id_to_aln(files, &ids);
        self.write_unique_id(&ids)
            .expect("Failed writing unique IDs to file");
        self.write_mapped_id(&ids, &mapped_ids, output_id)
            .expect("Failed writing mapped ID to file");
        spin.finish_with_message("DONE!\n");
        self.print_output(ids.len());
        log::info!("{:18}: {}", "Mapped ID output", output_id.display());
    }

    fn get_unique_id(&self, files: &[PathBuf]) -> IndexSet<String> {
        let mut id = IDs::new(files, self.input_fmt, self.datatype).get_id_unique();
        id.sort();
        id
    }

    // fn map_id_to_aln(
    //     &self,
    //     files: &[PathBuf],
    //     ids: &IndexSet<String>,
    // ) -> IndexMap<String, Vec<bool>> {
    //     let mut mapped_ids: IndexMap<String, Vec<bool>> = IndexMap::new();
    //     files.iter().for_each(|file| {
    //         let (seq, _) = Sequence::new(file, self.datatype).get(self.input_fmt);
    //         let mut is_id = Vec::with_capacity(ids.len());
    //         ids.iter().for_each(|id| {
    //             is_id.push(seq.contains_key(id));
    //         });
    //         let fstem = file
    //             .file_stem()
    //             .and_then(OsStr::to_str)
    //             .expect("Failed getting file stem for mapping IDs")
    //             .to_string();
    //         mapped_ids.insert(fstem, is_id);
    //     });
    //     mapped_ids
    // }

    fn map_id_to_aln(&self, files: &[PathBuf], ids: &IndexSet<String>) -> Vec<IdRecords> {
        let mut mapped_ids = Vec::new();
        files.iter().for_each(|file| {
            let fstem = self.get_aln_name(file);
            let mut rec = IdRecords::new(fstem, ids.len());
            let (seq, _) = Sequence::new(file, self.datatype).get(self.input_fmt);
            ids.iter().for_each(|id| {
                rec.records.push(seq.contains_key(id));
            });
            mapped_ids.push(rec);
        });
        mapped_ids
    }

    fn get_aln_name(&self, file: &Path) -> String {
        file.file_stem()
            .and_then(OsStr::to_str)
            .expect("Failed getting file stem for mapping IDs")
            .to_string()
    }

    fn write_unique_id(&self, ids: &IndexSet<String>) -> Result<()> {
        let mut writer = self.write_file(self.output);
        ids.iter().for_each(|id| {
            writeln!(writer, "{}", id).unwrap();
        });
        writer.flush()?;
        Ok(())
    }

    fn write_mapped_id(
        &self,
        ids: &IndexSet<String>,
        mapped_ids: &[IdRecords],
        output: &Path,
    ) -> Result<()> {
        let mut writer = self.write_file(output);
        write!(writer, "Alignments")?;
        ids.iter().for_each(|id| {
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

    fn write_file(&self, output: &Path) -> BufWriter<File> {
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output)
            .expect("Failed writing id results");
        BufWriter::new(file)
    }

    fn print_output(&self, ids: usize) {
        log::info!("{}", Yellow.paint("Output"));
        log::info!("{:18}: {}", "Total unique IDs", ids);
        log::info!("{:18}: {}", "File output", self.output.display());
    }
}

struct IdRecords {
    name: String,
    records: Vec<bool>,
}

impl IdRecords {
    fn new(name: String, size: usize) -> Self {
        Self {
            name,
            records: Vec::with_capacity(size),
        }
    }
}
