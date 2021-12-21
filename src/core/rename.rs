use std::path::{Path, PathBuf};

use ansi_term::Colour::Yellow;
use rayon::prelude::*;

use crate::helper::filenames;
use crate::helper::sequence::Sequence;
use crate::helper::types::{DataType, InputFmt, OutputFmt, PartitionFmt};
use crate::helper::utils;
use crate::parser::delimited;
use crate::writer::sequences::SeqWriter;

#[allow(dead_code)]
pub struct Rename<'a> {
    input_fmt: &'a InputFmt,
    datatype: &'a DataType,
    ids: &'a Path,
}

impl<'a> Rename<'a> {
    pub fn new(input_fmt: &'a InputFmt, datatype: &'a DataType, ids: &'a Path) -> Self {
        Self {
            input_fmt,
            datatype,
            ids,
        }
    }

    pub fn dry_run(&self) {
        log::info!("{:18}: Dry run", "Status");
        log::info!("Results:");
        let names = self.get_names();
        names
            .iter()
            .for_each(|(origin, destination)| log::info!("{} -> {}", origin, destination));
    }

    pub fn rename(&self, files: &[PathBuf], outdir: &Path, output_fmt: &OutputFmt) {
        let spin = utils::set_spinner();
        spin.set_message("Batch renaming dna sequence IDs...");
        let names = self.get_names();
        files.par_iter().for_each(|file| {
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
            let outpath = filenames::create_output_fname(outdir, file, output_fmt);
            let mut writer = SeqWriter::new(&outpath, &seq, &header, None, &PartitionFmt::None);
            writer
                .write_sequence(output_fmt)
                .expect("Failed writing output sequence");
        });
        spin.finish_with_message("Finished batch renaming dna sequence IDs!\n");
        self.print_output_info(outdir);
    }

    fn get_names(&self) -> Vec<(String, String)> {
        delimited::parse_delimited_text(self.ids)
    }

    fn print_output_info(&self, output: &Path) {
        log::info!("{}", Yellow.paint("Output"));
        log::info!("{:18}: {}", "Output dir", output.display());
    }
}
