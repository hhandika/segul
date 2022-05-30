use std::path::{Path, PathBuf};

use ansi_term::Colour::Yellow;
use rayon::prelude::*;

use crate::handler::OutputPrint;
use crate::helper::filenames;
use crate::helper::sequence::Sequence;
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix};
use crate::helper::utils;
use crate::parser::delimited;
use crate::writer::sequences::SeqWriter;

impl OutputPrint for Rename<'_> {}

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
        let names = self.get_names();
        log::info!("{:18}: {}", "New ID count", names.len());
        log::info!("{:18}: Dry run\n", "Status");
        log::info!("{}", Yellow.paint("Results"));
        names
            .iter()
            .for_each(|(origin, destination)| log::info!("{} -> {}", origin, destination));
        println!();
    }

    pub fn rename(&self, files: &[PathBuf], outdir: &Path, output_fmt: &OutputFmt) {
        let names = self.get_names();
        log::info!("{:18}: {}\n", "New ID count", names.len());
        let spin = utils::set_spinner();
        spin.set_message("Batch renaming dna sequence IDs...");
        files.par_iter().for_each(|file| {
            let (seq, header) = self.rename_seq_id(file, &names);
            let outpath = filenames::create_output_fname(outdir, file, output_fmt);
            let mut writer = SeqWriter::new(&outpath, &seq, &header);
            writer
                .write_sequence(output_fmt)
                .expect("Failed writing output sequence");
        });
        spin.finish_with_message("Finished batch renaming dna sequence IDs!\n");
        self.print_output_info(outdir, output_fmt);
    }

    fn rename_seq_id(&self, file: &Path, names: &[(String, String)]) -> (SeqMatrix, Header) {
        let (mut seq, header) = Sequence::new(file, self.datatype).get(self.input_fmt);
        let original_size = seq.len();
        names.iter().for_each(|(origin, destination)| {
            let values = seq.remove(origin);
            if let Some(value) = values {
                seq.insert(destination.to_string(), value);
            }
        });

        assert_eq!(
            original_size,
            seq.len(),
            "Failed renaming files. New ID counts does not match original ID counts. \
         Original ID counts: {}. New ID counts: {}",
            original_size,
            seq.len()
        );
        (seq, header)
    }

    fn get_names(&self) -> Vec<(String, String)> {
        delimited::parse_delimited_text(self.ids)
    }

    fn print_output_info(&self, output: &Path, output_fmt: &OutputFmt) {
        log::info!("{}", Yellow.paint("Output"));
        log::info!("{:18}: {}", "Output dir", output.display());
        self.print_output_fmt(output_fmt);
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_rename_seq() {
        let input_fmt = InputFmt::Fasta;
        let datatype = DataType::Dna;
        let file = Path::new("tests/files/simple.fas");
        let names = [(String::from("ABCD"), String::from("WXYZ"))];
        let ids = Path::new("tests/delimited.tsv");
        let rename = Rename::new(&input_fmt, &datatype, &ids);
        let (seq, _) = rename.rename_seq_id(&file, &names);
        assert_eq!(seq.len(), 2);
        assert_eq!(seq.get("WXYZ"), Some(&String::from("AGTATG")));
        assert_eq!(seq.get("ABCD"), None);
    }
}
