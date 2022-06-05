use std::path::{Path, PathBuf};

use ansi_term::Colour::Yellow;
use rayon::prelude::*;

use crate::handler::OutputPrint;
use crate::helper::filenames;
use crate::helper::sequence::Sequence;
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix};
use crate::helper::utils;
use crate::writer::sequences::SeqWriter;

impl OutputPrint for Rename<'_> {}

#[allow(dead_code)]
pub enum RenameOpts {
    RnId(Vec<(String, String)>), // Rename ID using tabulated file
    RmStr(String),               // Remove characters in seq id using string input
    RmRegex(String),             // Similar to RmStr but using regex as input
}

pub struct Rename<'a> {
    input_fmt: &'a InputFmt,
    datatype: &'a DataType,
    outdir: &'a Path,
    output_fmt: &'a OutputFmt,
}

impl<'a> Rename<'a> {
    pub fn new(
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        outdir: &'a Path,
        output_fmt: &'a OutputFmt,
    ) -> Self {
        Self {
            input_fmt,
            datatype,
            outdir,
            output_fmt,
        }
    }

    pub fn dry_run(&self, opts: &RenameOpts) {
        let names = match opts {
            RenameOpts::RnId(names) => names,
            _ => unimplemented!(),
        };
        log::info!("{:18}: {}", "New ID count", names.len());
        log::info!("{:18}: Dry run\n", "Status");
        log::info!("{}", Yellow.paint("Results"));
        names
            .iter()
            .for_each(|(origin, destination)| log::info!("{} -> {}", origin, destination));
        println!();
    }

    pub fn rename(&self, files: &[PathBuf], opts: &RenameOpts) {
        let spin = utils::set_spinner();
        spin.set_message("Batch renaming dna sequence IDs...");
        files.par_iter().for_each(|file| {
            let (seq, header) = self.parse_opts(file, opts);
            let outpath = filenames::create_output_fname(self.outdir, file, self.output_fmt);
            let mut writer = SeqWriter::new(&outpath, &seq, &header);
            writer
                .write_sequence(self.output_fmt)
                .expect("Failed writing output sequence");
        });
        spin.finish_with_message("Finished batch renaming dna sequence IDs!\n");
        self.print_output_info(self.outdir, self.output_fmt);
    }

    pub fn parse_opts(&self, file: &Path, opts: &RenameOpts) -> (SeqMatrix, Header) {
        match opts {
            RenameOpts::RnId(names) => self.rename_seq_id(file, names),
            RenameOpts::RmStr(_) => unimplemented!(),
            RenameOpts::RmRegex(_) => unimplemented!(),
        }
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
        let rename = Rename::new(&input_fmt, &datatype, Path::new("."), &OutputFmt::Nexus);
        let (seq, _) = rename.rename_seq_id(&file, &names);
        assert_eq!(seq.len(), 2);
        assert_eq!(seq.get("WXYZ"), Some(&String::from("AGTATG")));
        assert_eq!(seq.get("ABCD"), None);
    }
}
