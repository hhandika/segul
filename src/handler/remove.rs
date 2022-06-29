use std::path::{Path, PathBuf};

use ansi_term::Colour::Yellow;

use crate::handler::OutputPrint;
use crate::helper::filenames;
use crate::helper::sequence::{SeqCheck, SeqParser};
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix};
use crate::helper::utils;
use crate::writer::sequences::SeqWriter;

impl OutputPrint for Remove<'_> {}

pub enum RemoveOpts {
    Id(Vec<String>),
    Regex(String),
}

pub struct Remove<'a> {
    input_fmt: &'a InputFmt,
    datatype: &'a DataType,
    outdir: &'a Path,
    output_fmt: &'a OutputFmt,
    opts: &'a RemoveOpts,
}

impl<'a> Remove<'a> {
    pub fn new(
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        outdir: &'a Path,
        output_fmt: &'a OutputFmt,
        opts: &'a RemoveOpts,
    ) -> Self {
        Self {
            input_fmt,
            datatype,
            outdir,
            output_fmt,
            opts,
        }
    }

    pub fn remove(&self, files: &[PathBuf]) {
        let spin = utils::set_spinner();
        spin.set_message("Removing sequences...");
        match self.opts {
            RemoveOpts::Id(ids) => files.iter().for_each(|file| {
                let (matrix, header) = self.remove_sequence(file, ids);
                self.write_output(&matrix, &header, file);
            }),
            RemoveOpts::Regex(input_re) => println!("Regex {} is working!", input_re),
        }
        spin.finish_with_message("Finished removing sequences!\n");
        self.print_output_info();
    }

    fn write_output(&self, matrix: &SeqMatrix, header: &Header, file: &Path) {
        let outpath = filenames::create_output_fname(self.outdir, file, self.output_fmt);
        let mut writer = SeqWriter::new(&outpath, matrix, header);
        writer
            .write_sequence(self.output_fmt)
            .expect("Failed writing output sequence");
    }

    fn remove_sequence(&self, fpath: &Path, ids: &[String]) -> (SeqMatrix, Header) {
        let (mut matrix, _) = SeqParser::new(fpath, self.datatype).get(self.input_fmt);
        ids.iter().for_each(|id| match matrix.remove(id) {
            Some(_) => (),
            None => (),
        });

        let header = self.get_header(&matrix);
        (matrix, header)
    }

    fn get_header(&self, matrix: &SeqMatrix) -> Header {
        let mut seq_info = SeqCheck::new();
        seq_info.check(matrix);
        let mut header = Header::new();
        header.aligned = seq_info.is_alignment;
        header.nchar = seq_info.longest;
        header.ntax = matrix.len();
        header
    }

    fn print_output_info(&self) {
        log::info!("{}", Yellow.paint("Output"));
        log::info!("{:18}: {}", "Output dir", self.outdir.display());
        self.print_output_fmt(self.output_fmt);
    }
}
