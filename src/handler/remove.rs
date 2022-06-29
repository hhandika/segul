use std::path::{Path, PathBuf};

use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix};
use crate::helper::utils;

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

    pub fn remove(&self, file: &[PathBuf]) {
        let spin = utils::set_spinner();
        spin.set_message("Removing sequences...");
        match self.opts {
            RemoveOpts::Id(names) => {}
            RemoveOpts::Regex(input_re) => {}
        }
        spin.finish_with_message("Finished removing sequences!\n");
    }
}
