use std::fs;
use std::path::{Path, PathBuf};

use crate::helper::types::{
    DataType, GeneticCodes, Header, InputFmt, OutputFmt, PartitionFmt, SeqMatrix,
};
use crate::parser::delimited;

#[allow(dead_code)]
pub struct Rename<'a> {
    input_fmt: &'a InputFmt,
    datatype: &'a DataType,
}

#[allow(dead_code)]
impl<'a> Rename<'a> {
    pub fn new(input_fmt: &'a InputFmt, datatype: &'a DataType) -> Self {
        Self {
            input_fmt,
            datatype,
        }
    }

    #[allow(unused_variables)]
    pub fn rename<P: AsRef<Path>>(&self, ids: &P, outdir: &Path, output_fmt: &OutputFmt) {
        let names = delimited::parse_delimited_text(ids.as_ref());
        names
            .iter()
            .for_each(|(origin, destination)| log::info!("{} -> {}", origin, destination));
    }
}
