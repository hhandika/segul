#[allow(dead_code, unused_imports)]
use std::ffi::OsStr;
use std::path::{Path, PathBuf};

#[allow(dead_code, unused_imports)]
use crate::helper::finder::IDs;
use crate::helper::sequence::Sequence;
use crate::helper::types::{
    DataType, Header, InputFmt, OutputFmt, Partition, PartitionFmt, SeqMatrix,
};
use crate::helper::utils;
use crate::writer::sequences::SeqWriter;

#[allow(dead_code)]
pub struct Splitter<'a> {
    input_fmt: &'a InputFmt,
    output: &'a Path,
    output_fmt: &'a OutputFmt,
    partition: &'a Path,
}

impl<'a> Splitter<'a> {
    pub fn new(
        input_fmt: &'a InputFmt,
        output: &'a Path,
        output_fmt: &'a OutputFmt,
        partition: &'a Path,
    ) -> Self {
        Self {
            input_fmt,
            output,
            output_fmt,
            partition,
        }
    }

    #[allow(dead_code)]
    fn parse_partition_raxml(&self) {
        let mut partition = Partition::new();
    }
}
