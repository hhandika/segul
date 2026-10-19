//! Helper functions for converting alignment partition.
use std::{
    ffi::OsStr,
    path::{Path, PathBuf},
};

use super::types::PartitionFmt;

pub fn construct_partition_path(input: &Path, out_part_fmt: &PartitionFmt) -> PathBuf {
    let file_stem = input
        .file_stem()
        .and_then(OsStr::to_str)
        .expect("Failed to parse input file stem");
    let mut fname = PathBuf::from(format!("{}_partition", file_stem));
    match *out_part_fmt {
        PartitionFmt::Nexus | PartitionFmt::NexusCodon => {
            fname.set_extension("nex");
        }
        PartitionFmt::Raxml | PartitionFmt::RaxmlCodon => {
            fname.set_extension("txt");
        }
        _ => unreachable!("Failed to parse partition format"),
    }
    let parent_path = input.parent().expect("Failed to parse input parent path");
    parent_path.join(fname)
}
