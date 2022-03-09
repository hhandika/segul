pub mod concat;
pub mod convert;
pub mod extract;
pub mod filter;
pub mod id;
pub mod partition;
pub mod rename;
pub mod split;
pub mod summarize;
pub mod translate;

use std::path::Path;

use ansi_term::Colour::Yellow;

use crate::helper::types::OutputFmt;
use crate::helper::utils;

macro_rules! log_output_fmt {
    ($ftype: expr) => {
        log::info!("{:18}: {}", "Output format", $ftype)
    };
}

trait OutputPrint {
    fn print_output_fmt(&self, output_fmt: &OutputFmt) {
        match output_fmt {
            OutputFmt::Fasta => log_output_fmt!("Fasta Sequential"),
            OutputFmt::Nexus => log_output_fmt!("Nexus Sequential"),
            OutputFmt::Phylip => log_output_fmt!("Phylip Sequential"),
            OutputFmt::FastaInt => log_output_fmt!("Fasta Interleaved"),
            OutputFmt::NexusInt => log_output_fmt!("Nexus Interleaved"),
            OutputFmt::PhylipInt => log_output_fmt!("Phylip Interleaved"),
        }
    }
}

trait PartitionPrint {
    fn print_partition_info(&self, part_path: &Path, part_counts: &usize) {
        log::info!("{}", Yellow.paint("Partitions"));
        log::info!("{:18}: {}", "Partition counts", utils::fmt_num(part_counts));
        log::info!("{:18}: {}\n", "File path", part_path.display());
    }
}
