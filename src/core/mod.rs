//! Handle segul main functions.
pub mod align;
pub mod contig;
pub mod read;
pub mod sequence;

use std::path::Path;

use colored::Colorize;

use crate::helper::types::OutputFmt;
use crate::helper::utils;

macro_rules! log_output_fmt {
    ($filetype: expr) => {
        log::info!("{:18}: {}", "Output format", $filetype)
    };
}

trait OutputPrint {
    fn print_output_fmt(&self, output_fmt: &OutputFmt) {
        match output_fmt {
            OutputFmt::Fasta => log_output_fmt!("FASTA sequential"),
            OutputFmt::Nexus => log_output_fmt!("NEXUS sequential"),
            OutputFmt::Phylip => log_output_fmt!("PHYLIP sequential"),
            OutputFmt::FastaInt => log_output_fmt!("FASTA interleaved"),
            OutputFmt::NexusInt => log_output_fmt!("NEXUS interleaved"),
            OutputFmt::PhylipInt => log_output_fmt!("PHYLIP interleaved"),
        }
    }
}

trait PartitionPrint {
    fn print_partition_info(&self, part_path: &Path, part_counts: &usize) {
        log::info!("{}", "Partitions".yellow());
        log::info!("{:18}: {}", "Partition counts", utils::fmt_num(part_counts));
        log::info!("{:18}: {}\n", "File path", part_path.display());
    }
}
