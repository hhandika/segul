//! Handle segul main functions. Mostly used by the Cli app.
pub mod concat;
pub mod convert;
pub mod extract;
pub mod filter;
pub mod id;
pub mod partition;
pub mod remove;
pub mod rename;
pub mod split;
pub mod summarize;
pub mod translate;

use std::path::Path;

use colored::Colorize;

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
