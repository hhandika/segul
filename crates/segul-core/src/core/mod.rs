//! Handle segul main functions.
#[cfg(feature = "align")]
pub mod align;
#[cfg(feature = "contig")]
pub mod contig;
#[cfg(feature = "maf")]
pub mod maf;
#[cfg(feature = "read")]
pub mod read;
#[cfg(feature = "sequence")]
pub mod sequence;

#[cfg(any(feature = "align", feature = "sequence", feature = "maf"))]
use crate::helper::types::OutputFmt;

#[cfg(any(feature = "align", feature = "sequence", feature = "maf"))]
macro_rules! log_output_fmt {
    ($filetype: expr_2021) => {
        log::info!("{:18}: {}", "Output format", $filetype)
    };
}

#[cfg(any(feature = "align", feature = "sequence", feature = "maf"))]
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

#[cfg(feature = "align")]
trait PartitionPrint {
    fn print_partition_info(&self, part_path: &std::path::Path, part_counts: &usize) {
        use crate::helper::utils;
        use colored::Colorize;

        log::info!("{}", "Partitions".yellow());
        log::info!("{:18}: {}", "Partition counts", utils::fmt_num(part_counts));
        log::info!("{:18}: {}\n", "File path", part_path.display());
    }
}
