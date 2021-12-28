pub mod concat;
pub mod convert;
pub mod extract;
pub mod filter;
pub mod id;
pub mod rename;
pub mod summarize;
pub mod translate;

use crate::helper::types::OutputFmt;

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
