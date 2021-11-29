use std::path::{Path, PathBuf};

use crate::helper::types::OutputFmt;

pub fn create_output_fname(path: &Path, output_fmt: &OutputFmt) -> PathBuf {
    match output_fmt {
        OutputFmt::Fasta | OutputFmt::FastaInt => path.with_extension("fas"),
        OutputFmt::Nexus | OutputFmt::NexusInt => path.with_extension("nex"),
        OutputFmt::Phylip | OutputFmt::PhylipInt => path.with_extension("phy"),
    }
}
