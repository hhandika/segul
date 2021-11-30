use std::path::{Path, PathBuf};

use crate::helper::types::OutputFmt;

pub fn create_output_fname(path: &Path, output_fmt: &OutputFmt) -> PathBuf {
    match output_fmt {
        OutputFmt::Fasta | OutputFmt::FastaInt => path.with_extension("fas"),
        OutputFmt::Nexus | OutputFmt::NexusInt => path.with_extension("nex"),
        OutputFmt::Phylip | OutputFmt::PhylipInt => path.with_extension("phy"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_output_fname() {
        let path = Path::new("tests/test_create_output_fname.nex");

        assert_eq!(
            create_output_fname(path, &OutputFmt::Fasta),
            Path::new("tests/test_create_output_fname.fas")
        );
    }
}
