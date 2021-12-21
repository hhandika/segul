use std::path::{Path, PathBuf};

use crate::helper::types::OutputFmt;

pub fn create_output_fname(dir: &Path, file: &Path, output_fmt: &OutputFmt) -> PathBuf {
    let path = dir.join(
        file.file_name()
            .expect("Failed parsing filename for output file"),
    );
    create_output_fname_from_path(&path, output_fmt)
}

pub fn create_output_fname_from_path(path: &Path, output_fmt: &OutputFmt) -> PathBuf {
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
        let dir = Path::new("tests");
        assert_eq!(
            create_output_fname(dir, path, &OutputFmt::Fasta),
            Path::new("tests/test_create_output_fname.fas")
        );
    }
}
