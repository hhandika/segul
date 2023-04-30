//! Helper functions for creating output filenames
use std::path::{Path, PathBuf};
use std::{fs::File, io::BufReader};

use crate::helper::types::OutputFmt;

use flate2::read::MultiGzDecoder;

/// Decode gzip compressed files
/// Returns a BufReader of the decoded file
pub fn decode_gzip(path: &Path) -> BufReader<MultiGzDecoder<File>> {
    let file = File::open(path).expect("Failed to open file");
    let decoder = MultiGzDecoder::new(file);
    BufReader::new(decoder)
}

/// Open a file
/// Returns a BufReader of the file
pub fn open_file(path: &Path) -> BufReader<File> {
    let file = File::open(path).expect("Failed opening fastq file");
    BufReader::new(file)
}

/// Combine the output directory and the input filename
/// Returns a PathBuf of the output path
/// # Arguments
/// * `dir` - Output directory
/// * `file` - Input filename
/// * `output_fmt` - Output format
/// # Example
/// ```
/// use std::path::Path;
/// use segul::helper::filenames;
/// use segul::helper::types::OutputFmt;
///
/// let dir = Path::new("output");
/// let file = Path::new("input.fas");
/// let output_fmt = OutputFmt::Fasta;
/// let output = filenames::create_output_fname(&dir, &file, &output_fmt);
/// assert_eq!(output, Path::new("output/input.fas"));
/// ```
pub fn create_output_fname(dir: &Path, file: &Path, output_fmt: &OutputFmt) -> PathBuf {
    let path = dir.join(
        file.file_name()
            .expect("Failed parsing filename for output file"),
    );
    create_output_fname_from_path(&path, output_fmt)
}

/// Create output filename from input filename
/// Returns a PathBuf of the output filename
/// # Arguments
/// * `path` - Input filename
/// * `output_fmt` - Output format
/// # Example
/// ```
/// use std::path::Path;
/// use segul::helper::filenames;
/// use segul::helper::types::OutputFmt;
///
/// let path = Path::new("input.fas");
/// let output_fmt = OutputFmt::Fasta;
/// let output = filenames::create_output_fname_from_path(&path, &output_fmt);
/// assert_eq!(output, Path::new("input.fas"));
/// ```
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
