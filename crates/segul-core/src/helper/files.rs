//! Helper functions for creating output filenames
use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::{fs::File, io::BufReader};

use anyhow::Context;

use crate::helper::types::OutputFmt;

use flate2::read::MultiGzDecoder;

/// Check whether a file is gzip compressed by its extension (`.gz` or `.gzip`).
/// # Example
/// ```
/// use std::path::Path;
/// use segul_core::helper::files;
///
/// assert!(files::is_gzip(Path::new("genome.maf.gz")));
/// assert!(!files::is_gzip(Path::new("genome.maf")));
/// ```
pub fn is_gzip(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.eq_ignore_ascii_case("gz") || ext.eq_ignore_ascii_case("gzip"))
        .unwrap_or(false)
}

/// Open a plain text or gzip compressed file for buffered reading.
/// Gzip files are detected by their extension.
pub fn open_reader(path: &Path) -> anyhow::Result<Box<dyn BufRead + Send>> {
    let file = File::open(path).with_context(|| format!("Failed opening {}", path.display()))?;
    if is_gzip(path) {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

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
    let file = File::open(path).expect("Failed opening a file");
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
/// use segul_core::helper::files;
/// use segul_core::helper::types::OutputFmt;
///
/// let dir = Path::new("output");
/// let file = Path::new("input.fas");
/// let output_fmt = OutputFmt::Fasta;
/// let output = files::create_output_fname(&dir, &file, &output_fmt);
/// assert_eq!(output, Path::new("output/input.fas"));
/// ```
pub fn create_output_fname(dir: &Path, file: &Path, output_fmt: &OutputFmt) -> PathBuf {
    let path = dir.join(
        file.file_name()
            .expect("Failed parsing filename for output file"),
    );
    create_output_fname_from_path(&path, output_fmt)
}

/// Combine the output directory and a file stem.
/// Unlike [`create_output_fname`], dots in the stem are kept.
/// # Example
/// ```
/// use std::path::Path;
/// use segul_core::helper::files;
/// use segul_core::helper::types::OutputFmt;
///
/// let output = files::create_output_fname_from_stem(
///     Path::new("output"),
///     "hg18.chr7-10-20",
///     &OutputFmt::Fasta,
/// );
/// assert_eq!(output, Path::new("output").join("hg18.chr7-10-20.fas"));
/// ```
pub fn create_output_fname_from_stem(dir: &Path, stem: &str, output_fmt: &OutputFmt) -> PathBuf {
    dir.join(format!("{}.{}", stem, output_extension(output_fmt)))
}

/// File extension for an output format, without the dot.
pub fn output_extension(output_fmt: &OutputFmt) -> &'static str {
    match output_fmt {
        OutputFmt::Fasta | OutputFmt::FastaInt => "fas",
        OutputFmt::Nexus | OutputFmt::NexusInt => "nex",
        OutputFmt::Phylip | OutputFmt::PhylipInt => "phy",
    }
}

pub fn create_output_fname_for_text(dir: &Path, prefix: &Path) -> PathBuf {
    dir.join(prefix.with_extension("txt"))
}

/// Create output filename from input filename
/// Returns a PathBuf of the output filename
/// # Arguments
/// * `path` - Input filename
/// * `output_fmt` - Output format
/// # Example
/// ```
/// use std::path::Path;
/// use segul_core::helper::files;
/// use segul_core::helper::types::OutputFmt;
///
/// let path = Path::new("input.fas");
/// let output_fmt = OutputFmt::Fasta;
/// let output = files::create_output_fname_from_path(&path, &output_fmt);
/// assert_eq!(output, Path::new("input.fas"));
/// ```
pub fn create_output_fname_from_path(path: &Path, output_fmt: &OutputFmt) -> PathBuf {
    path.with_extension(output_extension(output_fmt))
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

    #[test]
    fn test_is_gzip() {
        assert!(is_gzip(Path::new("a.maf.gz")));
        assert!(is_gzip(Path::new("a.BED.GZ")));
        assert!(is_gzip(Path::new("a.gzip")));
        assert!(!is_gzip(Path::new("a.maf")));
        assert!(!is_gzip(Path::new("gz")));
    }

    #[test]
    fn test_output_fname_from_stem_keeps_dots() {
        let output = create_output_fname_from_stem(Path::new("out"), "a.b", &OutputFmt::Nexus);
        assert_eq!(output, Path::new("out").join("a.b.nex"));
    }
}
