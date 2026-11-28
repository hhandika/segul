//! High-performance phylogenetic tree parsers for Newick, Extended Newick, NHX, NEXUS, PhyloXML, and NeXML.

pub mod comments;
pub mod model;
pub mod newick;
pub mod nexml;
pub mod nexus;
pub mod phyloxml;
pub(crate) mod xml_pull;

use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::str::FromStr;

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use walkdir::WalkDir;

use crate::helper::files;
use crate::parser::IO_BUF_SIZE;

pub use model::{
    Annotation, AnnotationFormat, AnnotationTarget, AnnotationValue, HybridInfo, Node, NodeId,
    Rootedness, TaxonId, TaxonTable, Tree, TreeMetadata,
};
pub use newick::{NewickParser, NewickReader, ParseOptions, parse_newick_bytes_par};
pub use nexml::NeXmlReader;
pub use nexus::{NexusTreeReader, TranslateTable, parse_nexus_bytes_par};
pub use phyloxml::PhyloXmlReader;

/// Input tree format selector.
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum TreeFmt {
    /// Detect format from file content (after transparent gzip decoding).
    Auto,
    /// Standard Newick and NHX-annotated Newick.
    Newick,
    /// Extended Newick (eNewick) phylogenetic networks.
    ExtendedNewick,
    /// NEXUS `TREES` block (skips `DATA` and other non-tree blocks).
    Nexus,
    /// PhyloXML XML tree format.
    PhyloXml,
    /// NeXML XML tree format.
    NeXml,
}

impl fmt::Display for TreeFmt {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Auto => write!(f, "auto"),
            Self::Newick => write!(f, "newick"),
            Self::ExtendedNewick => write!(f, "enewick"),
            Self::Nexus => write!(f, "nexus"),
            Self::PhyloXml => write!(f, "phyloxml"),
            Self::NeXml => write!(f, "nexml"),
        }
    }
}

impl FromStr for TreeFmt {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "auto" => Ok(Self::Auto),
            "newick" | "nwk" => Ok(Self::Newick),
            "enewick" | "extended-newick" | "enwk" => Ok(Self::ExtendedNewick),
            "nexus" | "nex" => Ok(Self::Nexus),
            "phyloxml" | "pxml" => Ok(Self::PhyloXml),
            "nexml" => Ok(Self::NeXml),
            _ => Err(format!("{} is not a valid tree input format", s)),
        }
    }
}

/// Output tree format selector.
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum TreeOutputFmt {
    /// Standard Newick (plus annotations when enabled).
    Newick,
    /// Extended Newick (eNewick) with reticulation tags and `:len:support:gamma` fields.
    ExtendedNewick,
    /// NEXUS `TREES` block.
    Nexus,
    /// PhyloXML XML format.
    PhyloXml,
    /// NeXML XML format.
    NeXml,
}

impl TreeOutputFmt {
    /// Default file extension for this output format, without the leading dot.
    pub fn extension(self) -> &'static str {
        match self {
            Self::Newick => "nwk",
            Self::ExtendedNewick => "enwk",
            Self::Nexus => "nex",
            Self::PhyloXml => "phyloxml",
            Self::NeXml => "nexml",
        }
    }
}

impl fmt::Display for TreeOutputFmt {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Newick => write!(f, "newick"),
            Self::ExtendedNewick => write!(f, "enewick"),
            Self::Nexus => write!(f, "nexus"),
            Self::PhyloXml => write!(f, "phyloxml"),
            Self::NeXml => write!(f, "nexml"),
        }
    }
}

impl FromStr for TreeOutputFmt {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "newick" | "nwk" => Ok(Self::Newick),
            "enewick" | "extended-newick" | "enwk" => Ok(Self::ExtendedNewick),
            "nexus" | "nex" => Ok(Self::Nexus),
            "phyloxml" | "pxml" => Ok(Self::PhyloXml),
            "nexml" => Ok(Self::NeXml),
            _ => Err(format!("{} is not a valid tree output format", s)),
        }
    }
}

/// Sniff tree format from the beginning of a decoded byte buffer.
///
/// Skips UTF-8 BOM, leading ASCII whitespace, and leading `[...]` comments
/// before checking for `#NEXUS` or XML root elements (`phyloxml`/`nexml`).
pub fn detect_tree_format(buf: &[u8]) -> TreeFmt {
    let mut pos = if buf.starts_with(&[0xEF, 0xBB, 0xBF]) {
        3
    } else {
        0
    };

    loop {
        while pos < buf.len() && buf[pos].is_ascii_whitespace() {
            pos += 1;
        }
        if pos < buf.len()
            && buf[pos] == b'['
            && let Ok((_, next)) = comments::scan_bracket_comment(buf, pos)
        {
            pos = next;
            continue;
        }
        break;
    }

    let rest = &buf[pos..];
    if rest.len() >= 6 && rest[..6].eq_ignore_ascii_case(b"#nexus") {
        TreeFmt::Nexus
    } else if rest.starts_with(b"<") {
        detect_xml_tree_format(rest).unwrap_or(TreeFmt::PhyloXml)
    } else {
        TreeFmt::Newick
    }
}

fn detect_xml_tree_format(buf: &[u8]) -> Option<TreeFmt> {
    let mut reader = xml_pull::XmlPullReader::new(std::io::Cursor::new(buf));
    loop {
        match reader.next_event().ok()? {
            xml_pull::XmlEvent::Start { name, .. } => {
                if name.eq_ignore_ascii_case("nexml") {
                    return Some(TreeFmt::NeXml);
                }
                if name.eq_ignore_ascii_case("phyloxml") || name.eq_ignore_ascii_case("phylogeny") {
                    return Some(TreeFmt::PhyloXml);
                }
                return None;
            }
            xml_pull::XmlEvent::Text(_) => {}
            xml_pull::XmlEvent::End(_) | xml_pull::XmlEvent::Eof => return None,
        }
    }
}

enum ReaderKind<R: BufRead> {
    Newick(NewickReader<R>),
    Nexus(NexusTreeReader<R>),
    PhyloXml(PhyloXmlReader<R>),
    NeXml(NeXmlReader<R>),
}

/// Unified streaming tree reader supporting Newick, Extended Newick, NHX, NEXUS, PhyloXML, NeXML, and gzip.
pub struct TreeReader<R: BufRead> {
    inner: ReaderKind<R>,
}

impl<R: BufRead> TreeReader<R> {
    /// Construct a streaming reader from a buffered source.
    ///
    /// When `fmt` is [`TreeFmt::Auto`], peeks the buffered header bytes without
    /// consuming them to select the format.
    pub fn from_reader(mut reader: R, fmt: TreeFmt, opts: ParseOptions) -> Result<Self> {
        let resolved = match fmt {
            TreeFmt::Auto => {
                let peek = reader.fill_buf()?;
                detect_tree_format(peek)
            }
            other => other,
        };
        let inner = match resolved {
            TreeFmt::Nexus => ReaderKind::Nexus(NexusTreeReader::new(reader, opts)),
            TreeFmt::PhyloXml => ReaderKind::PhyloXml(PhyloXmlReader::new(reader, opts)),
            TreeFmt::NeXml => ReaderKind::NeXml(NeXmlReader::new(reader, opts)),
            TreeFmt::Newick | TreeFmt::ExtendedNewick | TreeFmt::Auto => {
                ReaderKind::Newick(NewickReader::new(reader, opts))
            }
        };
        Ok(Self { inner })
    }

    /// Skip up to `count` trees without parsing topology, returning how many were skipped.
    pub fn skip_trees(&mut self, count: usize) -> Result<usize> {
        match &mut self.inner {
            ReaderKind::Newick(r) => r.skip_trees(count),
            ReaderKind::Nexus(r) => r.skip_trees(count),
            ReaderKind::PhyloXml(r) => r.skip_trees(count),
            ReaderKind::NeXml(r) => r.skip_trees(count),
        }
    }
}

impl TreeReader<Box<dyn BufRead + Send>> {
    /// Open a tree file from disk, transparently decoding `.gz`/`.gzip` or gzip magic bytes.
    pub fn open(path: &Path, fmt: TreeFmt, opts: ParseOptions) -> Result<Self> {
        let reader = open_tree_reader(path)?;
        Self::from_reader(reader, fmt, opts)
    }
}

impl<R: BufRead> Iterator for TreeReader<R> {
    type Item = Result<Tree>;

    fn next(&mut self) -> Option<Self::Item> {
        match &mut self.inner {
            ReaderKind::Newick(r) => r.next(),
            ReaderKind::Nexus(r) => r.next(),
            ReaderKind::PhyloXml(r) => r.next(),
            ReaderKind::NeXml(r) => r.next(),
        }
    }
}

fn open_tree_reader(path: &Path) -> Result<Box<dyn BufRead + Send>> {
    let file = File::open(path).with_context(|| format!("Failed opening {}", path.display()))?;
    let mut buffered = BufReader::with_capacity(IO_BUF_SIZE, file);
    let is_gz = files::is_gzip(path) || buffered.fill_buf()?.starts_with(&[0x1F, 0x8B]);
    if is_gz {
        Ok(Box::new(BufReader::with_capacity(
            IO_BUF_SIZE,
            MultiGzDecoder::new(buffered),
        )))
    } else {
        Ok(Box::new(buffered))
    }
}

/// Parse a single tree from a Newick or NEXUS string (auto-detected).
pub fn parse_tree_str(input: &str) -> Result<Tree> {
    let mut reader = TreeReader::from_reader(
        std::io::Cursor::new(input.as_bytes()),
        TreeFmt::Auto,
        ParseOptions::default(),
    )?;
    reader
        .next()
        .transpose()?
        .context("No trees found in input string")
}

/// Parse all trees from a Newick or NEXUS string (auto-detected).
pub fn parse_trees_str(input: &str) -> Result<Vec<Tree>> {
    TreeReader::from_reader(
        std::io::Cursor::new(input.as_bytes()),
        TreeFmt::Auto,
        ParseOptions::default(),
    )?
    .collect()
}

/// Parse the first tree from a tree file.
pub fn parse_tree_file(path: &Path, fmt: TreeFmt) -> Result<Tree> {
    let mut reader = TreeReader::open(path, fmt, ParseOptions::default())?;
    reader
        .next()
        .transpose()?
        .with_context(|| format!("No trees found in {}", path.display()))
}

/// Parse all trees from a tree file using streaming I/O.
pub fn parse_trees_file(path: &Path, fmt: TreeFmt) -> Result<Vec<Tree>> {
    TreeReader::open(path, fmt, ParseOptions::default())?.collect()
}

/// Parse all trees from a tree file in parallel across threads using `rayon`.
pub fn parse_trees_file_par(path: &Path, fmt: TreeFmt) -> Result<Vec<Tree>> {
    let mut reader = open_tree_reader(path)?;
    let mut bytes = Vec::new();
    reader
        .read_to_end(&mut bytes)
        .with_context(|| format!("Failed reading {}", path.display()))?;
    let resolved = match fmt {
        TreeFmt::Auto => detect_tree_format(&bytes),
        other => other,
    };
    let opts = ParseOptions::default();
    match resolved {
        TreeFmt::Nexus => parse_nexus_bytes_par(&bytes, &opts),
        TreeFmt::PhyloXml => PhyloXmlReader::new(std::io::Cursor::new(&bytes), opts).collect(),
        TreeFmt::NeXml => NeXmlReader::new(std::io::Cursor::new(&bytes), opts).collect(),
        TreeFmt::Newick | TreeFmt::ExtendedNewick | TreeFmt::Auto => {
            parse_newick_bytes_par(&bytes, &opts)
        }
    }
}

/// Discover tree files in a directory independently of [`crate::helper::finder::SeqFileFinder`].
pub struct TreeFileFinder<'a> {
    dir: &'a Path,
}

impl<'a> TreeFileFinder<'a> {
    /// Create a new finder rooted at `dir`.
    pub fn new(dir: &'a Path) -> Self {
        Self { dir }
    }

    /// Find matching tree files directly inside `dir`, sorted alphanumerically.
    pub fn find(&self, fmt: &TreeFmt) -> Result<Vec<PathBuf>> {
        let mut files: Vec<PathBuf> = std::fs::read_dir(self.dir)
            .with_context(|| format!("Failed reading tree directory {}", self.dir.display()))?
            .filter_map(|entry| entry.ok())
            .map(|entry| entry.path())
            .filter(|path| path.is_file() && matches_tree_extension(path, fmt))
            .collect();
        alphanumeric_sort::sort_path_slice(&mut files);
        if files.is_empty() {
            bail!(
                "No tree files found in {} for format '{}'",
                self.dir.display(),
                fmt
            );
        }
        Ok(files)
    }

    /// Find matching tree files recursively under `dir`, sorted alphanumerically.
    pub fn find_recursive(&self, fmt: &TreeFmt) -> Vec<PathBuf> {
        let mut files: Vec<PathBuf> = WalkDir::new(self.dir)
            .into_iter()
            .filter_map(|entry| entry.ok())
            .filter(|entry| entry.file_type().is_file())
            .map(|entry| entry.into_path())
            .filter(|path| matches_tree_extension(path, fmt))
            .collect();
        alphanumeric_sort::sort_path_slice(&mut files);
        files
    }
}

fn matches_tree_extension(path: &Path, fmt: &TreeFmt) -> bool {
    let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
        return false;
    };
    let lower = name.to_ascii_lowercase();
    let stripped = lower
        .strip_suffix(".gz")
        .or_else(|| lower.strip_suffix(".gzip"))
        .unwrap_or(&lower);
    let Some((_, ext)) = stripped.rsplit_once('.') else {
        return false;
    };
    match fmt {
        TreeFmt::Newick => matches!(
            ext,
            "nwk" | "newick" | "tre" | "tree" | "trees" | "treefile" | "nhx"
        ),
        TreeFmt::ExtendedNewick => matches!(
            ext,
            "enwk" | "enewick" | "nwk" | "newick" | "tre" | "tree" | "net"
        ),
        TreeFmt::Nexus => matches!(ext, "nex" | "nxs" | "nexus" | "tre" | "trees" | "t"),
        TreeFmt::PhyloXml => matches!(ext, "phyloxml" | "pxml" | "xml"),
        TreeFmt::NeXml => matches!(ext, "nexml" | "xml"),
        TreeFmt::Auto => matches!(
            ext,
            "nwk"
                | "newick"
                | "enwk"
                | "enewick"
                | "net"
                | "tre"
                | "tree"
                | "trees"
                | "treefile"
                | "nhx"
                | "nex"
                | "nxs"
                | "nexus"
                | "t"
                | "phyloxml"
                | "pxml"
                | "nexml"
                | "xml"
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    use flate2::{Compression, write::GzEncoder};

    use crate::helper::types::DataType;
    use crate::parser::nexus::Nexus;

    #[test]
    fn test_detect_tree_format() {
        assert_eq!(detect_tree_format(b"(A,B);"), TreeFmt::Newick);
        assert_eq!(detect_tree_format(b"  [&R] (A,B);"), TreeFmt::Newick);
        assert_eq!(detect_tree_format(b"#NEXUS\nBEGIN TREES;"), TreeFmt::Nexus);
        assert_eq!(
            detect_tree_format(b"\xEF\xBB\xBF  [comment] #nexus\nBEGIN TREES;"),
            TreeFmt::Nexus
        );
        assert_eq!(
            detect_tree_format(b"<?xml version=\"1.0\"?><!-- c --><phyloxml>"),
            TreeFmt::PhyloXml
        );
        assert_eq!(
            detect_tree_format(
                b"<?xml version=\"1.0\"?><nex:nexml xmlns:nex=\"http://www.nexml.org/2009\">"
            ),
            TreeFmt::NeXml
        );
    }

    #[test]
    fn test_gzip_tree_round_trip_read() {
        let dir = tempfile::TempDir::with_prefix("segul_tree_gz").unwrap();
        let nwk_gz = dir.path().join("sample.nwk.gz");
        let mut enc = GzEncoder::new(File::create(&nwk_gz).unwrap(), Compression::fast());
        enc.write_all(b"(A:0.1,B:0.2);\n(C:0.3,D:0.4);\n").unwrap();
        enc.finish().unwrap();

        let trees = parse_trees_file(&nwk_gz, TreeFmt::Auto).unwrap();
        assert_eq!(trees.len(), 2);
        let par_trees = parse_trees_file_par(&nwk_gz, TreeFmt::Auto).unwrap();
        assert_eq!(par_trees.len(), 2);
        assert_eq!(
            par_trees[1]
                .leaves()
                .filter_map(|id| par_trees[1].label(id))
                .collect::<Vec<_>>(),
            vec!["C", "D"]
        );
    }

    #[test]
    fn test_mixed_data_and_trees_nexus_isolation() {
        let dir = tempfile::TempDir::with_prefix("segul_mixed_nex").unwrap();
        let path = dir.path().join("mixed.nex");
        let content = "#NEXUS\n\
                       BEGIN DATA;\n\
                         DIMENSIONS NTAX=3 NCHAR=4;\n\
                         FORMAT DATATYPE=DNA MISSING=? GAP=-;\n\
                         MATRIX\n\
                           TaxonA ACGT\n\
                           TaxonB TCGA\n\
                           TaxonC AAAA\n\
                         ;\n\
                       END;\n\
                       BEGIN TREES;\n\
                         TRANSLATE\n\
                           1 TaxonA,\n\
                           2 TaxonB,\n\
                           3 TaxonC;\n\
                         TREE t1 = [&R] ((1:0.1,2:0.2):0.3,3:0.4);\n\
                       END;\n";
        std::fs::write(&path, content).unwrap();

        // Sequence NEXUS parser remains completely unaffected.
        let mut seq_nex = Nexus::new(&path, &DataType::Dna);
        seq_nex.parse();
        assert_eq!(seq_nex.matrix.len(), 3);
        assert_eq!(seq_nex.header.ntax, 3);
        assert_eq!(seq_nex.header.nchar, 4);
        assert_eq!(seq_nex.matrix.get("TaxonA"), Some(b"ACGT".as_slice()));

        // Tree NEXUS parser skips DATA and reads TREES cleanly.
        let trees = parse_trees_file(&path, TreeFmt::Auto).unwrap();
        assert_eq!(trees.len(), 1);
        assert_eq!(
            trees[0]
                .leaves()
                .filter_map(|id| trees[0].label(id))
                .collect::<Vec<_>>(),
            vec!["TaxonA", "TaxonB", "TaxonC"]
        );
    }

    #[test]
    fn test_tree_file_finder() {
        let dir = tempfile::TempDir::with_prefix("segul_tree_find").unwrap();
        std::fs::write(dir.path().join("a.nwk"), "(A,B);").unwrap();
        std::fs::write(dir.path().join("b.tre.gz"), "").unwrap();
        std::fs::write(dir.path().join("c.nex"), "#NEXUS").unwrap();
        std::fs::write(dir.path().join("ignore.txt"), "hi").unwrap();

        let finder = TreeFileFinder::new(dir.path());
        assert_eq!(finder.find(&TreeFmt::Auto).unwrap().len(), 3);
        assert_eq!(finder.find(&TreeFmt::Newick).unwrap().len(), 2);
        // `.tre.gz` is shared by Newick and NEXUS tools, plus `c.nex`.
        assert_eq!(finder.find(&TreeFmt::Nexus).unwrap().len(), 2);
    }
}
