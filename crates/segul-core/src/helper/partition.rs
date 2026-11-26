//! Helper functions for converting alignment partition.
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};

use super::types::PartitionFmt;

/// Construct an output partition file path alongside `input`.
pub fn construct_partition_path(input: &Path, out_part_fmt: &PartitionFmt) -> PathBuf {
    let file_stem = input
        .file_stem()
        .and_then(OsStr::to_str)
        .expect("Failed to parse input file stem");
    let mut fname = PathBuf::from(format!("{}_partition", file_stem));
    match *out_part_fmt {
        PartitionFmt::Charset
        | PartitionFmt::CharsetCodon
        | PartitionFmt::Nexus
        | PartitionFmt::NexusCodon => {
            fname.set_extension("nex");
        }
        PartitionFmt::Raxml | PartitionFmt::RaxmlCodon => {
            fname.set_extension("txt");
        }
    }
    let parent_path = input.parent().expect("Failed to parse input parent path");
    parent_path.join(fname)
}

/// Infer the partition file format from extension and, when ambiguous, file contents.
pub fn infer_partition_fmt(path: &Path) -> Result<PartitionFmt> {
    let ext = path
        .extension()
        .and_then(OsStr::to_str)
        .unwrap_or("")
        .to_ascii_lowercase();
    match ext.as_str() {
        "nex" | "nexus" | "nxs" | "charset" => Ok(PartitionFmt::Nexus),
        "raxml" => Ok(PartitionFmt::Raxml),
        _ => sniff_partition_fmt(path),
    }
}

fn sniff_partition_fmt(path: &Path) -> Result<PartitionFmt> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open partition file {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut has_equals = false;
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with("//") {
            continue;
        }
        let lower = trimmed.to_ascii_lowercase();
        if lower.starts_with("#nexus")
            || lower.starts_with("begin sets")
            || lower.starts_with("begin assumptions")
            || starts_with_keyword(&lower, "charset")
            || starts_with_keyword(&lower, "charpartition")
        {
            return Ok(PartitionFmt::Nexus);
        }
        if trimmed.starts_with('#') {
            continue;
        }
        if trimmed.contains('=') {
            has_equals = true;
        }
    }
    if has_equals {
        Ok(PartitionFmt::Raxml)
    } else {
        bail!(
            "Cannot infer partition format for {}. Please specify the partition format.",
            path.display()
        )
    }
}

fn starts_with_keyword(lower: &str, keyword: &str) -> bool {
    lower
        .strip_prefix(keyword)
        .is_some_and(|rest| rest.starts_with(char::is_whitespace))
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_construct_partition_path_all_formats() {
        let input = Path::new("dir/input.nex");
        assert_eq!(
            construct_partition_path(input, &PartitionFmt::Nexus),
            PathBuf::from("dir/input_partition.nex")
        );
        assert_eq!(
            construct_partition_path(input, &PartitionFmt::Charset),
            PathBuf::from("dir/input_partition.nex")
        );
        assert_eq!(
            construct_partition_path(input, &PartitionFmt::CharsetCodon),
            PathBuf::from("dir/input_partition.nex")
        );
        assert_eq!(
            construct_partition_path(input, &PartitionFmt::RaxmlCodon),
            PathBuf::from("dir/input_partition.txt")
        );
    }

    #[test]
    fn test_infer_partition_fmt() {
        assert_eq!(
            infer_partition_fmt(Path::new("tests/files/partition/partition.nex")).unwrap(),
            PartitionFmt::Nexus
        );
        assert_eq!(
            infer_partition_fmt(Path::new("tests/files/partition/partition.txt")).unwrap(),
            PartitionFmt::Raxml
        );
    }
}
