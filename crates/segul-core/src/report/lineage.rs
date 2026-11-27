//! Automatic parent-report lineage discovery for chained analyses.

use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use super::RunReport;
use super::checksum::sha256_file;
use super::paths::{relative_posix, resolve_portable_path};
use super::spec::TaskKind;

/// Reference to an upstream `report.json` whose outputs fed into the current run.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ParentReportRef {
    /// Portable POSIX path from the current run's `output_dir` to the parent report.
    pub path: String,
    /// Task kind of the parent run.
    pub kind: TaskKind,
    /// UTC timestamp of the parent run.
    pub timestamp_utc: String,
    /// SHA-256 of the parent `report.json` file itself.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sha256: Option<String>,
    /// Input files (relative to current `base_dir` or filename) matched to the
    /// parent report's outputs.
    pub matched_files: Vec<String>,
}

/// Inspect parent directories of `primary` and `auxiliary` input files for
/// upstream `*report.json` files whose recorded outputs match any input file.
pub fn discover_parent_reports(
    primary: &[PathBuf],
    auxiliary: &[PathBuf],
    output_dir: &Path,
    base_dir: &Path,
    input_hashes: &HashMap<PathBuf, String>,
) -> Vec<ParentReportRef> {
    let mut dirs: BTreeSet<PathBuf> = BTreeSet::new();
    for path in primary.iter().chain(auxiliary.iter()) {
        if let Some(parent) = path.parent() {
            let dir = if parent.as_os_str().is_empty() {
                PathBuf::from(".")
            } else {
                parent.to_path_buf()
            };
            dirs.insert(dir.clone());
            // Also check one level up (e.g. `trimmed_alignments/` inside output dir).
            if let Some(grandparent) = dir.parent()
                && !grandparent.as_os_str().is_empty()
            {
                dirs.insert(grandparent.to_path_buf());
            }
        }
    }

    let mut input_map: HashMap<PathBuf, String> = HashMap::new();
    for path in primary.iter().chain(auxiliary.iter()) {
        let key = canonical_or_normalized(path);
        input_map.insert(key, relative_posix(path, base_dir));
    }

    let mut seen_reports: HashSet<PathBuf> = HashSet::new();
    let mut parents: Vec<ParentReportRef> = Vec::new();
    let canonical_out = canonical_or_normalized(output_dir);

    for dir in dirs {
        let Ok(entries) = fs::read_dir(&dir) else {
            continue;
        };
        for entry in entries.flatten() {
            let candidate = entry.path();
            if !candidate.is_file() || !is_report_filename(&candidate) {
                continue;
            }
            let canon = canonical_or_normalized(&candidate);
            if !seen_reports.insert(canon.clone()) {
                continue;
            }
            // Do not self-link to a report inside the current output directory.
            if canon.parent() == Some(canonical_out.as_path()) {
                continue;
            }
            let Ok(parent_report) = RunReport::from_path(&candidate) else {
                continue;
            };
            let Some(parent_dir) = candidate.parent() else {
                continue;
            };
            let mut matched: BTreeSet<String> = BTreeSet::new();
            for out_rec in &parent_report.outputs.files {
                let resolved_out = resolve_portable_path(parent_dir, &out_rec.path);
                let out_key = canonical_or_normalized(&resolved_out);
                if let Some(rel_input) = input_map.get(&out_key) {
                    // If both have SHA-256 digests, require them to agree.
                    let hash_ok = match (out_rec.sha256.as_ref(), input_hashes.get(&out_key)) {
                        (Some(expected), Some(actual)) => expected == actual,
                        _ => true,
                    };
                    if hash_ok {
                        matched.insert(rel_input.clone());
                    }
                }
            }
            if !matched.is_empty() {
                let sha256 = sha256_file(&candidate).ok();
                parents.push(ParentReportRef {
                    path: relative_posix(&candidate, output_dir),
                    kind: parent_report.kind,
                    timestamp_utc: parent_report.run.timestamp_utc,
                    sha256,
                    matched_files: matched.into_iter().collect(),
                });
            }
        }
    }

    parents.sort_by(|a, b| {
        a.timestamp_utc
            .cmp(&b.timestamp_utc)
            .then_with(|| a.path.cmp(&b.path))
    });
    parents
}

pub(super) fn is_report_filename(path: &Path) -> bool {
    let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
        return false;
    };
    name == "report.json" || name.ends_with("_report.json")
}

pub(super) fn canonical_or_normalized(path: &Path) -> PathBuf {
    fs::canonicalize(path).unwrap_or_else(|_| {
        if path.is_absolute() {
            path.to_path_buf()
        } else {
            std::env::current_dir()
                .map(|cwd| cwd.join(path))
                .unwrap_or_else(|_| path.to_path_buf())
        }
    })
}
