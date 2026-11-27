//! SHA-256 and file-size verification for `report.json`.

use std::fs;
use std::path::Path;
use std::str::FromStr;

use anyhow::Result;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use super::RunReport;
use super::checksum::sha256_file;
use super::paths::{resolve_base_dir, resolve_portable_path, to_posix_string};

/// Which files inside a `RunReport` to verify.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum VerifyScope {
    #[default]
    All,
    Inputs,
    Outputs,
}

impl std::fmt::Display for VerifyScope {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::All => write!(f, "all"),
            Self::Inputs => write!(f, "inputs"),
            Self::Outputs => write!(f, "outputs"),
        }
    }
}

impl FromStr for VerifyScope {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "all" => Ok(Self::All),
            "inputs" => Ok(Self::Inputs),
            "outputs" => Ok(Self::Outputs),
            _ => Err(format!(
                "'{}' is not a valid verification scope. Valid options: all, inputs, outputs",
                s
            )),
        }
    }
}

/// Verification outcome for a single recorded file.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case")]
pub enum FileVerifyStatus {
    /// File exists, size matches, and SHA-256 matches recorded digest.
    Ok,
    /// File exists but SHA-256 differs from recorded digest.
    Mismatch { expected: String, actual: String },
    /// File does not exist at the resolved path.
    Missing,
    /// File exists (without a recorded SHA-256) but byte length differs.
    SizeMismatch { expected: u64, actual: u64 },
    /// File exists and byte length matches, but no SHA-256 was recorded.
    Unchecked,
}

/// Verification details for one entry from a report manifest.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FileVerifyResult {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub role: Option<String>,
    pub path: String,
    pub resolved_path: String,
    pub expected_size_bytes: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub actual_size_bytes: Option<u64>,
    #[serde(flatten)]
    pub status: FileVerifyStatus,
}

/// Aggregated verification result across inputs, auxiliary files, outputs,
/// and parent report links.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct VerifyReport {
    pub report_path: String,
    pub scope: VerifyScope,
    pub passed: bool,
    pub checked: usize,
    pub ok: usize,
    pub mismatched: usize,
    pub missing: usize,
    pub unchecked: usize,
    pub inputs: Vec<FileVerifyResult>,
    pub auxiliary: Vec<FileVerifyResult>,
    pub outputs: Vec<FileVerifyResult>,
    pub parents: Vec<FileVerifyResult>,
}

impl VerifyReport {
    /// Parse a serialized [`VerifyReport`] from JSON.
    pub fn from_json(json: &str) -> Result<Self> {
        Ok(serde_json::from_str(json)?)
    }

    /// Serialize this verification report as compact JSON.
    pub fn to_json(&self) -> Result<String> {
        Ok(serde_json::to_string(self)?)
    }

    /// Serialize this verification report as pretty-printed JSON.
    pub fn to_json_pretty(&self) -> Result<String> {
        Ok(serde_json::to_string_pretty(self)?)
    }
}

/// Verify files referenced by `report_path`.
///
/// - `input_base`: overrides `report.inputs.base_dir` when resolving inputs.
/// - `output_dir`: overrides `report_path.parent()` when resolving outputs.
pub fn verify_report(
    report_path: &Path,
    input_base: Option<&Path>,
    output_dir: Option<&Path>,
    scope: VerifyScope,
) -> Result<VerifyReport> {
    let report = RunReport::from_path(report_path)?;
    let report_dir = report_path.parent().unwrap_or_else(|| Path::new("."));
    let resolved_input_base = match input_base {
        Some(dir) => dir.to_path_buf(),
        None => resolve_base_dir(&report.inputs.base_dir, None),
    };
    let resolved_output_dir = output_dir
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| report_dir.to_path_buf());

    let check_inputs = matches!(scope, VerifyScope::All | VerifyScope::Inputs);
    let check_outputs = matches!(scope, VerifyScope::All | VerifyScope::Outputs);

    let inputs = if check_inputs {
        let specs: Vec<_> = report
            .inputs
            .primary
            .iter()
            .map(|f| {
                (
                    None,
                    f.path.clone(),
                    resolve_portable_path(&resolved_input_base, &f.path),
                    Some(f.size_bytes),
                    f.sha256.clone(),
                )
            })
            .collect();
        verify_entries(&specs)
    } else {
        Vec::new()
    };

    let auxiliary = if check_inputs {
        let specs: Vec<_> = report
            .inputs
            .auxiliary
            .iter()
            .map(|f| {
                (
                    Some(f.role.clone()),
                    f.path.clone(),
                    resolve_portable_path(&resolved_input_base, &f.path),
                    Some(f.size_bytes),
                    f.sha256.clone(),
                )
            })
            .collect();
        verify_entries(&specs)
    } else {
        Vec::new()
    };

    let outputs = if check_outputs {
        let specs: Vec<_> = report
            .outputs
            .files
            .iter()
            .map(|f| {
                (
                    Some(f.role.clone()),
                    f.path.clone(),
                    resolve_portable_path(&resolved_output_dir, &f.path),
                    Some(f.size_bytes),
                    f.sha256.clone(),
                )
            })
            .collect();
        verify_entries(&specs)
    } else {
        Vec::new()
    };

    let parents = if check_inputs {
        let specs: Vec<_> = report
            .run
            .parent_reports
            .iter()
            .map(|p| {
                (
                    Some(p.kind.to_string()),
                    p.path.clone(),
                    resolve_portable_path(&resolved_output_dir, &p.path),
                    None,
                    p.sha256.clone(),
                )
            })
            .collect();
        verify_entries(&specs)
    } else {
        Vec::new()
    };

    let mut checked = 0usize;
    let mut ok = 0usize;
    let mut mismatched = 0usize;
    let mut missing = 0usize;
    let mut unchecked = 0usize;

    for res in inputs
        .iter()
        .chain(auxiliary.iter())
        .chain(outputs.iter())
        .chain(parents.iter())
    {
        checked += 1;
        match &res.status {
            FileVerifyStatus::Ok => ok += 1,
            FileVerifyStatus::Mismatch { .. } | FileVerifyStatus::SizeMismatch { .. } => {
                mismatched += 1
            }
            FileVerifyStatus::Missing => missing += 1,
            FileVerifyStatus::Unchecked => unchecked += 1,
        }
    }

    Ok(VerifyReport {
        report_path: to_posix_string(report_path),
        scope,
        passed: mismatched == 0 && missing == 0,
        checked,
        ok,
        mismatched,
        missing,
        unchecked,
        inputs,
        auxiliary,
        outputs,
        parents,
    })
}

type EntrySpec = (
    Option<String>,
    String,
    std::path::PathBuf,
    Option<u64>,
    Option<String>,
);

fn verify_entries(specs: &[EntrySpec]) -> Vec<FileVerifyResult> {
    specs
        .par_iter()
        .map(
            |(role, path, resolved, expected_size, expected_hash)| match fs::metadata(resolved) {
                Err(_) => FileVerifyResult {
                    role: role.clone(),
                    path: path.clone(),
                    resolved_path: to_posix_string(resolved),
                    expected_size_bytes: expected_size.unwrap_or(0),
                    actual_size_bytes: None,
                    status: FileVerifyStatus::Missing,
                },
                Ok(meta) => {
                    let actual_size = meta.len();
                    let status = match expected_hash {
                        Some(expected) => match sha256_file(resolved) {
                            Ok(actual) if &actual == expected => FileVerifyStatus::Ok,
                            Ok(actual) => FileVerifyStatus::Mismatch {
                                expected: expected.clone(),
                                actual,
                            },
                            Err(_) => FileVerifyStatus::Missing,
                        },
                        None => match expected_size {
                            Some(exp) if *exp != actual_size => FileVerifyStatus::SizeMismatch {
                                expected: *exp,
                                actual: actual_size,
                            },
                            _ => FileVerifyStatus::Unchecked,
                        },
                    };
                    FileVerifyResult {
                        role: role.clone(),
                        path: path.clone(),
                        resolved_path: to_posix_string(resolved),
                        expected_size_bytes: expected_size.unwrap_or(actual_size),
                        actual_size_bytes: Some(actual_size),
                        status,
                    }
                }
            },
        )
        .collect()
}
