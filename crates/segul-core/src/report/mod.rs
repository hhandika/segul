//! Unified, reproducible `report.json` provenance, statistics, verification,
//! replay, workflow lineage, and publication methods export.

pub mod checksum;
pub mod compat;
pub mod lineage;
pub mod methods;
pub mod paths;
pub mod replay;
pub mod spec;
pub mod stats;
pub mod verify;
pub mod workflow;

#[cfg(test)]
mod tests;

use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime};

use anyhow::{Context, Result};
use chrono::{SecondsFormat, Utc};
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::helper::utils;

pub use checksum::ChecksumMode;
pub use lineage::ParentReportRef;
pub use spec::{TaskKind, TaskParameters};
pub use stats::ReportStatistics;

/// Current schema version for the unified `report.json` envelope.
pub const SCHEMA_VERSION: u32 = 2;

/// Default file name of the unified JSON report.
pub const DEFAULT_REPORT_NAME: &str = "report.json";

/// Build the report path inside `output_dir`, prepending `prefix` if present.
pub fn report_path(output_dir: &Path, prefix: Option<&str>) -> PathBuf {
    match prefix.filter(|p| !p.is_empty()) {
        Some(prefix) => output_dir.join(format!("{}_{}", prefix, DEFAULT_REPORT_NAME)),
        None => output_dir.join(DEFAULT_REPORT_NAME),
    }
}

/// Frontend application that invoked `segul-core`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FrontendInfo {
    pub name: String,
    pub version: String,
}

impl FrontendInfo {
    pub fn cli(version: impl Into<String>) -> Self {
        Self {
            name: String::from("segul-cli"),
            version: version.into(),
        }
    }

    pub fn gui(version: impl Into<String>) -> Self {
        Self {
            name: String::from("segui"),
            version: version.into(),
        }
    }

    pub fn python(version: impl Into<String>) -> Self {
        Self {
            name: String::from("pysegul"),
            version: version.into(),
        }
    }

    pub fn core() -> Self {
        Self {
            name: String::from("segul-core"),
            version: utils::get_crate_version(),
        }
    }
}

impl Default for FrontendInfo {
    fn default() -> Self {
        Self::core()
    }
}

/// Host platform summary recorded for reproducibility.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PlatformInfo {
    pub os: String,
    pub arch: String,
    pub threads: usize,
}

impl PlatformInfo {
    pub fn current() -> Self {
        Self {
            os: std::env::consts::OS.to_string(),
            arch: std::env::consts::ARCH.to_string(),
            threads: rayon::current_num_threads(),
        }
    }
}

/// How primary inputs were supplied by the user.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum SourceMode {
    #[default]
    Files,
    Dir,
}

/// Context passed by CLI, SEGUI, or PySEGUL into `segul-core` tasks.
#[derive(Debug, Clone)]
pub struct ReportContext {
    pub frontend: FrontendInfo,
    pub checksum_mode: ChecksumMode,
    pub source_mode: SourceMode,
    pub input_dir: Option<PathBuf>,
}

impl Default for ReportContext {
    fn default() -> Self {
        Self {
            frontend: FrontendInfo::core(),
            checksum_mode: ChecksumMode::from_env(),
            source_mode: SourceMode::Files,
            input_dir: None,
        }
    }
}

impl ReportContext {
    pub fn new(frontend: FrontendInfo, checksum_mode: ChecksumMode) -> Self {
        Self {
            frontend,
            checksum_mode,
            source_mode: SourceMode::Files,
            input_dir: None,
        }
    }

    pub fn with_input_dir(mut self, dir: Option<PathBuf>) -> Self {
        if let Some(d) = dir {
            self.source_mode = SourceMode::Dir;
            self.input_dir = Some(d);
        }
        self
    }
}

/// Metadata for a primary input file.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FileRecord {
    pub path: String,
    pub size_bytes: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sha256: Option<String>,
}

/// Metadata for a secondary input file (BED, partition, group map, ID list, etc.).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AuxFileRecord {
    pub role: String,
    pub path: String,
    pub size_bytes: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sha256: Option<String>,
}

/// Metadata for a written output file.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct OutputFileRecord {
    pub role: String,
    pub path: String,
    pub size_bytes: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sha256: Option<String>,
}

/// Input manifest inside `report.json`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct InputManifest {
    pub source_mode: SourceMode,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub input_dir: Option<String>,
    pub base_dir: String,
    pub primary: Vec<FileRecord>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub auxiliary: Vec<AuxFileRecord>,
}

/// Output manifest inside `report.json`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct OutputManifest {
    pub output_dir: String,
    pub files: Vec<OutputFileRecord>,
}

/// Provenance and execution metadata inside `report.json`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RunMetadata {
    pub segul_core_version: String,
    pub frontend: FrontendInfo,
    pub timestamp_utc: String,
    pub elapsed_ms: u64,
    pub platform: PlatformInfo,
    pub checksum_mode: ChecksumMode,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub checksum_algorithm: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub parent_reports: Vec<ParentReportRef>,
}

/// Unified Schema v2 JSON report emitted by every SEGUL operation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RunReport {
    pub schema_version: u32,
    pub kind: TaskKind,
    pub run: RunMetadata,
    pub inputs: InputManifest,
    pub parameters: Value,
    pub outputs: OutputManifest,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub statistics: Option<Value>,
}

impl RunReport {
    /// Load a v2 `RunReport` from disk.
    pub fn from_path(path: &Path) -> Result<Self> {
        let text = fs::read_to_string(path)
            .with_context(|| format!("Failed reading {}", path.display()))?;
        Self::from_json(&text).with_context(|| format!("Failed parsing {}", path.display()))
    }

    /// Parse a v2 `RunReport` from a JSON string.
    pub fn from_json(json: &str) -> Result<Self> {
        let report: Self = serde_json::from_str(json)?;
        if report.schema_version < 2 {
            anyhow::bail!(
                "Unsupported report schema version {}",
                report.schema_version
            );
        }
        Ok(report)
    }

    pub fn to_json(&self) -> Result<String> {
        Ok(serde_json::to_string(self)?)
    }

    pub fn to_json_pretty(&self) -> Result<String> {
        Ok(serde_json::to_string_pretty(self)?)
    }

    /// Write formatted JSON to `path`, creating parent directories if needed.
    pub fn write(&self, path: &Path) -> Result<()> {
        if let Some(parent) = path.parent()
            && !parent.as_os_str().is_empty()
        {
            fs::create_dir_all(parent)
                .with_context(|| format!("Failed creating {}", parent.display()))?;
        }
        fs::write(path, self.to_json_pretty()?)
            .with_context(|| format!("Failed writing {}", path.display()))
    }

    /// Parse `parameters` into a strongly typed [`TaskParameters`] variant.
    pub fn typed_parameters(&self) -> Result<TaskParameters> {
        TaskParameters::from_json_value(self.kind, &self.parameters)
    }

    /// Parse `statistics` into a strongly typed [`ReportStatistics`] variant.
    pub fn typed_statistics(&self) -> Result<Option<ReportStatistics>> {
        match &self.statistics {
            Some(val) if !val.is_null() => {
                ReportStatistics::from_json_value(self.kind, val).map(Some)
            }
            _ => Ok(None),
        }
    }

    /// Append additional output files (such as CLI/PySEGUL `--plot` charts) and
    /// preserve deterministic ordering.
    pub fn append_outputs(&mut self, extra: &[(&str, PathBuf)], output_dir: &Path) -> Result<()> {
        if extra.is_empty() {
            return Ok(());
        }
        let hash_outputs = self.run.checksum_mode.includes_outputs();
        let records = build_output_records(extra, output_dir, hash_outputs)?;
        for rec in records {
            if !self
                .outputs
                .files
                .iter()
                .any(|existing| existing.path == rec.path && existing.role == rec.role)
            {
                self.outputs.files.push(rec);
            }
        }
        self.outputs
            .files
            .sort_by(|a, b| alphanumeric_sort::compare_str(&a.path, &b.path));
        Ok(())
    }

    /// Check whether `path` is a v2 `RunReport` JSON file.
    pub fn is_run_report_file(path: &Path) -> bool {
        Self::from_path(path).is_ok()
    }
}

/// Builder that collects inputs, parameters, outputs, checksums, lineage, and
/// statistics to write a [`RunReport`].
pub struct RunReportBuilder {
    kind: TaskKind,
    ctx: ReportContext,
    started_at: Instant,
    primary: Vec<PathBuf>,
    auxiliary: Vec<(String, PathBuf)>,
    parameters: Value,
    output_dir: PathBuf,
    outputs: Vec<(String, PathBuf)>,
    statistics: Option<Value>,
}

impl RunReportBuilder {
    pub fn new(kind: TaskKind, ctx: &ReportContext) -> Self {
        Self {
            kind,
            ctx: ctx.clone(),
            started_at: Instant::now(),
            primary: Vec::new(),
            auxiliary: Vec::new(),
            parameters: Value::Object(Default::default()),
            output_dir: PathBuf::from("."),
            outputs: Vec::new(),
            statistics: None,
        }
    }

    pub fn primary_inputs(mut self, files: &[PathBuf]) -> Self {
        self.primary = files.to_vec();
        alphanumeric_sort::sort_path_slice(&mut self.primary);
        self
    }

    pub fn auxiliary_inputs<I, S>(mut self, aux: I) -> Self
    where
        I: IntoIterator<Item = (S, PathBuf)>,
        S: Into<String>,
    {
        self.auxiliary = aux.into_iter().map(|(r, p)| (r.into(), p)).collect();
        self.auxiliary.sort_by(|a, b| {
            a.0.cmp(&b.0)
                .then_with(|| alphanumeric_sort::compare_path(&a.1, &b.1))
        });
        self
    }

    pub fn parameters<T: Serialize>(mut self, params: &T) -> Result<Self> {
        self.parameters = serde_json::to_value(params).context("Failed serializing parameters")?;
        Ok(self)
    }

    pub fn outputs<I, S>(mut self, output_dir: &Path, files: I) -> Self
    where
        I: IntoIterator<Item = (S, PathBuf)>,
        S: Into<String>,
    {
        self.output_dir = output_dir.to_path_buf();
        self.outputs = files.into_iter().map(|(r, p)| (r.into(), p)).collect();
        self
    }

    pub fn statistics<T: Serialize>(mut self, stats: &T) -> Result<Self> {
        self.statistics = Some(serde_json::to_value(stats).context("Failed serializing stats")?);
        Ok(self)
    }

    pub fn statistics_value(mut self, stats: Value) -> Self {
        self.statistics = Some(stats);
        self
    }

    /// Build the [`RunReport`], computing any requested SHA-256 digests and
    /// parent lineage links.
    pub fn build(self) -> Result<RunReport> {
        let elapsed_ms = self.started_at.elapsed().as_millis() as u64;
        let timestamp_utc = Utc::now().to_rfc3339_opts(SecondsFormat::Secs, true);
        let base_dir = paths::common_base_dir(&self.primary);
        let hash_inputs = self.ctx.checksum_mode.includes_inputs();
        let hash_outputs = self.ctx.checksum_mode.includes_outputs();

        let primary_hashes = if hash_inputs {
            Some(checksum::sha256_files(&self.primary)?)
        } else {
            None
        };
        let mut input_hash_map: HashMap<PathBuf, String> = HashMap::new();
        let mut primary_records = Vec::with_capacity(self.primary.len());
        for (idx, path) in self.primary.iter().enumerate() {
            let size_bytes = file_size(path);
            let sha256 = primary_hashes.as_ref().map(|h| h[idx].clone());
            if let Some(hash) = &sha256 {
                input_hash_map.insert(lineage::canonical_or_normalized(path), hash.clone());
            }
            primary_records.push(FileRecord {
                path: paths::relative_posix(path, &base_dir),
                size_bytes,
                sha256,
            });
        }

        let aux_paths: Vec<PathBuf> = self.auxiliary.iter().map(|(_, p)| p.clone()).collect();
        let aux_hashes = if hash_inputs {
            Some(checksum::sha256_files(&aux_paths)?)
        } else {
            None
        };
        let mut aux_records = Vec::with_capacity(self.auxiliary.len());
        for (idx, (role, path)) in self.auxiliary.iter().enumerate() {
            let size_bytes = file_size(path);
            let sha256 = aux_hashes.as_ref().map(|h| h[idx].clone());
            if let Some(hash) = &sha256 {
                input_hash_map.insert(lineage::canonical_or_normalized(path), hash.clone());
            }
            aux_records.push(AuxFileRecord {
                role: role.clone(),
                path: paths::relative_posix(path, &base_dir),
                size_bytes,
                sha256,
            });
        }

        let output_pairs: Vec<(&str, PathBuf)> = self
            .outputs
            .iter()
            .filter(|(_, p)| !lineage::is_report_filename(p))
            .map(|(r, p)| (r.as_str(), p.clone()))
            .collect();
        let mut output_records =
            build_output_records(&output_pairs, &self.output_dir, hash_outputs)?;
        output_records.sort_by(|a, b| alphanumeric_sort::compare_str(&a.path, &b.path));

        let parent_reports = lineage::discover_parent_reports(
            &self.primary,
            &aux_paths,
            &self.output_dir,
            &base_dir,
            &input_hash_map,
        );

        let inputs = InputManifest {
            source_mode: self.ctx.source_mode,
            input_dir: self.ctx.input_dir.as_deref().map(paths::to_posix_string),
            base_dir: paths::to_posix_string(&base_dir),
            primary: primary_records,
            auxiliary: aux_records,
        };
        let outputs = OutputManifest {
            output_dir: paths::to_posix_string(&self.output_dir),
            files: output_records,
        };
        let run = RunMetadata {
            segul_core_version: utils::get_crate_version(),
            frontend: self.ctx.frontend,
            timestamp_utc,
            elapsed_ms,
            platform: PlatformInfo::current(),
            checksum_mode: self.ctx.checksum_mode,
            checksum_algorithm: self.ctx.checksum_mode.algorithm_name().map(str::to_string),
            parent_reports,
        };

        Ok(RunReport {
            schema_version: SCHEMA_VERSION,
            kind: self.kind,
            run,
            inputs,
            parameters: self.parameters,
            outputs,
            statistics: self.statistics,
        })
    }

    /// Build the report and write it to `report_path`.
    pub fn write_to(self, report_path: &Path) -> Result<RunReport> {
        let report = self.build()?;
        report.write(report_path)?;
        Ok(report)
    }
}

/// Snapshot of files inside an output directory before running a transform task,
/// so only files created or modified by the run are recorded in `outputs.files`.
#[derive(Debug, Default)]
pub struct DirSnapshot {
    entries: BTreeMap<PathBuf, (u64, Option<SystemTime>)>,
}

impl DirSnapshot {
    /// Capture existing files in `dir` (empty snapshot if `dir` does not exist yet).
    pub fn capture(dir: &Path) -> Self {
        let mut entries = BTreeMap::new();
        if dir.is_dir() {
            for path in walk_files(dir) {
                if let Ok(meta) = fs::metadata(&path) {
                    entries.insert(path, (meta.len(), meta.modified().ok()));
                }
            }
        }
        Self { entries }
    }

    /// Return sorted files in `dir` that were created or modified since `capture`,
    /// ignoring any report JSON files.
    pub fn written_files(&self, dir: &Path) -> Vec<PathBuf> {
        let mut written = Vec::new();
        if !dir.is_dir() {
            return written;
        }
        for path in walk_files(dir) {
            if lineage::is_report_filename(&path) {
                continue;
            }
            let Ok(meta) = fs::metadata(&path) else {
                continue;
            };
            let current = (meta.len(), meta.modified().ok());
            let changed = match self.entries.get(&path) {
                None => true,
                Some(prev) => prev != &current,
            };
            if changed {
                written.push(path);
            }
        }
        alphanumeric_sort::sort_path_slice(&mut written);
        written
    }
}

/// Convenience helper for transform/utility tasks to write `report.json`
/// after capturing outputs via [`DirSnapshot`] or explicit file lists.
#[allow(clippy::too_many_arguments)]
pub fn write_task_report<P: Serialize, S: Serialize>(
    ctx: &ReportContext,
    kind: TaskKind,
    primary_inputs: &[PathBuf],
    auxiliary_inputs: Vec<(&str, PathBuf)>,
    params: &P,
    output_dir: &Path,
    output_files: Vec<(&str, PathBuf)>,
    stats: Option<&S>,
    prefix: Option<&str>,
) -> Result<(PathBuf, RunReport)> {
    let mut builder = RunReportBuilder::new(kind, ctx)
        .primary_inputs(primary_inputs)
        .auxiliary_inputs(auxiliary_inputs)
        .parameters(params)?
        .outputs(output_dir, output_files);
    if let Some(stats) = stats {
        builder = builder.statistics(stats)?;
    }
    let path = report_path(output_dir, prefix);
    let report = builder.write_to(&path)?;
    Ok((path, report))
}

fn build_output_records(
    files: &[(&str, PathBuf)],
    output_dir: &Path,
    compute_sha256: bool,
) -> Result<Vec<OutputFileRecord>> {
    let paths: Vec<PathBuf> = files.iter().map(|(_, p)| p.clone()).collect();
    let hashes = if compute_sha256 {
        Some(checksum::sha256_files(&paths)?)
    } else {
        None
    };
    let mut records = Vec::with_capacity(files.len());
    for (idx, (role, path)) in files.iter().enumerate() {
        records.push(OutputFileRecord {
            role: (*role).to_string(),
            path: paths::relative_posix(path, output_dir),
            size_bytes: file_size(path),
            sha256: hashes.as_ref().map(|h| h[idx].clone()),
        });
    }
    Ok(records)
}

fn file_size(path: &Path) -> u64 {
    fs::metadata(path).map(|m| m.len()).unwrap_or(0)
}

fn walk_files(dir: &Path) -> Vec<PathBuf> {
    walkdir::WalkDir::new(dir)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().is_file())
        .map(|e| e.into_path())
        .collect()
}
