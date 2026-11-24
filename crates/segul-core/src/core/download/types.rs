//! Shared types, errors, and progress contracts for archive downloads.

use std::fmt;
use std::path::PathBuf;
use std::sync::LazyLock;

use regex::Regex;

static SUFFIX_VERSION_RE: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"\.\d+$").expect("Valid version regex"));
static UNSAFE_NAME_RE: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"[^A-Za-z0-9_.]+").expect("Valid safe name regex"));
static MULTI_UNDER_RE: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"_{2,}").expect("Valid underscore regex"));
static EDGE_PUNCT_RE: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^[_.]+|[_.]+$").expect("Valid edge trim regex"));

/// Public sequence archive supported by the downloader.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DatabaseSource {
    GenBank,
    Ena,
}

impl DatabaseSource {
    /// Human-readable archive name used in logs and `sources.tsv`.
    pub fn label(self) -> &'static str {
        match self {
            Self::GenBank => "GenBank",
            Self::Ena => "ENA",
        }
    }
}

/// Kind of biological data downloaded from an archive.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DownloadContent {
    Nucleotide,
    Coding,
    Protein,
    Reads,
    Assembly,
}

impl DownloadContent {
    /// Lowercase identifier matching SEGUI's record keys.
    pub fn key_prefix(self) -> &'static str {
        match self {
            Self::Nucleotide => "nucleotide",
            Self::Coding => "coding",
            Self::Protein => "protein",
            Self::Reads => "reads",
            Self::Assembly => "assembly",
        }
    }

    /// Fallback label for output folder naming.
    pub fn folder_label(self) -> &'static str {
        match self {
            Self::Nucleotide => "Whole record",
            Self::Coding => "Coding sequence",
            Self::Protein => "Protein",
            Self::Reads => "Reads",
            Self::Assembly => "Assembly",
        }
    }

    /// Whether records are written as FASTA with rewritten headers.
    pub fn is_sequence(self) -> bool {
        matches!(self, Self::Nucleotide | Self::Coding | Self::Protein)
    }
}

/// How downloaded sequences are split across output files.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SequenceGrouping {
    SingleFile,
    PerRecord,
}

/// How downloaded FASTA headers are constructed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SequenceNaming {
    SpeciesAccession,
    Accession,
}

/// Current phase of an active download job.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DownloadStage {
    Connecting,
    Downloading,
    Retrying,
    Verifying,
    Decompressing,
    Preparing,
}

/// Per-record status reported during a download.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecordDownloadStatus {
    Pending,
    Active,
    Done,
    Failed,
}

/// A remote file associated with a read run or assembly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RemoteFile {
    pub url: String,
    pub bytes: Option<u64>,
    pub md5: Option<String>,
}

/// A selected archive record queued for download.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DownloadRecord {
    pub source: DatabaseSource,
    pub content: DownloadContent,
    pub accession: String,
    pub organism: String,
    pub title: String,
    pub length: Option<u64>,
    pub codon_start: u32,
    pub gene: String,
    pub files: Vec<RemoteFile>,
}

impl DownloadRecord {
    /// Unique selection key matching SEGUI's `DatabaseRecord.key`.
    pub fn key(&self) -> String {
        format!(
            "{}:{}",
            self.content.key_prefix(),
            base_accession(&self.accession)
        )
    }

    /// Estimated download size when known.
    pub fn estimated_bytes(&self) -> Option<u64> {
        if !self.files.is_empty() && self.files.iter().all(|f| f.bytes.is_some()) {
            return Some(self.files.iter().filter_map(|f| f.bytes).sum());
        }
        if self.content == DownloadContent::Coding && self.source != DatabaseSource::Ena {
            return None;
        }
        self.length
    }
}

/// User options controlling output grouping, naming, and gene filtering.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct DownloadOptions {
    pub grouping: Option<SequenceGrouping>,
    pub naming: Option<SequenceNaming>,
    pub gene: String,
    pub file_name: String,
}

impl DownloadOptions {
    pub fn effective_grouping(&self) -> SequenceGrouping {
        self.grouping.unwrap_or(SequenceGrouping::SingleFile)
    }

    pub fn effective_naming(&self) -> SequenceNaming {
        self.naming.unwrap_or(SequenceNaming::SpeciesAccession)
    }
}

/// Client and NCBI identification settings for a download job.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DownloadConfig {
    pub ncbi_api_key: Option<String>,
    pub tool_name: String,
    pub contact_email: String,
}

impl Default for DownloadConfig {
    fn default() -> Self {
        Self {
            ncbi_api_key: None,
            tool_name: "segui".to_string(),
            contact_email: "herubiolog@gmail.com".to_string(),
        }
    }
}

/// Snapshot of download progress, throughput, stage, and per-record states.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DownloadProgress {
    pub done: u32,
    pub total: u32,
    pub bytes: u64,
    pub bytes_total: u64,
    pub bytes_per_second: Option<u64>,
    pub stage: DownloadStage,
    pub retry_attempt: Option<u32>,
    pub retry_max: Option<u32>,
    pub statuses: Vec<(String, RecordDownloadStatus)>,
}

/// Reason a single record failed to download.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DownloadFailure {
    pub record_key: String,
    pub accession: String,
    pub reason: String,
}

/// Summary returned after `DownloadJob::run`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DownloadReport {
    pub completed: u32,
    pub failures: Vec<DownloadFailure>,
}

impl DownloadReport {
    pub fn is_complete(&self) -> bool {
        self.failures.is_empty()
    }
}

/// Committed output files and inferred SEGUL data type.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompletedDownload {
    pub files: Vec<PathBuf>,
    pub data_type: Option<String>,
}

/// Error returned when a download is cancelled via [`DownloadProgressReporter`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DownloadCancelled;

impl fmt::Display for DownloadCancelled {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Stopped.")
    }
}

impl std::error::Error for DownloadCancelled {}

/// User-actionable failure from an archive request or validation step.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DownloadError(pub String);

impl fmt::Display for DownloadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for DownloadError {}

/// Receives progress snapshots and cancellation queries from [`crate::core::download::DownloadJob`].
pub trait DownloadProgressReporter: Send + Sync {
    fn on_progress(&self, _progress: DownloadProgress) {}

    fn is_cancelled(&self) -> bool {
        false
    }
}

/// Returns `accession` without its `.N` version suffix, in uppercase.
pub fn base_accession(accession: &str) -> String {
    SUFFIX_VERSION_RE
        .replace(accession.trim(), "")
        .to_uppercase()
}

/// Returns the `.N` version suffix of `accession`, if present.
pub fn accession_version(accession: &str) -> Option<u32> {
    let trimmed = accession.trim();
    let dot = trimmed.rfind('.')?;
    let suffix = &trimmed[dot + 1..];
    if suffix.is_empty() || !suffix.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    suffix.parse().ok()
}

/// Whether `found` matches `requested`. A request without a version matches any version.
pub fn accession_matches(requested: &str, found: &str) -> bool {
    if base_accession(requested) != base_accession(found) {
        return false;
    }
    match accession_version(requested) {
        Some(version) => Some(version) == accession_version(found),
        None => true,
    }
}

/// Sanitizes `value` so it is safe as a FASTA/NEXUS/PHYLIP sequence ID and filename stem.
pub fn safe_name(value: &str) -> String {
    let step1 = UNSAFE_NAME_RE.replace_all(value, "_");
    let step2 = MULTI_UNDER_RE.replace_all(&step1, "_");
    let step3 = EDGE_PUNCT_RE.replace_all(&step2, "");
    if step3.is_empty() {
        "sequence".to_string()
    } else {
        step3.into_owned()
    }
}

/// Splits `records` into batches bounded by both `max_count` and cumulative `length` `max_bytes`.
///
/// Always includes at least one record per batch, even if its length exceeds `max_bytes`.
pub fn adaptive_batches(
    records: &[DownloadRecord],
    max_count: usize,
    max_bytes: u64,
) -> Vec<Vec<DownloadRecord>> {
    let max_count = max_count.max(1);
    let mut batches = Vec::new();
    let mut current = Vec::new();
    let mut current_bytes: u64 = 0;

    for record in records {
        let len = record.length.unwrap_or(1_000);
        let exceeds_count = current.len() >= max_count;
        let exceeds_bytes = !current.is_empty() && current_bytes.saturating_add(len) > max_bytes;
        if exceeds_count || exceeds_bytes {
            batches.push(std::mem::take(&mut current));
            current_bytes = 0;
        }
        current_bytes = current_bytes.saturating_add(len);
        current.push(record.clone());
    }

    if !current.is_empty() {
        batches.push(current);
    }
    batches
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_record(acc: &str, length: Option<u64>) -> DownloadRecord {
        DownloadRecord {
            source: DatabaseSource::GenBank,
            content: DownloadContent::Nucleotide,
            accession: acc.to_string(),
            organism: "Mus musculus".to_string(),
            title: String::new(),
            length,
            codon_start: 1,
            gene: String::new(),
            files: Vec::new(),
        }
    }

    #[test]
    fn test_accession_helpers_and_safe_name() {
        assert_eq!(base_accession(" u49845.1 "), "U49845");
        assert_eq!(accession_version("U49845.12"), Some(12));
        assert_eq!(accession_version("U49845"), None);
        assert!(accession_matches("U49845", "U49845.1"));
        assert!(accession_matches("U49845.1", "U49845.1"));
        assert!(!accession_matches("U49845.2", "U49845.1"));
        assert_eq!(safe_name("Rattus sp. 'NH-2147'"), "Rattus_sp._NH_2147");
        assert_eq!(safe_name("../../etc/passwd"), "etc_passwd");
        assert_eq!(safe_name("///"), "sequence");
    }

    #[test]
    fn test_adaptive_batches_packs_small_and_splits_large() {
        let small: Vec<_> = (0..5)
            .map(|i| sample_record(&format!("A{i}.1"), Some(500)))
            .collect();
        let small_batches = adaptive_batches(&small, 3, 5_000_000);
        assert_eq!(small_batches.len(), 2);
        assert_eq!(small_batches[0].len(), 3);
        assert_eq!(small_batches[1].len(), 2);

        let mixed = vec![
            sample_record("BIG1.1", Some(4_000_000)),
            sample_record("BIG2.1", Some(3_000_000)),
            sample_record("HUGE.1", Some(10_000_000)),
            sample_record("TINY.1", Some(500)),
        ];
        let split = adaptive_batches(&mixed, 200, 5_000_000);
        assert_eq!(split.len(), 4);
        assert_eq!(split[0][0].accession, "BIG1.1");
        assert_eq!(split[1][0].accession, "BIG2.1");
        assert_eq!(split[2][0].accession, "HUGE.1");
        assert_eq!(split[3][0].accession, "TINY.1");
    }
}
