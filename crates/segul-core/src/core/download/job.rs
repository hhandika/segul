//! Concurrent, resumable download job orchestration for GenBank and ENA.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use anyhow::Result;
use chrono::Utc;
use regex::Regex;
use reqwest::Url;
use std::sync::LazyLock;
use tokio::sync::Semaphore;
use tokio::task::JoinSet;

use super::assembly::{decompress_fasta, parse_sequence_report, write_ena_assembly};
use super::client::{ArchiveClient, CancelSignal, SaveOptions};
use super::fasta::{PreparedBatch, prepare_batch};
use super::types::{
    CompletedDownload, DatabaseSource, DownloadCancelled, DownloadConfig, DownloadContent,
    DownloadError, DownloadFailure, DownloadOptions, DownloadProgress, DownloadProgressReporter,
    DownloadRecord, DownloadReport, DownloadStage, RecordDownloadStatus, RemoteFile,
    SequenceGrouping, adaptive_batches, base_accession, safe_name,
};

const GENBANK_BATCH_COUNT: usize = 200;
const ENA_BATCH_COUNT: usize = 100;
const MAX_BATCH_BYTES: u64 = 5 * 1024 * 1024;

static FASTQ_EXT_RE: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"\.(fastq|fq)\.gz$").expect("Valid fastq suffix regex"));
static WHITESPACE_RE: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"\s").expect("Valid whitespace regex"));

const SOURCE_COLUMNS: [&str; 8] = [
    "database",
    "accession",
    "file",
    "id",
    "original_header",
    "trimmed_bases",
    "url",
    "downloaded",
];

#[derive(Debug, Clone)]
struct Endpoints {
    genbank_efetch: String,
    ena_fasta_base: String,
    ena_report_base: String,
}

impl Default for Endpoints {
    fn default() -> Self {
        Self {
            genbank_efetch: "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi".to_string(),
            ena_fasta_base: "https://www.ebi.ac.uk/ena/browser/api/fasta".to_string(),
            ena_report_base: "https://ftp.ebi.ac.uk/pub/databases/ena/assembly".to_string(),
        }
    }
}

/// Manages a staged download of selected GenBank and ENA records.
///
/// [`DownloadJob::run`] can be repeated to retry only failed records. Output files
/// remain hidden in a temporary staging folder until [`DownloadJob::commit`].
pub struct DownloadJob {
    client: Arc<ArchiveClient>,
    records: Vec<DownloadRecord>,
    options: DownloadOptions,
    config: DownloadConfig,
    trims_to_frame: bool,
    root: PathBuf,
    taken_names: HashSet<String>,
    endpoints: Endpoints,
    staging: Option<PathBuf>,
    done: HashSet<String>,
    used_ids: HashSet<String>,
    used_names: HashSet<String>,
    file_names: Vec<String>,
    name_of: HashMap<String, String>,
    sources: Vec<[String; 8]>,
    single_name: Option<String>,
    has_unknown_residues: bool,
    temporary_seq: u64,
}

impl DownloadJob {
    /// Creates a new download job using the production NCBI and ENA endpoints.
    pub fn new(
        records: Vec<DownloadRecord>,
        options: DownloadOptions,
        config: DownloadConfig,
        trims_to_frame: bool,
        root: PathBuf,
        taken_names: HashSet<String>,
    ) -> Result<Self> {
        let has_key = config
            .ncbi_api_key
            .as_deref()
            .is_some_and(|k| !k.trim().is_empty());
        let client = Arc::new(ArchiveClient::with_profile(has_key, config.concurrency)?);
        Ok(Self::with_client(
            client,
            records,
            options,
            config,
            trims_to_frame,
            root,
            taken_names,
            Endpoints::default(),
        ))
    }

    #[allow(clippy::too_many_arguments)]
    fn with_client(
        client: Arc<ArchiveClient>,
        records: Vec<DownloadRecord>,
        options: DownloadOptions,
        config: DownloadConfig,
        trims_to_frame: bool,
        root: PathBuf,
        taken_names: HashSet<String>,
        endpoints: Endpoints,
    ) -> Self {
        Self {
            client,
            records,
            options,
            config,
            trims_to_frame,
            root,
            taken_names,
            endpoints,
            staging: None,
            done: HashSet::new(),
            used_ids: HashSet::new(),
            used_names: HashSet::new(),
            file_names: Vec::new(),
            name_of: HashMap::new(),
            sources: Vec::new(),
            single_name: None,
            has_unknown_residues: false,
            temporary_seq: 0,
        }
    }

    /// Executes the download for all pending records, reporting live progress.
    pub async fn run(
        &mut self,
        reporter: Option<Arc<dyn DownloadProgressReporter>>,
    ) -> Result<DownloadReport> {
        let cancel = CancelSignal::new();
        let watcher = spawn_cancel_watcher(reporter.clone(), &cancel);
        let outcome = self.run_inner(reporter, &cancel).await;
        if let Some(handle) = watcher {
            handle.abort();
        }
        outcome
    }

    /// Finalizes `sources.tsv` and renames the staging folder into `root`.
    pub fn commit(&mut self) -> Result<CompletedDownload> {
        let Some(staging) = self.staging.as_ref() else {
            return Err(DownloadError("Nothing was downloaded.".to_string()).into());
        };
        if self.file_names.is_empty() {
            return Err(DownloadError("Nothing was downloaded.".to_string()).into());
        }
        let sources_path = staging.join("sources.tsv");
        {
            let mut writer = BufWriter::new(File::create(&sources_path)?);
            writeln!(writer, "{}", SOURCE_COLUMNS.join("\t"))?;
            for row in &self.sources {
                let clean: Vec<String> = row
                    .iter()
                    .map(|cell| WHITESPACE_RE.replace_all(cell, " ").into_owned())
                    .collect();
                writeln!(writer, "{}", clean.join("\t"))?;
            }
            writer.flush()?;
        }

        let stamp = Utc::now().format("%Y%m%d-%H%M%S").to_string();
        let fallback = self
            .records
            .first()
            .map(|r| r.content.folder_label())
            .unwrap_or("download");
        let raw_label = self.single_name.as_deref().unwrap_or(fallback);
        let stem = Path::new(raw_label)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or(raw_label);
        let label = safe_name(stem);

        let mut target = self.root.join(format!("{label}-{stamp}"));
        let mut idx = 2u32;
        while target.exists() {
            target = self.root.join(format!("{label}-{stamp}-{idx}"));
            idx += 1;
        }
        fs::rename(staging, &target)?;
        self.staging = None;

        let files = self.file_names.iter().map(|n| target.join(n)).collect();
        let data_type = match self.records.first().map(|r| r.content) {
            Some(c) if !c.is_sequence() => None,
            _ if self.has_unknown_residues => Some("Ignore".to_string()),
            Some(DownloadContent::Protein) => Some("AA".to_string()),
            _ => Some("DNA".to_string()),
        };
        Ok(CompletedDownload { files, data_type })
    }

    /// Removes the staging directory and clears job state.
    pub fn discard(&mut self) {
        let staging = self.staging.take();
        self.done.clear();
        self.used_ids.clear();
        self.used_names.clear();
        self.file_names.clear();
        self.name_of.clear();
        self.sources.clear();
        self.single_name = None;
        if let Some(dir) = staging {
            let _ = fs::remove_dir_all(dir);
        }
    }

    async fn run_inner(
        &mut self,
        reporter: Option<Arc<dyn DownloadProgressReporter>>,
        cancel: &CancelSignal,
    ) -> Result<DownloadReport> {
        let result = self.execute_pending(reporter, cancel).await;
        if result.is_err() {
            self.discard();
        }
        result
    }

    async fn execute_pending(
        &mut self,
        reporter: Option<Arc<dyn DownloadProgressReporter>>,
        cancel: &CancelSignal,
    ) -> Result<DownloadReport> {
        self.prepare_staging()?;
        let pending: Vec<DownloadRecord> = self
            .records
            .iter()
            .filter(|r| !self.done.contains(&r.key()))
            .cloned()
            .collect();

        let known_total = known_bytes(&pending);
        let tracker = Arc::new(ProgressTracker::new(
            &self.records,
            &self.done,
            known_total,
            reporter,
        ));
        tracker.emit();

        let heartbeat = spawn_heartbeat(Arc::clone(&tracker));
        let mut reasons: HashMap<String, String> = HashMap::new();

        let seq_res = self
            .run_sequences(&pending, &mut reasons, &tracker, cancel)
            .await;
        let whole_res = if seq_res.is_ok() {
            self.run_whole_files(&pending, &mut reasons, &tracker, cancel)
                .await
        } else {
            Ok(())
        };

        heartbeat.abort();
        seq_res?;
        whole_res?;
        cancel.throw_if_cancelled()?;
        tracker.emit();

        let failures = self
            .records
            .iter()
            .filter(|r| !self.done.contains(&r.key()))
            .map(|r| DownloadFailure {
                record_key: r.key(),
                accession: r.accession.clone(),
                reason: reasons
                    .get(&r.key())
                    .cloned()
                    .unwrap_or_else(|| "Not downloaded.".to_string()),
            })
            .collect();

        Ok(DownloadReport {
            completed: self.done.len() as u32,
            failures,
        })
    }

    async fn run_sequences(
        &mut self,
        pending: &[DownloadRecord],
        reasons: &mut HashMap<String, String>,
        tracker: &Arc<ProgressTracker>,
        cancel: &CancelSignal,
    ) -> Result<()> {
        type BatchOutcome = (
            Vec<DownloadRecord>,
            PathBuf,
            PathBuf,
            Result<(String, u64, PreparedBatch)>,
        );

        let staging = self.staging.clone().expect("Staging prepared");
        let mut join_set: JoinSet<Result<(usize, usize, BatchOutcome)>> = JoinSet::new();

        for (source_idx, source) in [DatabaseSource::GenBank, DatabaseSource::Ena]
            .into_iter()
            .enumerate()
        {
            let records: Vec<DownloadRecord> = pending
                .iter()
                .filter(|r| r.source == source && r.content.is_sequence())
                .cloned()
                .collect();
            if records.is_empty() {
                continue;
            }
            let max_count = if source == DatabaseSource::GenBank {
                GENBANK_BATCH_COUNT
            } else {
                ENA_BATCH_COUNT
            };
            let host = match source {
                DatabaseSource::GenBank => endpoint_host(&self.endpoints.genbank_efetch),
                DatabaseSource::Ena => endpoint_host(&self.endpoints.ena_fasta_base),
            }
            .unwrap_or_else(|| "eutils.ncbi.nlm.nih.gov".to_string());
            let concurrency = self.client.limiter().max_limit(&host);
            let sem = Arc::new(Semaphore::new(concurrency));
            let batches = adaptive_batches(&records, max_count, MAX_BATCH_BYTES);

            for (batch_idx, batch) in batches.into_iter().enumerate() {
                let raw = self.next_temp_file(&staging);
                let parts_dir = self.next_temp_dir(&staging, "batch");
                let sem = Arc::clone(&sem);
                let client = Arc::clone(&self.client);
                let cancel = cancel.clone();
                let tracker = Arc::clone(tracker);
                let options = self.options.clone();
                let config = self.config.clone();
                let endpoints = self.endpoints.clone();
                let trims = self.trims_to_frame;

                join_set.spawn(async move {
                    let _permit = sem.acquire_owned().await?;
                    cancel.throw_if_cancelled()?;
                    for r in &batch {
                        tracker.set_status(&r.key(), RecordDownloadStatus::Active);
                    }
                    tracker.set_stage(DownloadStage::Downloading, None, None);

                    let stream_id = tracker.alloc_stream();
                    let res = fetch_and_prepare_batch(
                        batch_idx,
                        source,
                        batch.clone(),
                        raw.clone(),
                        parts_dir.clone(),
                        client,
                        options,
                        config,
                        trims,
                        endpoints,
                        stream_id,
                        Arc::clone(&tracker),
                        cancel.clone(),
                    )
                    .await;
                    tracker.finish_stream(stream_id);
                    Ok((source_idx, batch_idx, (batch, raw, parts_dir, res)))
                });
            }
        }

        let mut next_idx = [0usize; 2];
        let mut ready: [BTreeMap<usize, BatchOutcome>; 2] = Default::default();

        while let Some(joined) = join_set.join_next().await {
            let (source_idx, batch_idx, item) = joined??;
            ready[source_idx].insert(batch_idx, item);

            while let Some((batch, raw, parts_dir, outcome)) =
                ready[source_idx].remove(&next_idx[source_idx])
            {
                next_idx[source_idx] += 1;
                let cleanup = |raw_path: &Path, dir_path: &Path| {
                    let _ = fs::remove_file(raw_path);
                    let _ = fs::remove_dir_all(dir_path);
                };
                match outcome {
                    Ok((url, bytes_len, prepared)) => {
                        let apply_res = self.apply_prepared_batch(
                            &batch, &parts_dir, &url, bytes_len, prepared, reasons, tracker, cancel,
                        );
                        cleanup(&raw, &parts_dir);
                        if apply_res.is_err() {
                            join_set.abort_all();
                            for map in &ready {
                                for (_, rem_raw, rem_dir, _) in map.values() {
                                    cleanup(rem_raw, rem_dir);
                                }
                            }
                            apply_res?;
                        }
                    }
                    Err(err) => {
                        cleanup(&raw, &parts_dir);
                        if err.downcast_ref::<DownloadCancelled>().is_some() {
                            join_set.abort_all();
                            for map in &ready {
                                for (_, rem_raw, rem_dir, _) in map.values() {
                                    cleanup(rem_raw, rem_dir);
                                }
                            }
                            return Err(err);
                        }
                        let msg = err.to_string();
                        for r in &batch {
                            reasons.insert(r.key(), msg.clone());
                            tracker.set_status(&r.key(), RecordDownloadStatus::Failed);
                        }
                    }
                }
                tracker.emit();
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn apply_prepared_batch(
        &mut self,
        batch: &[DownloadRecord],
        parts_dir: &Path,
        url: &str,
        bytes_len: u64,
        prepared: PreparedBatch,
        reasons: &mut HashMap<String, String>,
        tracker: &Arc<ProgressTracker>,
        cancel: &CancelSignal,
    ) -> Result<()> {
        let staging = self.staging.clone().expect("Staging prepared");
        let is_single = self.options.effective_grouping() == SequenceGrouping::SingleFile;

        let part_map: HashMap<String, PathBuf> = batch
            .iter()
            .enumerate()
            .map(|(i, r)| {
                let name = if is_single {
                    "all".to_string()
                } else {
                    format!("{i}")
                };
                (r.key(), parts_dir.join(name))
            })
            .collect();

        // Re-unique IDs across concurrent batches and flush into staging in batch order.
        let mut renamed_ids: HashMap<String, String> = HashMap::new();
        let mut final_entries = Vec::with_capacity(prepared.sequences.len());
        for seq in prepared.sequences {
            cancel.throw_if_cancelled()?;
            let final_id = ensure_unique_id(&seq.id, &mut self.used_ids);
            if final_id != seq.id {
                renamed_ids.insert(seq.id.clone(), final_id.clone());
            }
            final_entries.push((seq, final_id));
        }

        let mut appended = HashSet::new();
        for (seq, final_id) in final_entries {
            cancel.throw_if_cancelled()?;
            let name = self.sequence_file_name(&seq.record);
            let target = staging.join(&name);
            let part = part_map.get(&seq.record.key()).expect("Part exists");

            if appended.insert(part.clone()) {
                if renamed_ids.is_empty() && !target.exists() && !is_single {
                    if fs::rename(part, &target).is_err() {
                        append_file(part, &target)?;
                    }
                } else if renamed_ids.is_empty() {
                    append_file(part, &target)?;
                } else {
                    append_file_renaming(part, &target, &renamed_ids)?;
                }
            }
            self.add_source(&seq.record, &name, url, &final_id, &seq.header, seq.trimmed);
            self.done.insert(seq.record.key());
            tracker.mark_done(&seq.record.key());
        }

        self.has_unknown_residues |= prepared.has_unknown_residues;
        tracker.add_completed_bytes(bytes_len);
        for record in batch {
            if !self.done.contains(&record.key()) {
                reasons.insert(record.key(), self.absent_reason(record));
                tracker.set_status(&record.key(), RecordDownloadStatus::Failed);
            }
        }
        Ok(())
    }

    async fn run_whole_files(
        &mut self,
        pending: &[DownloadRecord],
        reasons: &mut HashMap<String, String>,
        tracker: &Arc<ProgressTracker>,
        cancel: &CancelSignal,
    ) -> Result<()> {
        let whole: Vec<DownloadRecord> = pending
            .iter()
            .filter(|r| !r.content.is_sequence())
            .cloned()
            .collect();
        if whole.is_empty() {
            return Ok(());
        }

        let is_reads = whole[0].content == DownloadContent::Reads;
        if is_reads {
            self.run_reads(&whole, reasons, tracker, cancel).await
        } else {
            self.run_assemblies(&whole, reasons, tracker, cancel).await
        }
    }

    async fn run_reads(
        &mut self,
        records: &[DownloadRecord],
        reasons: &mut HashMap<String, String>,
        tracker: &Arc<ProgressTracker>,
        cancel: &CancelSignal,
    ) -> Result<()> {
        let staging = self.staging.clone().expect("Staging prepared");
        let read_host = records
            .iter()
            .flat_map(|r| r.files.iter())
            .find_map(|f| endpoint_host(&f.url))
            .unwrap_or_else(|| "ftp.sra.ebi.ac.uk".to_string());
        let host_limit = self.client.limiter().max_limit(&read_host);
        let total_files: usize = records.iter().map(|r| r.files.len()).sum();
        let concurrent_files = total_files.min(host_limit).max(1);
        let sem = Arc::new(Semaphore::new(host_limit));
        let mut run_handles = JoinSet::new();

        for record in records {
            if record.files.is_empty() {
                reasons.insert(
                    record.key(),
                    "ENA has no FASTQ files for this run.".to_string(),
                );
                tracker.set_status(&record.key(), RecordDownloadStatus::Failed);
                continue;
            }
            let record = record.clone();
            let mut file_jobs = Vec::with_capacity(record.files.len());
            for remote in &record.files {
                file_jobs.push((remote.clone(), self.next_temp_file(&staging)));
            }
            let sem = Arc::clone(&sem);
            let client = Arc::clone(&self.client);
            let tracker = Arc::clone(tracker);
            let cancel = cancel.clone();

            run_handles.spawn(async move {
                let mut f_handles = Vec::with_capacity(file_jobs.len());
                for (remote, part) in &file_jobs {
                    let sem = Arc::clone(&sem);
                    let client = Arc::clone(&client);
                    let tracker = Arc::clone(&tracker);
                    let cancel = cancel.clone();
                    let remote = remote.clone();
                    let part = part.clone();
                    let key = record.key();
                    f_handles.push(tokio::spawn(async move {
                        let _permit = sem.acquire_owned().await?;
                        tracker.set_status(&key, RecordDownloadStatus::Active);
                        download_and_verify_read(
                            &client,
                            &remote,
                            &part,
                            concurrent_files,
                            &tracker,
                            &cancel,
                        )
                        .await
                    }));
                }
                let mut saved = Vec::with_capacity(file_jobs.len());
                let mut first_err: Option<anyhow::Error> = None;
                for (handle, (remote, part)) in f_handles.into_iter().zip(file_jobs) {
                    match handle.await? {
                        Ok(size) => saved.push((part, remote, size)),
                        Err(err) => {
                            let _ = fs::remove_file(&part);
                            if first_err.is_none() {
                                first_err = Some(err);
                            }
                        }
                    }
                }
                if let Some(err) = first_err {
                    for (part, _, _) in saved {
                        let _ = fs::remove_file(part);
                    }
                    return Ok::<_, anyhow::Error>((record, Err(err)));
                }
                Ok((record, Ok(saved)))
            });
        }

        while let Some(joined) = run_handles.join_next().await {
            let (record, outcome) = joined??;
            match outcome {
                Ok(saved) => {
                    for (part, remote, size) in saved {
                        let remote_name = Url::parse(&remote.url)
                            .ok()
                            .and_then(|u| {
                                u.path_segments()
                                    .and_then(|mut s| s.next_back().map(str::to_string))
                            })
                            .unwrap_or_else(|| record.accession.clone());
                        let stem = safe_name(&FASTQ_EXT_RE.replace(&remote_name, ""));
                        let name = self.unique_name(&stem, ".fastq.gz", "-");
                        fs::rename(&part, staging.join(&name))?;
                        self.file_names.push(name.clone());
                        self.add_source(&record, &name, &remote.url, "", "", 0);
                        tracker.add_completed_bytes(size);
                    }
                    self.done.insert(record.key());
                    tracker.mark_done(&record.key());
                }
                Err(err) => {
                    if err.downcast_ref::<DownloadCancelled>().is_some() {
                        run_handles.abort_all();
                        return Err(err);
                    }
                    reasons.insert(record.key(), err.to_string());
                    tracker.set_status(&record.key(), RecordDownloadStatus::Failed);
                }
            }
            tracker.emit();
        }
        Ok(())
    }

    async fn run_assemblies(
        &mut self,
        records: &[DownloadRecord],
        reasons: &mut HashMap<String, String>,
        tracker: &Arc<ProgressTracker>,
        cancel: &CancelSignal,
    ) -> Result<()> {
        let staging = self.staging.clone().expect("Staging prepared");
        let gb_host = records
            .iter()
            .filter(|r| r.source == DatabaseSource::GenBank)
            .find_map(|r| r.files.first().and_then(|f| endpoint_host(&f.url)))
            .unwrap_or_else(|| "ftp.ncbi.nlm.nih.gov".to_string());
        let ena_host = endpoint_host(&self.endpoints.ena_fasta_base)
            .unwrap_or_else(|| "www.ebi.ac.uk".to_string());

        let gb_limit = self.client.limiter().max_limit(&gb_host);
        let ena_limit = self.client.limiter().max_limit(&ena_host);
        let gb_count = records
            .iter()
            .filter(|r| r.source == DatabaseSource::GenBank)
            .count();
        let gb_concurrent = gb_count.min(gb_limit).max(1);

        let gb_sem = Arc::new(Semaphore::new(gb_limit));
        let ena_sem = Arc::new(Semaphore::new(ena_limit));
        let mut handles = JoinSet::new();

        for record in records {
            let record = record.clone();
            let raw = self.next_temp_file(&staging);
            let ready = self.next_temp_file(&staging);
            let sem = if record.source == DatabaseSource::GenBank {
                Arc::clone(&gb_sem)
            } else {
                Arc::clone(&ena_sem)
            };
            let client = Arc::clone(&self.client);
            let tracker = Arc::clone(tracker);
            let cancel = cancel.clone();
            let endpoints = self.endpoints.clone();

            handles.spawn(async move {
                let _permit = sem.acquire_owned().await?;
                cancel.throw_if_cancelled()?;
                tracker.set_status(&record.key(), RecordDownloadStatus::Active);
                let res = fetch_and_prepare_assembly(
                    &record,
                    &raw,
                    &ready,
                    gb_concurrent,
                    &client,
                    &endpoints,
                    &tracker,
                    &cancel,
                )
                .await;
                let _ = fs::remove_file(&raw);
                if res.is_err() {
                    let _ = fs::remove_file(&ready);
                }
                Ok::<_, anyhow::Error>((record, ready, res))
            });
        }

        while let Some(joined) = handles.join_next().await {
            let (record, ready, outcome) = joined??;
            match outcome {
                Ok((url, raw_size)) => {
                    let stem = stem_of(&record);
                    let name = self.unique_name(&stem, ".fna", "_");
                    fs::rename(&ready, staging.join(&name))?;
                    self.file_names.push(name.clone());
                    self.add_source(&record, &name, &url, "", "", 0);
                    self.done.insert(record.key());
                    tracker.add_completed_bytes(raw_size);
                    tracker.mark_done(&record.key());
                }
                Err(err) => {
                    if err.downcast_ref::<DownloadCancelled>().is_some() {
                        handles.abort_all();
                        return Err(err);
                    }
                    reasons.insert(record.key(), err.to_string());
                    tracker.set_status(&record.key(), RecordDownloadStatus::Failed);
                }
            }
            tracker.emit();
        }
        Ok(())
    }

    fn prepare_staging(&mut self) -> Result<()> {
        if self.staging.is_some() {
            return Ok(());
        }
        fs::create_dir_all(&self.root)?;
        let pid = std::process::id();
        let mut attempt = 0u64;
        loop {
            attempt += 1;
            let name = format!(".partial-{pid}-{attempt}");
            let candidate = self.root.join(name);
            if !candidate.exists() {
                fs::create_dir(&candidate)?;
                self.staging = Some(candidate);
                return Ok(());
            }
        }
    }

    fn next_temp_file(&mut self, staging: &Path) -> PathBuf {
        self.temporary_seq += 1;
        staging.join(format!("download-{}.part", self.temporary_seq))
    }

    fn next_temp_dir(&mut self, staging: &Path, prefix: &str) -> PathBuf {
        self.temporary_seq += 1;
        let dir = staging.join(format!("{prefix}-{}", self.temporary_seq));
        let _ = fs::create_dir_all(&dir);
        dir
    }

    fn sequence_file_name(&mut self, record: &DownloadRecord) -> String {
        if self.options.effective_grouping() == SequenceGrouping::PerRecord {
            let key = record.key();
            if let Some(existing) = self.name_of.get(&key) {
                return existing.clone();
            }
            let name = self.unique_name(&stem_of(record), ".fasta", "_");
            self.file_names.push(name.clone());
            self.name_of.insert(key, name.clone());
            return name;
        }
        if let Some(name) = self.single_name.as_ref() {
            return name.clone();
        }
        let mut stem = self.options.file_name.trim();
        if stem.is_empty() {
            stem = self.options.gene.trim();
        }
        if stem.is_empty() {
            stem = "sequences";
        }
        let clean = safe_name(stem).replace('.', "_");
        let name = self.unique_name(&clean, ".fasta", "_");
        self.file_names.push(name.clone());
        self.single_name = Some(name.clone());
        name
    }

    fn unique_name(&mut self, stem: &str, extension: &str, separator: &str) -> String {
        let mut candidate = stem.to_string();
        let mut idx = 2u32;
        while self.taken_names.contains(&candidate.to_lowercase())
            || self.used_names.contains(&candidate.to_lowercase())
        {
            candidate = format!("{stem}{separator}{idx}");
            idx += 1;
        }
        self.used_names.insert(candidate.to_lowercase());
        format!("{candidate}{extension}")
    }

    fn add_source(
        &mut self,
        record: &DownloadRecord,
        file: &str,
        url: &str,
        id: &str,
        header: &str,
        trimmed: u32,
    ) {
        let stamp = Utc::now().to_rfc3339_opts(chrono::SecondsFormat::Micros, true);
        self.sources.push([
            record.source.label().to_string(),
            record.accession.clone(),
            file.to_string(),
            id.to_string(),
            header.to_string(),
            trimmed.to_string(),
            url.to_string(),
            stamp,
        ]);
    }

    fn absent_reason(&self, record: &DownloadRecord) -> String {
        if record.content != DownloadContent::Coding {
            return format!("{} did not return this record.", record.source.label());
        }
        let gene = self.options.gene.trim();
        if gene.is_empty() {
            "No coding sequence in this record.".to_string()
        } else {
            format!("No {gene} coding sequence in this record.")
        }
    }
}

#[allow(clippy::too_many_arguments)]
async fn fetch_and_prepare_batch(
    _idx: usize,
    source: DatabaseSource,
    batch: Vec<DownloadRecord>,
    raw: PathBuf,
    parts_dir: PathBuf,
    client: Arc<ArchiveClient>,
    options: DownloadOptions,
    config: DownloadConfig,
    trims_to_frame: bool,
    endpoints: Endpoints,
    stream_id: u64,
    tracker: Arc<ProgressTracker>,
    cancel: CancelSignal,
) -> Result<(String, u64, PreparedBatch)> {
    let is_single = options.effective_grouping() == SequenceGrouping::SingleFile;
    let seq_req = build_sequence_request(source, &batch, &config, &endpoints)?;

    let tracker_bytes = Arc::clone(&tracker);
    let tracker_retry = Arc::clone(&tracker);
    let save_res = client
        .save_file(
            &seq_req.url,
            &raw,
            &SaveOptions {
                form: seq_req.form,
                allow_resume: false,
                allow_parallel_ranges: false,
                known_size: None,
                compute_md5: false,
                concurrent_files: 1,
            },
            &cancel,
            &move |received, total| {
                tracker_bytes.update_stream(stream_id, received, total);
            },
            &move |attempt, max| {
                tracker_retry.set_stage(DownloadStage::Retrying, Some(attempt), Some(max));
            },
        )
        .await?;

    tracker.set_stage(DownloadStage::Preparing, None, None);
    let part_paths: HashMap<String, PathBuf> = batch
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let name = if is_single {
                "all".to_string()
            } else {
                format!("{i}")
            };
            (r.key(), parts_dir.join(name))
        })
        .collect();

    let raw_clone = raw.clone();
    let cancel_clone = cancel.clone();
    let prepared = tokio::task::spawn_blocking(move || {
        prepare_batch(
            &raw_clone,
            source,
            &batch,
            &options,
            trims_to_frame,
            &part_paths,
            &cancel_clone,
        )
    })
    .await??;

    Ok((seq_req.recorded_url, save_res.bytes_written, prepared))
}

async fn download_and_verify_read(
    client: &ArchiveClient,
    remote: &RemoteFile,
    part: &Path,
    concurrent_files: usize,
    tracker: &Arc<ProgressTracker>,
    cancel: &CancelSignal,
) -> Result<u64> {
    let url = Url::parse(&remote.url)?;
    let stream_id = tracker.alloc_stream();
    tracker.set_stage(DownloadStage::Downloading, None, None);

    let tracker_bytes = Arc::clone(tracker);
    let tracker_retry = Arc::clone(tracker);
    let res = client
        .save_file(
            &url,
            part,
            &SaveOptions {
                form: None,
                allow_resume: true,
                allow_parallel_ranges: true,
                known_size: remote.bytes,
                compute_md5: remote.md5.is_some(),
                concurrent_files,
            },
            cancel,
            &move |received, total| {
                tracker_bytes.update_stream(stream_id, received, total);
            },
            &move |attempt, max| {
                tracker_retry.set_stage(DownloadStage::Retrying, Some(attempt), Some(max));
            },
        )
        .await;
    tracker.finish_stream(stream_id);
    let save_res = res?;

    tracker.set_stage(DownloadStage::Verifying, None, None);
    let size = fs::metadata(part)?.len();
    verify_remote(remote, save_res.md5.as_deref(), size)?;
    Ok(size)
}

#[allow(clippy::too_many_arguments)]
async fn fetch_and_prepare_assembly(
    record: &DownloadRecord,
    raw: &Path,
    ready: &Path,
    gb_concurrent_files: usize,
    client: &ArchiveClient,
    endpoints: &Endpoints,
    tracker: &Arc<ProgressTracker>,
    cancel: &CancelSignal,
) -> Result<(String, u64)> {
    let stream_id = tracker.alloc_stream();
    tracker.set_stage(DownloadStage::Downloading, None, None);

    if record.source == DatabaseSource::GenBank {
        let Some(first) = record.files.first() else {
            tracker.finish_stream(stream_id);
            return Err(DownloadError("NCBI has no FASTA for this assembly.".to_string()).into());
        };
        let url = Url::parse(&first.url)?;
        let checksum_fut = fetch_ncbi_checksum(client, &url, cancel);
        let tracker_bytes = Arc::clone(tracker);
        let tracker_retry = Arc::clone(tracker);
        let save_opts = SaveOptions {
            form: None,
            allow_resume: true,
            allow_parallel_ranges: true,
            known_size: first.bytes,
            compute_md5: true,
            concurrent_files: gb_concurrent_files,
        };
        let on_bytes = move |received, total| {
            tracker_bytes.update_stream(stream_id, received, total);
        };
        let on_retry = move |attempt, max| {
            tracker_retry.set_stage(DownloadStage::Retrying, Some(attempt), Some(max));
        };
        let save_fut = client.save_file(&url, raw, &save_opts, cancel, &on_bytes, &on_retry);
        let (expected_res, save_res) = tokio::join!(checksum_fut, save_fut);
        tracker.finish_stream(stream_id);

        let expected = expected_res?;
        let save = save_res?;
        tracker.set_stage(DownloadStage::Verifying, None, None);
        let raw_size = fs::metadata(raw)?.len();
        let remote = RemoteFile {
            url: first.url.clone(),
            bytes: first.bytes,
            md5: expected,
        };
        verify_remote(&remote, save.md5.as_deref(), raw_size)?;

        tracker.set_stage(DownloadStage::Decompressing, None, None);
        let raw_path = raw.to_path_buf();
        let ready_path = ready.to_path_buf();
        let cancel_clone = cancel.clone();
        tokio::task::spawn_blocking(move || {
            decompress_fasta(&raw_path, &ready_path, &cancel_clone)
        })
        .await??;
        return Ok((first.url.clone(), raw_size));
    }

    let Some(expected_bases) = record.length else {
        tracker.finish_stream(stream_id);
        return Err(DownloadError(
            "ENA gives no size to check this assembly. Try GenBank.".to_string(),
        )
        .into());
    };
    let report_url = ena_assembly_report_url(&endpoints.ena_report_base, &record.accession)?;
    let fasta_url = ena_assembly_fasta_url(&endpoints.ena_fasta_base, &record.accession)?;

    let report_fut = fetch_ena_report(client, &report_url, cancel);
    let tracker_bytes = Arc::clone(tracker);
    let tracker_retry = Arc::clone(tracker);
    let save_opts = SaveOptions {
        form: None,
        allow_resume: false,
        allow_parallel_ranges: false,
        known_size: None,
        compute_md5: false,
        concurrent_files: 1,
    };
    let on_bytes = move |received, total| {
        tracker_bytes.update_stream(stream_id, received, total);
    };
    let on_retry = move |attempt, max| {
        tracker_retry.set_stage(DownloadStage::Retrying, Some(attempt), Some(max));
    };
    let save_fut = client.save_file(&fasta_url, raw, &save_opts, cancel, &on_bytes, &on_retry);
    let (report_res, save_res) = tokio::join!(report_fut, save_fut);
    tracker.finish_stream(stream_id);

    let members = report_res?;
    let _save = save_res?;
    let raw_size = fs::metadata(raw)?.len();
    tracker.set_stage(DownloadStage::Preparing, None, None);
    let raw_path = raw.to_path_buf();
    let ready_path = ready.to_path_buf();
    let cancel_clone = cancel.clone();
    tokio::task::spawn_blocking(move || {
        write_ena_assembly(
            &raw_path,
            &ready_path,
            members.as_ref(),
            expected_bases,
            &cancel_clone,
        )
    })
    .await??;

    Ok((fasta_url.to_string(), raw_size))
}

async fn fetch_ncbi_checksum(
    client: &ArchiveClient,
    file_url: &Url,
    cancel: &CancelSignal,
) -> Result<Option<String>> {
    let mut list_url = file_url.clone();
    let dir = Path::new(file_url.path())
        .parent()
        .map(|p| p.to_string_lossy().replace('\\', "/"))
        .unwrap_or_else(|| "/".to_string());
    let joined = if dir.ends_with('/') {
        format!("{dir}md5checksums.txt")
    } else {
        format!("{dir}/md5checksums.txt")
    };
    list_url.set_path(&joined);

    let (status, body) = client
        .read_text(&list_url, None, cancel, &|_, _| {})
        .await?;
    if status == 404 {
        return Ok(None);
    }
    if status != 200 {
        return Err(DownloadError("Could not read the checksum. Try again.".to_string()).into());
    }
    let file_name = file_url
        .path_segments()
        .and_then(|mut s| s.next_back())
        .unwrap_or("");
    let target = format!("./{file_name}");
    for line in body.lines() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() == 2 && parts[1] == target {
            return Ok(Some(parts[0].to_string()));
        }
    }
    Ok(None)
}

async fn fetch_ena_report(
    client: &ArchiveClient,
    report_url: &Url,
    cancel: &CancelSignal,
) -> Result<Option<HashSet<String>>> {
    let (status, body) = client
        .read_text(report_url, None, cancel, &|_, _| {})
        .await?;
    if status == 404 {
        return Ok(None);
    }
    if status != 200 {
        return Err(
            DownloadError("Could not read the assembly report. Try again.".to_string()).into(),
        );
    }
    Ok(Some(parse_sequence_report(&body)?))
}

struct SequenceRequest {
    url: Url,
    form: Option<Vec<(String, String)>>,
    recorded_url: String,
}

fn build_sequence_request(
    source: DatabaseSource,
    batch: &[DownloadRecord],
    config: &DownloadConfig,
    endpoints: &Endpoints,
) -> Result<SequenceRequest> {
    match source {
        DatabaseSource::GenBank => {
            let content = batch[0].content;
            let db = match content {
                DownloadContent::Nucleotide | DownloadContent::Coding => "nuccore",
                DownloadContent::Protein => "protein",
                _ => {
                    return Err(
                        DownloadError("GenBank has no FASTQ files. Use ENA.".to_string()).into(),
                    );
                }
            };
            let rettype = if content == DownloadContent::Coding {
                "fasta_cds_na"
            } else {
                "fasta"
            };
            let ids: String = batch
                .iter()
                .map(|r| r.accession.as_str())
                .collect::<Vec<_>>()
                .join(",");
            let base = Url::parse(&endpoints.genbank_efetch)?;
            let mut recorded = base.clone();
            recorded
                .query_pairs_mut()
                .append_pair("db", db)
                .append_pair("rettype", rettype)
                .append_pair("retmode", "text")
                .append_pair("id", &ids);

            let mut form = vec![
                ("db".to_string(), db.to_string()),
                ("rettype".to_string(), rettype.to_string()),
                ("retmode".to_string(), "text".to_string()),
                ("id".to_string(), ids),
                ("tool".to_string(), config.tool_name.clone()),
                ("email".to_string(), config.contact_email.clone()),
            ];
            if let Some(key) = config.ncbi_api_key.as_deref().map(str::trim)
                && !key.is_empty()
            {
                form.push(("api_key".to_string(), key.to_string()));
            }
            Ok(SequenceRequest {
                url: base,
                form: Some(form),
                recorded_url: recorded.to_string(),
            })
        }
        DatabaseSource::Ena => {
            let ids: String = batch
                .iter()
                .map(|r| base_accession(&r.accession))
                .collect::<Vec<_>>()
                .join(",");
            let mut url = Url::parse(&format!(
                "{}/{}",
                endpoints.ena_fasta_base.trim_end_matches('/'),
                ids
            ))?;
            let recorded = url.to_string();
            url.query_pairs_mut()
                .append_pair("download", "true")
                .append_pair("gzip", "true");
            Ok(SequenceRequest {
                url,
                form: None,
                recorded_url: recorded,
            })
        }
    }
}

fn ena_assembly_fasta_url(base: &str, accession: &str) -> Result<Url> {
    let mut url = Url::parse(&format!("{}/{}", base.trim_end_matches('/'), accession))?;
    url.query_pairs_mut()
        .append_pair("download", "true")
        .append_pair("gzip", "true");
    Ok(url)
}

fn ena_assembly_report_url(base: &str, accession: &str) -> Result<Url> {
    let root = base_accession(accession);
    let p1 = &root[..root.len().min(7)];
    let p2 = &root[..root.len().min(10)];
    Ok(Url::parse(&format!(
        "{}/{p1}/{p2}/{accession}_sequence_report.txt",
        base.trim_end_matches('/')
    ))?)
}

fn verify_remote(remote: &RemoteFile, md5: Option<&str>, size: u64) -> Result<()> {
    if let Some(expected_bytes) = remote.bytes
        && expected_bytes != size
    {
        return Err(
            DownloadError("The file size does not match. Retry the download.".to_string()).into(),
        );
    }
    if let Some(expected_md5) = remote.md5.as_deref() {
        let actual = md5.unwrap_or_default();
        if !expected_md5.eq_ignore_ascii_case(actual) {
            return Err(DownloadError(
                "The checksum does not match. Retry the download.".to_string(),
            )
            .into());
        }
    }
    Ok(())
}

fn append_file(source: &Path, dest: &Path) -> io::Result<()> {
    let mut reader = io::BufReader::new(File::open(source)?);
    let mut writer = BufWriter::new(OpenOptions::new().create(true).append(true).open(dest)?);
    io::copy(&mut reader, &mut writer)?;
    writer.flush()
}

fn append_file_renaming(
    source: &Path,
    dest: &Path,
    renamed: &HashMap<String, String>,
) -> io::Result<()> {
    use io::BufRead;
    let reader = io::BufReader::new(File::open(source)?);
    let mut writer = BufWriter::new(OpenOptions::new().create(true).append(true).open(dest)?);
    for line in reader.lines() {
        let line = line?;
        if let Some(id) = line.strip_prefix('>')
            && let Some(new_id) = renamed.get(id.trim())
        {
            writeln!(writer, ">{new_id}")?;
            continue;
        }
        writeln!(writer, "{line}")?;
    }
    writer.flush()
}

fn ensure_unique_id(base: &str, used: &mut HashSet<String>) -> String {
    if used.insert(base.to_string()) {
        return base.to_string();
    }
    for i in 2.. {
        let candidate = format!("{base}_{i}");
        if used.insert(candidate.clone()) {
            return candidate;
        }
    }
    unreachable!()
}

fn stem_of(record: &DownloadRecord) -> String {
    safe_name(&record.accession).replace('.', "_")
}

fn endpoint_host(url: &str) -> Option<String> {
    Url::parse(url)
        .ok()
        .and_then(|u| u.host_str().map(str::to_string))
}

fn known_bytes(pending: &[DownloadRecord]) -> u64 {
    if pending.iter().any(|r| r.content != DownloadContent::Reads) {
        return 0;
    }
    pending.iter().filter_map(|r| r.estimated_bytes()).sum()
}

fn spawn_cancel_watcher(
    reporter: Option<Arc<dyn DownloadProgressReporter>>,
    cancel: &CancelSignal,
) -> Option<tokio::task::JoinHandle<()>> {
    let reporter = reporter?;
    let cancel = cancel.clone();
    Some(tokio::spawn(async move {
        loop {
            if reporter.is_cancelled() {
                cancel.cancel();
                break;
            }
            tokio::time::sleep(Duration::from_millis(100)).await;
        }
    }))
}

struct TrackerInner {
    total: u32,
    done: u32,
    completed_bytes: u64,
    static_bytes_total: u64,
    in_flight: HashMap<u64, (u64, Option<u64>)>,
    stage: DownloadStage,
    retry_attempt: Option<u32>,
    retry_max: Option<u32>,
    statuses: HashMap<String, RecordDownloadStatus>,
    order: Vec<String>,
    last_emit: Instant,
    last_sample_at: Instant,
    last_sample_bytes: u64,
    ema_bps: Option<f64>,
}

struct ProgressTracker {
    inner: Mutex<TrackerInner>,
    next_stream_id: AtomicU64,
    reporter: Option<Arc<dyn DownloadProgressReporter>>,
}

impl ProgressTracker {
    fn new(
        records: &[DownloadRecord],
        done_keys: &HashSet<String>,
        static_bytes_total: u64,
        reporter: Option<Arc<dyn DownloadProgressReporter>>,
    ) -> Self {
        let mut statuses = HashMap::new();
        let mut order = Vec::with_capacity(records.len());
        for record in records {
            let key = record.key();
            let status = if done_keys.contains(&key) {
                RecordDownloadStatus::Done
            } else {
                RecordDownloadStatus::Pending
            };
            statuses.insert(key.clone(), status);
            order.push(key);
        }
        let now = Instant::now();
        Self {
            inner: Mutex::new(TrackerInner {
                total: records.len() as u32,
                done: done_keys.len() as u32,
                completed_bytes: 0,
                static_bytes_total,
                in_flight: HashMap::new(),
                stage: DownloadStage::Connecting,
                retry_attempt: None,
                retry_max: None,
                statuses,
                order,
                last_emit: now - Duration::from_secs(1),
                last_sample_at: now,
                last_sample_bytes: 0,
                ema_bps: None,
            }),
            next_stream_id: AtomicU64::new(1),
            reporter,
        }
    }

    fn alloc_stream(&self) -> u64 {
        self.next_stream_id.fetch_add(1, Ordering::Relaxed)
    }

    fn update_stream(&self, id: u64, received: u64, total: Option<u64>) {
        let snapshot = {
            let mut inner = self.inner.lock().expect("Tracker lock");
            inner.in_flight.insert(id, (received, total));
            if inner.last_emit.elapsed() < Duration::from_millis(120) {
                return;
            }
            Self::build_snapshot(&mut inner)
        };
        if let Some(reporter) = self.reporter.as_ref() {
            reporter.on_progress(snapshot);
        }
    }

    fn finish_stream(&self, id: u64) {
        let mut inner = self.inner.lock().expect("Tracker lock");
        inner.in_flight.remove(&id);
    }

    fn add_completed_bytes(&self, bytes: u64) {
        let mut inner = self.inner.lock().expect("Tracker lock");
        inner.completed_bytes = inner.completed_bytes.saturating_add(bytes);
    }

    fn set_stage(&self, stage: DownloadStage, attempt: Option<u32>, max: Option<u32>) {
        let snapshot = {
            let mut inner = self.inner.lock().expect("Tracker lock");
            if inner.stage == stage && inner.retry_attempt == attempt && inner.retry_max == max {
                return;
            }
            inner.stage = stage;
            inner.retry_attempt = attempt;
            inner.retry_max = max;
            Self::build_snapshot(&mut inner)
        };
        if let Some(reporter) = self.reporter.as_ref() {
            reporter.on_progress(snapshot);
        }
    }

    fn set_status(&self, key: &str, status: RecordDownloadStatus) {
        let snapshot = {
            let mut inner = self.inner.lock().expect("Tracker lock");
            let changed = inner.statuses.get(key) != Some(&status);
            inner.statuses.insert(key.to_string(), status);
            if !changed {
                return;
            }
            Self::build_snapshot(&mut inner)
        };
        if let Some(reporter) = self.reporter.as_ref() {
            reporter.on_progress(snapshot);
        }
    }

    fn mark_done(&self, key: &str) {
        let snapshot = {
            let mut inner = self.inner.lock().expect("Tracker lock");
            if inner.statuses.get(key) != Some(&RecordDownloadStatus::Done) {
                inner.done += 1;
            }
            inner
                .statuses
                .insert(key.to_string(), RecordDownloadStatus::Done);
            Self::build_snapshot(&mut inner)
        };
        if let Some(reporter) = self.reporter.as_ref() {
            reporter.on_progress(snapshot);
        }
    }

    fn emit(&self) {
        let snapshot = {
            let mut inner = self.inner.lock().expect("Tracker lock");
            Self::build_snapshot(&mut inner)
        };
        if let Some(reporter) = self.reporter.as_ref() {
            reporter.on_progress(snapshot);
        }
    }

    fn build_snapshot(inner: &mut TrackerInner) -> DownloadProgress {
        let now = Instant::now();
        inner.last_emit = now;
        let in_flight_sum: u64 = inner.in_flight.values().map(|(b, _)| *b).sum();
        let bytes = inner.completed_bytes.saturating_add(in_flight_sum);

        let bytes_total = if inner.static_bytes_total > 0 {
            inner.static_bytes_total
        } else if !inner.in_flight.is_empty() && inner.in_flight.values().all(|(_, t)| t.is_some())
        {
            inner
                .completed_bytes
                .saturating_add(inner.in_flight.values().filter_map(|(_, t)| *t).sum())
        } else {
            0
        };

        let dt = now.duration_since(inner.last_sample_at).as_secs_f64();
        if dt >= 0.5 {
            let delta = bytes.saturating_sub(inner.last_sample_bytes) as f64;
            let inst = delta / dt;
            inner.ema_bps = Some(match inner.ema_bps {
                Some(prev) => 0.3 * inst + 0.7 * prev,
                None => inst,
            });
            inner.last_sample_at = now;
            inner.last_sample_bytes = bytes;
        }

        let statuses = inner
            .order
            .iter()
            .filter_map(|k| inner.statuses.get(k).map(|&s| (k.clone(), s)))
            .collect();

        DownloadProgress {
            done: inner.done,
            total: inner.total,
            bytes,
            bytes_total,
            bytes_per_second: inner.ema_bps.map(|v| v.max(0.0) as u64),
            stage: inner.stage,
            retry_attempt: inner.retry_attempt,
            retry_max: inner.retry_max,
            statuses,
        }
    }
}

fn spawn_heartbeat(tracker: Arc<ProgressTracker>) -> tokio::task::JoinHandle<()> {
    tokio::spawn(async move {
        loop {
            tokio::time::sleep(Duration::from_millis(500)).await;
            tracker.emit();
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use md5::{Digest, Md5};
    use tempfile::TempDir;
    use tokio::io::{AsyncReadExt, AsyncWriteExt};
    use tokio::net::TcpListener;

    fn loopback_client(has_key: bool) -> Arc<ArchiveClient> {
        Arc::new(
            ArchiveClient::with_validator(
                has_key,
                Arc::new(|u| u.host_str() == Some("127.0.0.1")),
                Duration::from_millis(5),
            )
            .unwrap(),
        )
    }

    #[tokio::test]
    async fn test_genbank_api_key_sent_in_post_but_omitted_from_sources_tsv() {
        let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
        let addr = listener.local_addr().unwrap();
        let seen_body = Arc::new(Mutex::new(String::new()));
        let seen_body_srv = Arc::clone(&seen_body);

        tokio::spawn(async move {
            loop {
                let Ok((mut sock, _)) = listener.accept().await else {
                    break;
                };
                let seen = Arc::clone(&seen_body_srv);
                tokio::spawn(async move {
                    let mut buf = [0u8; 4096];
                    let n = sock.read(&mut buf).await.unwrap_or(0);
                    let req = String::from_utf8_lossy(&buf[..n]).into_owned();
                    *seen.lock().unwrap() = req;
                    let fasta = ">U49845.1 Yeast\nACGT\n";
                    let resp = format!(
                        "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{}",
                        fasta.len(),
                        fasta
                    );
                    let _ = sock.write_all(resp.as_bytes()).await;
                });
            }
        });

        let dir = TempDir::new().unwrap();
        let endpoints = Endpoints {
            genbank_efetch: format!("http://127.0.0.1:{}/efetch.fcgi", addr.port()),
            ..Default::default()
        };
        let record = DownloadRecord {
            source: DatabaseSource::GenBank,
            content: DownloadContent::Nucleotide,
            accession: "U49845.1".to_string(),
            organism: "Mus musculus".to_string(),
            title: String::new(),
            length: Some(4),
            codon_start: 1,
            gene: String::new(),
            files: Vec::new(),
        };
        let mut job = DownloadJob::with_client(
            loopback_client(true),
            vec![record],
            DownloadOptions::default(),
            DownloadConfig {
                ncbi_api_key: Some("secret_key_123".to_string()),
                ..Default::default()
            },
            false,
            dir.path().to_path_buf(),
            HashSet::new(),
            endpoints,
        );
        let report = job.run(None).await.unwrap();
        assert!(report.is_complete());
        let completed = job.commit().unwrap();
        let sources =
            fs::read_to_string(completed.files[0].parent().unwrap().join("sources.tsv")).unwrap();
        assert!(seen_body.lock().unwrap().contains("api_key=secret_key_123"));
        assert!(!sources.contains("secret_key_123"));
    }

    #[tokio::test]
    async fn test_concurrent_paired_reads_and_ena_gzip_assembly() {
        let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
        let addr = listener.local_addr().unwrap();

        let fastq1 = b"@r1\nAC\n+\nII\n".to_vec();
        let fastq2 = b"@r2\nGT\n+\nII\n".to_vec();
        let md5_1 = format!("{:x}", Md5::digest(&fastq1));
        let md5_2 = format!("{:x}", Md5::digest(&fastq2));

        let mut gz_buf = Vec::new();
        {
            let mut enc = GzEncoder::new(&mut gz_buf, Compression::default());
            enc.write_all(
                b">ENA|CM000001|CM000001.1 chr\nACGTACGT\n>ENA|JACO02|JACO02.1 contig\nGGCC\n",
            )
            .unwrap();
            enc.finish().unwrap();
        }
        let report_txt = "accession\tname\nCM000001.1\tI\n";

        let f1 = fastq1.clone();
        let f2 = fastq2.clone();
        let gz = gz_buf.clone();

        tokio::spawn(async move {
            loop {
                let Ok((mut sock, _)) = listener.accept().await else {
                    break;
                };
                let f1 = f1.clone();
                let f2 = f2.clone();
                let gz = gz.clone();
                tokio::spawn(async move {
                    let mut buf = [0u8; 2048];
                    let n = sock.read(&mut buf).await.unwrap_or(0);
                    let req = String::from_utf8_lossy(&buf[..n]);
                    let path = req.lines().next().unwrap_or("");
                    let body: &[u8] = if path.contains("_1.fastq.gz") {
                        &f1
                    } else if path.contains("_2.fastq.gz") {
                        &f2
                    } else if path.contains("sequence_report.txt") {
                        report_txt.as_bytes()
                    } else {
                        &gz
                    };
                    let hdr = format!(
                        "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                        body.len()
                    );
                    let _ = sock.write_all(hdr.as_bytes()).await;
                    let _ = sock.write_all(body).await;
                });
            }
        });

        let dir = TempDir::new().unwrap();
        let base = format!("http://127.0.0.1:{}", addr.port());
        let read_rec = DownloadRecord {
            source: DatabaseSource::Ena,
            content: DownloadContent::Reads,
            accession: "SRR1".to_string(),
            organism: "Mus".to_string(),
            title: String::new(),
            length: None,
            codon_start: 1,
            gene: String::new(),
            files: vec![
                RemoteFile {
                    url: format!("{base}/SRR1_1.fastq.gz"),
                    bytes: Some(fastq1.len() as u64),
                    md5: Some(md5_1),
                },
                RemoteFile {
                    url: format!("{base}/SRR1_2.fastq.gz"),
                    bytes: Some(fastq2.len() as u64),
                    md5: Some(md5_2),
                },
            ],
        };
        let mut job = DownloadJob::with_client(
            loopback_client(false),
            vec![read_rec],
            DownloadOptions::default(),
            DownloadConfig::default(),
            false,
            dir.path().to_path_buf(),
            HashSet::new(),
            Endpoints::default(),
        );
        assert!(job.run(None).await.unwrap().is_complete());
        let completed = job.commit().unwrap();
        assert_eq!(completed.files.len(), 2);

        let asm_rec = DownloadRecord {
            source: DatabaseSource::Ena,
            content: DownloadContent::Assembly,
            accession: "GCA_000733235.3".to_string(),
            organism: "Mus".to_string(),
            title: String::new(),
            length: Some(8),
            codon_start: 1,
            gene: String::new(),
            files: Vec::new(),
        };
        let mut asm_job = DownloadJob::with_client(
            loopback_client(false),
            vec![asm_rec],
            DownloadOptions::default(),
            DownloadConfig::default(),
            false,
            dir.path().to_path_buf(),
            HashSet::new(),
            Endpoints {
                ena_fasta_base: format!("{base}/fasta"),
                ena_report_base: format!("{base}/assembly"),
                ..Default::default()
            },
        );
        assert!(asm_job.run(None).await.unwrap().is_complete());
        let asm_done = asm_job.commit().unwrap();
        let text = fs::read_to_string(&asm_done.files[0]).unwrap();
        assert!(text.contains("CM000001.1"));
        assert!(!text.contains("JACO02"));
    }

    struct CapturingReporter {
        updates: Mutex<Vec<u32>>,
    }

    impl DownloadProgressReporter for CapturingReporter {
        fn on_progress(&self, progress: DownloadProgress) {
            self.updates.lock().unwrap().push(progress.done);
        }
    }

    #[tokio::test]
    async fn test_mixed_genbank_and_ena_sequences_stream_progress_concurrently() {
        use std::sync::atomic::AtomicUsize;

        let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
        let addr = listener.local_addr().unwrap();
        let in_flight = Arc::new(AtomicUsize::new(0));
        let peak_in_flight = Arc::new(AtomicUsize::new(0));
        let in_flight_srv = Arc::clone(&in_flight);
        let peak_srv = Arc::clone(&peak_in_flight);

        tokio::spawn(async move {
            loop {
                let Ok((mut sock, _)) = listener.accept().await else {
                    break;
                };
                let active = Arc::clone(&in_flight_srv);
                let peak = Arc::clone(&peak_srv);
                tokio::spawn(async move {
                    let cur = active.fetch_add(1, Ordering::SeqCst) + 1;
                    peak.fetch_max(cur, Ordering::SeqCst);
                    let mut buf = [0u8; 4096];
                    let n = sock.read(&mut buf).await.unwrap_or(0);
                    let req = String::from_utf8_lossy(&buf[..n]);
                    let path = req.lines().next().unwrap_or("");
                    let body = if path.contains("efetch.fcgi") {
                        tokio::time::sleep(Duration::from_millis(80)).await;
                        ">U49845.1 Mus musculus\nACGT\n"
                    } else {
                        tokio::time::sleep(Duration::from_millis(10)).await;
                        ">ENA|AY123456|AY123456.1 Mus musculus\nTGCA\n"
                    };
                    let hdr = format!(
                        "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{}",
                        body.len(),
                        body
                    );
                    let _ = sock.write_all(hdr.as_bytes()).await;
                    active.fetch_sub(1, Ordering::SeqCst);
                });
            }
        });

        let dir = TempDir::new().unwrap();
        let base = format!("http://127.0.0.1:{}", addr.port());
        let endpoints = Endpoints {
            genbank_efetch: format!("{base}/efetch.fcgi"),
            ena_fasta_base: format!("{base}/fasta"),
            ..Default::default()
        };
        let records = vec![
            DownloadRecord {
                source: DatabaseSource::GenBank,
                content: DownloadContent::Nucleotide,
                accession: "U49845.1".to_string(),
                organism: "Mus musculus".to_string(),
                title: String::new(),
                length: Some(4),
                codon_start: 1,
                gene: String::new(),
                files: Vec::new(),
            },
            DownloadRecord {
                source: DatabaseSource::Ena,
                content: DownloadContent::Nucleotide,
                accession: "AY123456.1".to_string(),
                organism: "Mus musculus".to_string(),
                title: String::new(),
                length: Some(4),
                codon_start: 1,
                gene: String::new(),
                files: Vec::new(),
            },
        ];
        let reporter = Arc::new(CapturingReporter {
            updates: Mutex::new(Vec::new()),
        });
        let mut job = DownloadJob::with_client(
            loopback_client(true),
            records,
            DownloadOptions::default(),
            DownloadConfig::default(),
            false,
            dir.path().to_path_buf(),
            HashSet::new(),
            endpoints,
        );
        let report = job.run(Some(reporter.clone())).await.unwrap();
        assert!(report.is_complete());
        let completed = job.commit().unwrap();
        let fasta = fs::read_to_string(&completed.files[0]).unwrap();
        assert!(fasta.contains("U49845.1"));
        assert!(fasta.contains("AY123456.1"));
        assert!(peak_in_flight.load(Ordering::SeqCst) >= 2);
        let seen = reporter.updates.lock().unwrap().clone();
        assert!(seen.contains(&1));
        assert!(seen.contains(&2));
    }

    #[tokio::test]
    async fn test_genbank_assembly_parallel_work_stealing() {
        use std::sync::atomic::AtomicUsize;

        let large_fasta: Vec<u8> = b">contig1\nACGTACGT\n".repeat(500_000);
        let mut gz_buf = Vec::new();
        {
            let mut enc = GzEncoder::new(&mut gz_buf, Compression::new(0));
            enc.write_all(&large_fasta).unwrap();
            enc.finish().unwrap();
        }
        assert!(gz_buf.len() >= 8 * 1024 * 1024);
        let expected_md5 = format!("{:x}", Md5::digest(&gz_buf));
        let md5_txt = format!("{expected_md5}  ./GCA_001.fna.gz\n");
        let gz_len = gz_buf.len();

        let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
        let addr = listener.local_addr().unwrap();
        let range_gets = Arc::new(AtomicUsize::new(0));
        let range_srv = Arc::clone(&range_gets);
        let gz_srv = Arc::new(gz_buf);

        tokio::spawn(async move {
            loop {
                let Ok((mut sock, _)) = listener.accept().await else {
                    break;
                };
                let gets = Arc::clone(&range_srv);
                let data = Arc::clone(&gz_srv);
                let md5_body = md5_txt.clone();
                tokio::spawn(async move {
                    let mut buf = [0u8; 2048];
                    let n = sock.read(&mut buf).await.unwrap_or(0);
                    let req = String::from_utf8_lossy(&buf[..n]);
                    let first = req.lines().next().unwrap_or("");
                    if first.contains("md5checksums.txt") {
                        let hdr = format!(
                            "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{}",
                            md5_body.len(),
                            md5_body
                        );
                        let _ = sock.write_all(hdr.as_bytes()).await;
                    } else if first.starts_with("HEAD ") {
                        let hdr = format!(
                            "HTTP/1.1 200 OK\r\nAccept-Ranges: bytes\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                            data.len()
                        );
                        let _ = sock.write_all(hdr.as_bytes()).await;
                    } else if let Some(line) = req
                        .lines()
                        .find(|l| l.to_ascii_lowercase().starts_with("range:"))
                    {
                        gets.fetch_add(1, Ordering::SeqCst);
                        let spec = line.split('=').nth(1).unwrap().trim();
                        let mut parts = spec.split('-');
                        let start: usize = parts.next().unwrap().parse().unwrap();
                        let end: usize = parts.next().unwrap().parse().unwrap();
                        let slice = &data[start..=end];
                        let hdr = format!(
                            "HTTP/1.1 206 Partial Content\r\nContent-Length: {}\r\nContent-Range: bytes {}-{}/{}\r\nConnection: close\r\n\r\n",
                            slice.len(),
                            start,
                            end,
                            data.len()
                        );
                        let _ = sock.write_all(hdr.as_bytes()).await;
                        let _ = sock.write_all(slice).await;
                    }
                });
            }
        });

        let dir = TempDir::new().unwrap();
        let url = format!("http://127.0.0.1:{}/GCA_001.fna.gz", addr.port());
        let rec = DownloadRecord {
            source: DatabaseSource::GenBank,
            content: DownloadContent::Assembly,
            accession: "GCA_001.1".to_string(),
            organism: "Mus".to_string(),
            title: String::new(),
            length: None,
            codon_start: 1,
            gene: String::new(),
            files: vec![RemoteFile {
                url,
                bytes: Some(gz_len as u64),
                md5: None,
            }],
        };
        let mut job = DownloadJob::with_client(
            loopback_client(false),
            vec![rec],
            DownloadOptions::default(),
            DownloadConfig::default(),
            false,
            dir.path().to_path_buf(),
            HashSet::new(),
            Endpoints::default(),
        );
        let report = job.run(None).await.unwrap();
        assert!(report.is_complete());
        let completed = job.commit().unwrap();
        assert!(range_gets.load(Ordering::SeqCst) >= 3);
        assert_eq!(fs::read(&completed.files[0]).unwrap(), large_fasta);
    }
}
