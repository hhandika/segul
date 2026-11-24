//! Pooled, resumable, multi-range HTTPS client for NCBI and EMBL-EBI.

use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::time::{Duration, Instant};

use anyhow::Result;
use md5::{Digest, Md5};
use reqwest::header::{ACCEPT_RANGES, CONTENT_LENGTH, CONTENT_RANGE, LOCATION, RANGE, RETRY_AFTER};
use reqwest::{Client, StatusCode, Url};
use tokio::sync::{Mutex, Notify};
use tokio::time::{sleep, timeout};

use super::types::{DownloadCancelled, DownloadError};

const CONNECT_TIMEOUT: Duration = Duration::from_secs(20);
const RESPONSE_TIMEOUT: Duration = Duration::from_secs(45);
const READ_TIMEOUT: Duration = Duration::from_secs(30);
const MAX_ATTEMPTS: u32 = 3;
const MAX_REDIRECTS: usize = 4;
const RETRY_BASE_DELAY: Duration = Duration::from_secs(2);
const PARALLEL_THRESHOLD_BYTES: u64 = 8 * 1024 * 1024;
const PARALLEL_PARTS: u64 = 4;
const WRITE_BUF_SIZE: usize = 256 * 1024;

/// Cooperative cancellation signal shared across concurrent tasks.
#[derive(Debug, Clone, Default)]
pub struct CancelSignal {
    cancelled: Arc<AtomicBool>,
    notify: Arc<Notify>,
}

impl CancelSignal {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn cancel(&self) {
        if !self.cancelled.swap(true, Ordering::SeqCst) {
            self.notify.notify_waiters();
        }
    }

    pub fn is_cancelled(&self) -> bool {
        self.cancelled.load(Ordering::Relaxed)
    }

    pub fn throw_if_cancelled(&self) -> Result<()> {
        if self.is_cancelled() {
            return Err(DownloadCancelled.into());
        }
        Ok(())
    }

    pub async fn wait(&self, duration: Duration) -> Result<()> {
        if self.is_cancelled() {
            return Err(DownloadCancelled.into());
        }
        tokio::select! {
            _ = sleep(duration) => self.throw_if_cancelled(),
            _ = self.notify.notified() => Err(DownloadCancelled.into()),
        }
    }

    async fn run_cancellable<T, F>(&self, fut: F) -> Result<T>
    where
        F: std::future::Future<Output = Result<T>>,
    {
        self.throw_if_cancelled()?;
        tokio::select! {
            res = fut => {
                self.throw_if_cancelled()?;
                res
            }
            _ = self.notify.notified() => Err(DownloadCancelled.into()),
        }
    }
}

/// Options controlling how [`ArchiveClient::save_file`] fetches and hashes a file.
#[derive(Debug, Clone, Default)]
pub struct SaveOptions {
    pub form: Option<Vec<(String, String)>>,
    pub allow_resume: bool,
    pub allow_parallel_ranges: bool,
    pub known_size: Option<u64>,
    pub compute_md5: bool,
}

/// Outcome of [`ArchiveClient::save_file`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SaveResult {
    pub bytes_written: u64,
    pub md5: Option<String>,
}

/// Whether `url` is an HTTPS address hosted by NCBI or EMBL-EBI.
pub fn is_archive_url(url: &Url) -> bool {
    if url.scheme() != "https" {
        return false;
    }
    let Some(host) = url.host_str() else {
        return false;
    };
    host.ends_with(".ncbi.nlm.nih.gov") || host.ends_with(".ebi.ac.uk")
}

struct HostThrottle {
    ncbi_interval: Duration,
    default_interval: Duration,
    next_slot: Mutex<HashMap<String, Instant>>,
}

impl HostThrottle {
    fn new(has_ncbi_key: bool) -> Self {
        Self {
            ncbi_interval: if has_ncbi_key {
                Duration::from_millis(110)
            } else {
                Duration::from_millis(350)
            },
            default_interval: Duration::from_millis(50),
            next_slot: Mutex::new(HashMap::new()),
        }
    }

    async fn wait_slot(&self, host: &str, cancel: &CancelSignal) -> Result<()> {
        let interval = if host.starts_with("eutils.ncbi") {
            self.ncbi_interval
        } else {
            self.default_interval
        };
        let sleep_for = {
            let mut map = self.next_slot.lock().await;
            let now = Instant::now();
            let slot = map.get(host).copied().filter(|t| *t > now).unwrap_or(now);
            map.insert(host.to_string(), slot + interval);
            slot.saturating_duration_since(now)
        };
        if !sleep_for.is_zero() {
            cancel.wait(sleep_for).await?;
        }
        Ok(())
    }

    async fn push_back(&self, host: &str, delay: Duration) {
        let mut map = self.next_slot.lock().await;
        let target = Instant::now() + delay;
        let entry = map.entry(host.to_string()).or_insert(target);
        if *entry < target {
            *entry = target;
        }
    }
}

/// Pooled HTTPS client with retry, `Range` resume, and parallel byte-range support.
pub struct ArchiveClient {
    http: Client,
    throttle: Arc<HostThrottle>,
    is_allowed: Arc<dyn Fn(&Url) -> bool + Send + Sync>,
    retry_base_delay: Duration,
}

impl ArchiveClient {
    /// Creates a production client validated against NCBI and EMBL-EBI domains.
    pub fn new(has_ncbi_key: bool) -> Result<Self> {
        Self::with_validator(has_ncbi_key, Arc::new(is_archive_url), RETRY_BASE_DELAY)
    }

    /// Creates a client with a custom URL validator and retry base delay (for tests).
    pub fn with_validator(
        has_ncbi_key: bool,
        is_allowed: Arc<dyn Fn(&Url) -> bool + Send + Sync>,
        retry_base_delay: Duration,
    ) -> Result<Self> {
        let http = Client::builder()
            .redirect(reqwest::redirect::Policy::none())
            .connect_timeout(CONNECT_TIMEOUT)
            .pool_max_idle_per_host(6)
            .pool_idle_timeout(Duration::from_secs(15))
            .user_agent("segul-core/0.25")
            .build()?;
        Ok(Self {
            http,
            throttle: Arc::new(HostThrottle::new(has_ncbi_key)),
            is_allowed,
            retry_base_delay,
        })
    }

    /// Reads an HTTP response body as text, returning `(status, body)`.
    pub async fn read_text(
        &self,
        url: &Url,
        form: Option<&[(String, String)]>,
        cancel: &CancelSignal,
        on_retry: &(dyn Fn(u32, u32) + Send + Sync),
    ) -> Result<(u16, String)> {
        self.with_retry(url, cancel, on_retry, |is_last| async move {
            let response = self.send_request(url, form, None, true, cancel).await?;
            let status = response.status();
            if !is_last && is_busy(status) {
                let retry_after = parse_retry_after(&response);
                let _ = response.bytes().await;
                return Err(AttemptErr::Busy {
                    status,
                    retry_after,
                });
            }
            let bytes = cancel
                .run_cancellable(async {
                    timeout(READ_TIMEOUT, response.bytes())
                        .await
                        .map_err(|_| anyhow::anyhow!("Timed out reading response"))?
                        .map_err(Into::into)
                })
                .await
                .map_err(AttemptErr::classify)?;
            Ok((
                status.as_u16(),
                String::from_utf8_lossy(&bytes).into_owned(),
            ))
        })
        .await
    }

    /// Downloads `url` to `dest`, optionally resuming, splitting ranges, or computing MD5.
    pub async fn save_file(
        &self,
        url: &Url,
        dest: &Path,
        options: &SaveOptions,
        cancel: &CancelSignal,
        on_bytes: &(dyn Fn(u64, Option<u64>) + Send + Sync),
        on_retry: &(dyn Fn(u32, u32) + Send + Sync),
    ) -> Result<SaveResult> {
        if options.allow_parallel_ranges
            && options.form.is_none()
            && let Some(total) = self
                .probe_parallel_size(url, options.known_size, cancel)
                .await?
            && total >= PARALLEL_THRESHOLD_BYTES
        {
            match self
                .save_parallel(url, dest, total, options.compute_md5, cancel, on_bytes)
                .await
            {
                Ok(result) => return Ok(result),
                Err(err) if err.downcast_ref::<DownloadCancelled>().is_some() => {
                    return Err(err);
                }
                Err(_) => {
                    // Fall back cleanly to single-stream download.
                    let _ = std::fs::remove_file(dest);
                }
            }
        }
        self.save_single(url, dest, options, cancel, on_bytes, on_retry)
            .await
    }

    async fn probe_parallel_size(
        &self,
        start: &Url,
        known_size: Option<u64>,
        cancel: &CancelSignal,
    ) -> Result<Option<u64>> {
        if known_size.is_some_and(|s| s < PARALLEL_THRESHOLD_BYTES) {
            return Ok(None);
        }
        let Ok(resp) = self.send_head(start, cancel).await else {
            return Ok(None);
        };
        if !resp.status().is_success() {
            return Ok(None);
        }
        let accepts = resp
            .headers()
            .get(ACCEPT_RANGES)
            .and_then(|v| v.to_str().ok())
            .is_some_and(|v| v.eq_ignore_ascii_case("bytes"));
        if !accepts {
            return Ok(None);
        }
        let len = resp
            .headers()
            .get(CONTENT_LENGTH)
            .and_then(|v| v.to_str().ok())
            .and_then(|v| v.parse::<u64>().ok())
            .or(known_size);
        Ok(len.filter(|&s| s >= PARALLEL_THRESHOLD_BYTES))
    }

    async fn save_parallel(
        &self,
        url: &Url,
        dest: &Path,
        total: u64,
        compute_md5: bool,
        cancel: &CancelSignal,
        on_bytes: &(dyn Fn(u64, Option<u64>) + Send + Sync),
    ) -> Result<SaveResult> {
        {
            let file = File::create(dest)?;
            file.set_len(total)?;
        }
        let received_total = Arc::new(AtomicU64::new(0));
        on_bytes(0, Some(total));

        let chunk = total / PARALLEL_PARTS;
        let mut handles = Vec::with_capacity(PARALLEL_PARTS as usize);
        for part in 0..PARALLEL_PARTS {
            let start = part * chunk;
            let end = if part + 1 == PARALLEL_PARTS {
                total - 1
            } else {
                (part + 1) * chunk - 1
            };
            let url = url.clone();
            let dest = dest.to_path_buf();
            let cancel = cancel.clone();
            let received = Arc::clone(&received_total);
            let http = self.http.clone();
            let is_allowed = Arc::clone(&self.is_allowed);
            let throttle = Arc::clone(&self.throttle);
            let retry_base = self.retry_base_delay;
            handles.push(tokio::spawn(async move {
                download_range_part(
                    &http,
                    &*is_allowed,
                    &throttle,
                    retry_base,
                    &url,
                    &dest,
                    start,
                    end,
                    total,
                    &cancel,
                    &received,
                )
                .await
            }));
        }

        // Poll progress ticks while parallel range tasks write directly to offsets.
        let poll_cancel = cancel.clone();
        let poll_bytes = Arc::clone(&received_total);
        let mut join_set = Box::pin(futures_all(handles));
        loop {
            tokio::select! {
                res = &mut join_set => {
                    res?;
                    break;
                }
                _ = sleep(Duration::from_millis(100)) => {
                    poll_cancel.throw_if_cancelled()?;
                    on_bytes(poll_bytes.load(Ordering::Relaxed), Some(total));
                }
            }
        }
        on_bytes(total, Some(total));
        let md5 = if compute_md5 {
            Some(compute_file_md5(dest, cancel).await?)
        } else {
            None
        };
        Ok(SaveResult {
            bytes_written: total,
            md5,
        })
    }

    async fn save_single(
        &self,
        url: &Url,
        dest: &Path,
        options: &SaveOptions,
        cancel: &CancelSignal,
        on_bytes: &(dyn Fn(u64, Option<u64>) + Send + Sync),
        on_retry: &(dyn Fn(u32, u32) + Send + Sync),
    ) -> Result<SaveResult> {
        let form = options.form.as_deref();
        self.with_retry(url, cancel, on_retry, |is_last| async move {
            let existing = if options.allow_resume && form.is_none() {
                std::fs::metadata(dest).map(|m| m.len()).unwrap_or(0)
            } else {
                0
            };
            let range_header = (existing > 0).then(|| format!("bytes={existing}-"));
            let mut response = self
                .send_request(url, form, range_header.as_deref(), true, cancel)
                .await?;
            let status = response.status();

            if status != StatusCode::OK && status != StatusCode::PARTIAL_CONTENT {
                let retry_after = parse_retry_after(&response);
                let _ = response.bytes().await;
                if !is_last && is_busy(status) {
                    return Err(AttemptErr::Busy {
                        status,
                        retry_after,
                    });
                }
                return Err(AttemptErr::Fatal(
                    DownloadError(status_message(status)).into(),
                ));
            }

            let is_resumed = status == StatusCode::PARTIAL_CONTENT && existing > 0;
            let total_size = if is_resumed {
                parse_content_range_total(&response)
                    .or_else(|| response.content_length().map(|len| len + existing))
                    .or(options.known_size)
            } else {
                response.content_length().or(options.known_size)
            };

            let file = if is_resumed {
                OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(dest)
                    .map_err(AttemptErr::Fs)?
            } else {
                File::create(dest).map_err(AttemptErr::Fs)?
            };
            let mut writer = BufWriter::with_capacity(WRITE_BUF_SIZE, file);
            let mut hasher = (options.compute_md5 && !is_resumed).then(Md5::new);
            let mut received = if is_resumed { existing } else { 0 };
            on_bytes(received, total_size);

            loop {
                let chunk = cancel
                    .run_cancellable(async {
                        timeout(READ_TIMEOUT, response.chunk())
                            .await
                            .map_err(|_| anyhow::anyhow!("Timed out reading stream"))?
                            .map_err(Into::into)
                    })
                    .await
                    .map_err(AttemptErr::classify)?;
                let Some(bytes) = chunk else {
                    break;
                };
                writer.write_all(&bytes).map_err(AttemptErr::Fs)?;
                if let Some(h) = hasher.as_mut() {
                    h.update(&bytes);
                }
                received = received.saturating_add(bytes.len() as u64);
                on_bytes(received, total_size);
            }
            writer.flush().map_err(AttemptErr::Fs)?;
            drop(writer);

            let md5 = match (options.compute_md5, hasher) {
                (false, _) => None,
                (true, Some(h)) => Some(format!("{:x}", h.finalize())),
                (true, None) => Some(
                    compute_file_md5(dest, cancel)
                        .await
                        .map_err(AttemptErr::classify)?,
                ),
            };
            Ok(SaveResult {
                bytes_written: received,
                md5,
            })
        })
        .await
    }

    async fn with_retry<T, F, Fut>(
        &self,
        url: &Url,
        cancel: &CancelSignal,
        on_retry: &(dyn Fn(u32, u32) + Send + Sync),
        mut attempt: F,
    ) -> Result<T>
    where
        F: FnMut(bool) -> Fut,
        Fut: std::future::Future<Output = std::result::Result<T, AttemptErr>>,
    {
        for count in 1..=MAX_ATTEMPTS {
            cancel.throw_if_cancelled()?;
            let is_last = count >= MAX_ATTEMPTS;
            match attempt(is_last).await {
                Ok(val) => return Ok(val),
                Err(AttemptErr::Cancelled) => return Err(DownloadCancelled.into()),
                Err(AttemptErr::Fatal(err)) => return Err(err),
                Err(AttemptErr::Fs(err)) => {
                    cancel.throw_if_cancelled()?;
                    return Err(
                        DownloadError(format!("Could not save the download ({err}).")).into(),
                    );
                }
                Err(AttemptErr::Busy {
                    status,
                    retry_after,
                }) => {
                    if is_last {
                        return Err(DownloadError(status_message(status)).into());
                    }
                    let delay = retry_after.unwrap_or(self.retry_base_delay * count);
                    if let Some(host) = url.host_str() {
                        self.throttle.push_back(host, delay).await;
                    }
                    on_retry(count + 1, MAX_ATTEMPTS);
                    cancel.wait(delay).await?;
                }
                Err(AttemptErr::Transient(msg)) => {
                    if is_last {
                        return Err(DownloadError(msg).into());
                    }
                    on_retry(count + 1, MAX_ATTEMPTS);
                    cancel.wait(self.retry_base_delay * count).await?;
                }
            }
        }
        unreachable!()
    }

    async fn send_request(
        &self,
        start: &Url,
        form: Option<&[(String, String)]>,
        range: Option<&str>,
        apply_throttle: bool,
        cancel: &CancelSignal,
    ) -> std::result::Result<reqwest::Response, AttemptErr> {
        let mut url = start.clone();
        let mut send_form = form;
        for _ in 0..=MAX_REDIRECTS {
            if !(self.is_allowed)(&url) {
                return Err(AttemptErr::Fatal(
                    DownloadError("The database gave an unsafe address.".to_string()).into(),
                ));
            }
            if apply_throttle && let Some(host) = url.host_str() {
                self.throttle
                    .wait_slot(host, cancel)
                    .await
                    .map_err(AttemptErr::classify)?;
            }
            cancel.throw_if_cancelled().map_err(AttemptErr::classify)?;

            let mut req = match send_form {
                Some(pairs) => self.http.post(url.clone()).form(pairs),
                None => self.http.get(url.clone()),
            };
            if let Some(range_val) = range {
                req = req.header(RANGE, range_val);
            }
            let response = cancel
                .run_cancellable(async {
                    timeout(RESPONSE_TIMEOUT, req.send())
                        .await
                        .map_err(|_| anyhow::anyhow!("Timed out waiting for response"))?
                        .map_err(Into::into)
                })
                .await
                .map_err(AttemptErr::classify)?;

            if response.status().is_redirection()
                && let Some(loc) = response
                    .headers()
                    .get(LOCATION)
                    .and_then(|v| v.to_str().ok())
            {
                let next_url = url.join(loc).map_err(|_| {
                    AttemptErr::Fatal(
                        DownloadError("The database gave an unsafe address.".to_string()).into(),
                    )
                })?;
                if response.status() == StatusCode::SEE_OTHER {
                    send_form = None;
                }
                let _ = response.bytes().await;
                url = next_url;
                continue;
            }
            return Ok(response);
        }
        Err(AttemptErr::Fatal(
            DownloadError("The database redirected too many times.".to_string()).into(),
        ))
    }

    async fn send_head(
        &self,
        start: &Url,
        cancel: &CancelSignal,
    ) -> std::result::Result<reqwest::Response, AttemptErr> {
        let mut url = start.clone();
        for _ in 0..=MAX_REDIRECTS {
            if !(self.is_allowed)(&url) {
                return Err(AttemptErr::Fatal(
                    DownloadError("The database gave an unsafe address.".to_string()).into(),
                ));
            }
            let response = cancel
                .run_cancellable(async {
                    timeout(CONNECT_TIMEOUT, self.http.head(url.clone()).send())
                        .await
                        .map_err(|_| anyhow::anyhow!("Timed out probing file size"))?
                        .map_err(Into::into)
                })
                .await
                .map_err(AttemptErr::classify)?;
            if response.status().is_redirection()
                && let Some(loc) = response
                    .headers()
                    .get(LOCATION)
                    .and_then(|v| v.to_str().ok())
            {
                let next_url = url.join(loc).map_err(|_| {
                    AttemptErr::Fatal(
                        DownloadError("The database gave an unsafe address.".to_string()).into(),
                    )
                })?;
                url = next_url;
                continue;
            }
            return Ok(response);
        }
        Err(AttemptErr::Fatal(
            DownloadError("The database redirected too many times.".to_string()).into(),
        ))
    }
}

/// Computes lowercase hex MD5 of `path` on a blocking thread with cancellation checks.
pub async fn compute_file_md5(path: &Path, cancel: &CancelSignal) -> Result<String> {
    let path = path.to_path_buf();
    let cancel = cancel.clone();
    tokio::task::spawn_blocking(move || {
        let mut file = File::open(&path)?;
        let mut hasher = Md5::new();
        let mut buf = vec![0u8; WRITE_BUF_SIZE];
        loop {
            cancel.throw_if_cancelled()?;
            let n = file.read(&mut buf)?;
            if n == 0 {
                break;
            }
            hasher.update(&buf[..n]);
        }
        Ok(format!("{:x}", hasher.finalize()))
    })
    .await?
}

#[allow(clippy::too_many_arguments)]
async fn download_range_part(
    http: &Client,
    is_allowed: &(dyn Fn(&Url) -> bool + Send + Sync),
    throttle: &HostThrottle,
    retry_base: Duration,
    start_url: &Url,
    dest: &Path,
    range_start: u64,
    range_end: u64,
    _total_size: u64,
    cancel: &CancelSignal,
    received_total: &AtomicU64,
) -> Result<()> {
    let mut offset = range_start;
    for attempt in 1..=MAX_ATTEMPTS {
        cancel.throw_if_cancelled()?;
        if offset > range_end {
            return Ok(());
        }
        let mut url = start_url.clone();
        let mut response_opt = None;
        for _ in 0..=MAX_REDIRECTS {
            if !is_allowed(&url) {
                return Err(
                    DownloadError("The database gave an unsafe address.".to_string()).into(),
                );
            }
            let req = http
                .get(url.clone())
                .header(RANGE, format!("bytes={offset}-{range_end}"));
            let resp = cancel
                .run_cancellable(async {
                    timeout(RESPONSE_TIMEOUT, req.send())
                        .await
                        .map_err(|_| anyhow::anyhow!("Timed out waiting for range"))?
                        .map_err(Into::into)
                })
                .await?;
            if resp.status().is_redirection()
                && let Some(loc) = resp.headers().get(LOCATION).and_then(|v| v.to_str().ok())
            {
                url = url.join(loc)?;
                continue;
            }
            response_opt = Some(resp);
            break;
        }

        let Some(mut response) = response_opt else {
            return Err(
                DownloadError("The database redirected too many times.".to_string()).into(),
            );
        };
        let status = response.status();
        if status != StatusCode::PARTIAL_CONTENT {
            if attempt < MAX_ATTEMPTS && is_busy(status) {
                let delay = parse_retry_after(&response).unwrap_or(retry_base * attempt);
                if let Some(host) = url.host_str() {
                    throttle.push_back(host, delay).await;
                }
                cancel.wait(delay).await?;
                continue;
            }
            return Err(anyhow::anyhow!(
                "Server did not serve partial range ({status})"
            ));
        }

        let mut file = OpenOptions::new().write(true).open(dest)?;
        file.seek(SeekFrom::Start(offset))?;
        let mut writer = BufWriter::with_capacity(WRITE_BUF_SIZE, file);
        let mut stream_err: Option<anyhow::Error> = None;
        loop {
            let next = cancel
                .run_cancellable(async {
                    timeout(READ_TIMEOUT, response.chunk())
                        .await
                        .map_err(|_| anyhow::anyhow!("Timed out reading range"))?
                        .map_err(Into::into)
                })
                .await;
            match next {
                Ok(Some(chunk)) => {
                    let remaining = (range_end - offset + 1) as usize;
                    let take = chunk.len().min(remaining);
                    writer.write_all(&chunk[..take])?;
                    offset += take as u64;
                    received_total.fetch_add(take as u64, Ordering::Relaxed);
                    if offset > range_end {
                        break;
                    }
                }
                Ok(None) => break,
                Err(err) => {
                    stream_err = Some(err);
                    break;
                }
            }
        }
        writer.flush()?;
        if let Some(err) = stream_err {
            if err.downcast_ref::<DownloadCancelled>().is_some() || attempt >= MAX_ATTEMPTS {
                return Err(err);
            }
            cancel.wait(retry_base * attempt).await?;
            continue;
        }
        if offset > range_end {
            return Ok(());
        }
        if attempt >= MAX_ATTEMPTS {
            return Err(
                DownloadError("The connection was interrupted. Try again.".to_string()).into(),
            );
        }
        cancel.wait(retry_base * attempt).await?;
    }
    Ok(())
}

async fn futures_all(handles: Vec<tokio::task::JoinHandle<Result<()>>>) -> Result<()> {
    for handle in handles {
        handle.await??;
    }
    Ok(())
}

enum AttemptErr {
    Cancelled,
    Fatal(anyhow::Error),
    Fs(std::io::Error),
    Busy {
        status: StatusCode,
        retry_after: Option<Duration>,
    },
    Transient(String),
}

impl AttemptErr {
    fn classify(err: anyhow::Error) -> Self {
        if err.downcast_ref::<DownloadCancelled>().is_some() {
            return Self::Cancelled;
        }
        if err.downcast_ref::<DownloadError>().is_some() {
            return Self::Fatal(err);
        }
        let text = err.to_string();
        if text.contains("Timed out") {
            return Self::Transient("The database took too long to answer.".to_string());
        }
        if let Some(req_err) = err.downcast_ref::<reqwest::Error>() {
            if req_err.is_timeout() {
                return Self::Transient("The database took too long to answer.".to_string());
            }
            if req_err.is_connect() {
                return Self::Transient(
                    "Could not connect. Check your internet connection.".to_string(),
                );
            }
        }
        Self::Transient("The connection was interrupted. Try again.".to_string())
    }
}

fn is_busy(status: StatusCode) -> bool {
    status == StatusCode::TOO_MANY_REQUESTS || status.is_server_error()
}

fn status_message(status: StatusCode) -> String {
    if status == StatusCode::NOT_FOUND {
        "Not found in the database.".to_string()
    } else if is_busy(status) {
        "The database is busy. Try again later.".to_string()
    } else {
        format!("The database refused the request ({}).", status.as_u16())
    }
}

fn parse_retry_after(response: &reqwest::Response) -> Option<Duration> {
    let val = response.headers().get(RETRY_AFTER)?.to_str().ok()?.trim();
    val.parse::<u64>()
        .ok()
        .map(|s| Duration::from_secs(s.min(60)))
}

fn parse_content_range_total(response: &reqwest::Response) -> Option<u64> {
    let val = response.headers().get(CONTENT_RANGE)?.to_str().ok()?;
    let slash = val.rfind('/')?;
    val[slash + 1..].trim().parse::<u64>().ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicUsize;
    use tempfile::TempDir;
    use tokio::io::{AsyncReadExt, AsyncWriteExt};
    use tokio::net::TcpListener;

    fn test_client() -> ArchiveClient {
        ArchiveClient::with_validator(
            false,
            Arc::new(|u| u.host_str() == Some("127.0.0.1")),
            Duration::from_millis(10),
        )
        .unwrap()
    }

    #[tokio::test]
    async fn test_range_resume_and_md5() {
        let payload: Vec<u8> = (0..120_000).map(|i| (i % 251) as u8).collect();
        let expected_md5 = format!("{:x}", Md5::digest(&payload));
        let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
        let addr = listener.local_addr().unwrap();
        let calls = Arc::new(AtomicUsize::new(0));
        let calls_srv = Arc::clone(&calls);
        let payload_srv = payload.clone();

        tokio::spawn(async move {
            loop {
                let Ok((mut sock, _)) = listener.accept().await else {
                    break;
                };
                let c = calls_srv.fetch_add(1, Ordering::SeqCst);
                let data = payload_srv.clone();
                tokio::spawn(async move {
                    let mut buf = [0u8; 2048];
                    let n = sock.read(&mut buf).await.unwrap_or(0);
                    let req = String::from_utf8_lossy(&buf[..n]);
                    if c == 0 {
                        let half = data.len() / 2;
                        let hdr = format!(
                            "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                            data.len()
                        );
                        let _ = sock.write_all(hdr.as_bytes()).await;
                        let _ = sock.write_all(&data[..half]).await;
                        let _ = sock.shutdown().await;
                    } else if let Some(line) = req
                        .lines()
                        .find(|l| l.to_ascii_lowercase().starts_with("range:"))
                    {
                        let start: usize = line
                            .split('=')
                            .nth(1)
                            .unwrap()
                            .trim_end_matches('-')
                            .trim()
                            .parse()
                            .unwrap();
                        let slice = &data[start..];
                        let hdr = format!(
                            "HTTP/1.1 206 Partial Content\r\nContent-Length: {}\r\nContent-Range: bytes {}-{}/{}\r\nConnection: close\r\n\r\n",
                            slice.len(),
                            start,
                            data.len() - 1,
                            data.len()
                        );
                        let _ = sock.write_all(hdr.as_bytes()).await;
                        let _ = sock.write_all(slice).await;
                    }
                });
            }
        });

        let dir = TempDir::new().unwrap();
        let dest = dir.path().join("resumed.bin");
        let url = Url::parse(&format!("http://127.0.0.1:{}/file", addr.port())).unwrap();
        let retries = Arc::new(AtomicUsize::new(0));
        let retries_cb = Arc::clone(&retries);
        let res = test_client()
            .save_file(
                &url,
                &dest,
                &SaveOptions {
                    allow_resume: true,
                    compute_md5: true,
                    ..Default::default()
                },
                &CancelSignal::new(),
                &|_, _| {},
                &move |_, _| {
                    retries_cb.fetch_add(1, Ordering::SeqCst);
                },
            )
            .await
            .unwrap();
        assert_eq!(calls.load(Ordering::SeqCst), 2);
        assert_eq!(retries.load(Ordering::SeqCst), 1);
        assert_eq!(res.bytes_written as usize, payload.len());
        assert_eq!(res.md5.as_deref(), Some(expected_md5.as_str()));
        assert_eq!(std::fs::read(&dest).unwrap(), payload);
    }

    #[tokio::test]
    async fn test_parallel_four_ranges_download() {
        let payload: Vec<u8> = (0..9_000_000).map(|i| (i % 241) as u8).collect();
        let expected_md5 = format!("{:x}", Md5::digest(&payload));
        let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
        let addr = listener.local_addr().unwrap();
        let range_gets = Arc::new(AtomicUsize::new(0));
        let range_gets_srv = Arc::clone(&range_gets);
        let payload_srv = Arc::new(payload.clone());

        tokio::spawn(async move {
            loop {
                let Ok((mut sock, _)) = listener.accept().await else {
                    break;
                };
                let gets = Arc::clone(&range_gets_srv);
                let data = Arc::clone(&payload_srv);
                tokio::spawn(async move {
                    let mut buf = [0u8; 2048];
                    let n = sock.read(&mut buf).await.unwrap_or(0);
                    let req = String::from_utf8_lossy(&buf[..n]);
                    if req.starts_with("HEAD ") {
                        let hdr = format!(
                            "HTTP/1.1 200 OK\r\nAccept-Ranges: bytes\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                            data.len()
                        );
                        let _ = sock.write_all(hdr.as_bytes()).await;
                        return;
                    }
                    if let Some(line) = req
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
        let dest = dir.path().join("parallel.bin");
        let url = Url::parse(&format!("http://127.0.0.1:{}/big.fastq.gz", addr.port())).unwrap();
        let res = test_client()
            .save_file(
                &url,
                &dest,
                &SaveOptions {
                    allow_resume: true,
                    allow_parallel_ranges: true,
                    compute_md5: true,
                    ..Default::default()
                },
                &CancelSignal::new(),
                &|_, _| {},
                &|_, _| {},
            )
            .await
            .unwrap();
        assert_eq!(range_gets.load(Ordering::SeqCst), 4);
        assert_eq!(res.md5.as_deref(), Some(expected_md5.as_str()));
        assert_eq!(std::fs::read(&dest).unwrap(), payload);
    }

    #[tokio::test]
    async fn test_retry_after_and_keep_alive_reuse() {
        let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
        let addr = listener.local_addr().unwrap();
        let connections = Arc::new(AtomicUsize::new(0));
        let connections_srv = Arc::clone(&connections);
        let requests = Arc::new(AtomicUsize::new(0));
        let requests_srv = Arc::clone(&requests);

        tokio::spawn(async move {
            loop {
                let Ok((mut sock, _)) = listener.accept().await else {
                    break;
                };
                connections_srv.fetch_add(1, Ordering::SeqCst);
                let reqs = Arc::clone(&requests_srv);
                tokio::spawn(async move {
                    let mut buf = [0u8; 2048];
                    while let Ok(n) = sock.read(&mut buf).await {
                        if n == 0 {
                            break;
                        }
                        let num = reqs.fetch_add(1, Ordering::SeqCst);
                        if num == 0 {
                            let resp = "HTTP/1.1 429 Too Many Requests\r\nRetry-After: 0\r\nContent-Length: 4\r\n\r\nbusy";
                            let _ = sock.write_all(resp.as_bytes()).await;
                        } else {
                            let resp = "HTTP/1.1 200 OK\r\nContent-Length: 2\r\n\r\nok";
                            let _ = sock.write_all(resp.as_bytes()).await;
                        }
                    }
                });
            }
        });

        let client = test_client();
        let url = Url::parse(&format!("http://127.0.0.1:{}/ping", addr.port())).unwrap();
        let retries = Arc::new(AtomicUsize::new(0));
        let retries_cb = Arc::clone(&retries);
        let (status, body) = client
            .read_text(&url, None, &CancelSignal::new(), &move |_, _| {
                retries_cb.fetch_add(1, Ordering::SeqCst);
            })
            .await
            .unwrap();
        assert_eq!(status, 200);
        assert_eq!(body, "ok");
        assert_eq!(retries.load(Ordering::SeqCst), 1);

        let (status2, body2) = client
            .read_text(&url, None, &CancelSignal::new(), &|_, _| {})
            .await
            .unwrap();
        assert_eq!((status2, body2.as_str()), (200, "ok"));
        assert_eq!(connections.load(Ordering::SeqCst), 1);
        assert_eq!(requests.load(Ordering::SeqCst), 3);
    }
}
