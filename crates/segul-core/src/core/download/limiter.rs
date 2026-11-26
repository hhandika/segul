//! Per-host concurrency budget, request pacing, and adaptive congestion control.

use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, MutexGuard};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use anyhow::Result;
use chrono::{DateTime, Utc};
use tokio::sync::Notify;

use super::client::CancelSignal;

const MAX_RETRY_AFTER_SECS: i64 = 60;
const MAX_BACKOFF: Duration = Duration::from_secs(30);
const MAX_ADAPTIVE_INTERVAL: Duration = Duration::from_secs(2);
const INTERVAL_SNAP_THRESHOLD: Duration = Duration::from_millis(5);

static JITTER_SEQ: AtomicU64 = AtomicU64::new(0);

/// Controls how aggressively the downloader opens parallel connections per host.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum ConcurrencyProfile {
    /// Uses up to the safe maximum allowed by each archive host.
    #[default]
    Auto,
    /// Uses half the safe host limits for shared IPs or slower connections.
    Conservative,
    /// Caps per-host concurrency at `n`, clamped to each host's safe maximum.
    Max(usize),
}

/// Archive host category determining safe ceiling and baseline spacing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HostClass {
    NcbiEutils,
    NcbiBulk,
    EnaApi,
    EnaBulk,
    Other,
}

impl HostClass {
    /// Classifies a hostname for NCBI and EMBL-EBI rate and connection limits.
    pub fn classify(host: &str) -> Self {
        let lower = host.trim().to_ascii_lowercase();
        if lower == "eutils.ncbi.nlm.nih.gov" || lower.starts_with("eutils.ncbi.") {
            Self::NcbiEutils
        } else if lower.ends_with(".ncbi.nlm.nih.gov") || lower == "ncbi.nlm.nih.gov" {
            Self::NcbiBulk
        } else if lower == "www.ebi.ac.uk" {
            Self::EnaApi
        } else if lower.ends_with(".ebi.ac.uk") || lower == "ebi.ac.uk" {
            Self::EnaBulk
        } else {
            Self::Other
        }
    }

    /// Returns `(safe_max_connections, base_request_interval)` for this host class.
    pub fn policy(self, has_ncbi_key: bool) -> (usize, Duration) {
        match self {
            Self::NcbiEutils if has_ncbi_key => (8, Duration::from_millis(115)),
            Self::NcbiEutils => (3, Duration::from_millis(350)),
            Self::NcbiBulk => (4, Duration::from_millis(250)),
            Self::EnaApi | Self::EnaBulk => (8, Duration::from_millis(100)),
            Self::Other => (8, Duration::from_millis(10)),
        }
    }
}

struct HostState {
    active: usize,
    limit: usize,
    max_limit: usize,
    base_interval: Duration,
    current_interval: Duration,
    next_slot: Instant,
    backoff_until: Option<Instant>,
    clean_streak: usize,
    transient_streak: usize,
}

struct HostGate {
    state: Mutex<HostState>,
    notify: Notify,
}

/// RAII connection permit returned by [`HostLimiter::acquire`].
pub struct HostPermit {
    gate: Arc<HostGate>,
}

impl Drop for HostPermit {
    fn drop(&mut self) {
        let mut state = lock_gate(&self.gate);
        state.active = state.active.saturating_sub(1);
        drop(state);
        self.gate.notify.notify_waiters();
    }
}

/// Adaptive per-host concurrency gate and request pacer shared across a job.
pub struct HostLimiter {
    has_ncbi_key: bool,
    profile: ConcurrencyProfile,
    gates: Mutex<HashMap<String, Arc<HostGate>>>,
}

impl HostLimiter {
    /// Creates a per-host limiter configured for `has_ncbi_key` and `profile`.
    pub fn new(has_ncbi_key: bool, profile: ConcurrencyProfile) -> Self {
        Self {
            has_ncbi_key,
            profile,
            gates: Mutex::new(HashMap::new()),
        }
    }

    /// Returns the hard-clamped safe concurrency ceiling for `host`.
    pub fn max_limit(&self, host: &str) -> usize {
        let gate = self.gate_for(host);
        let state = lock_gate(&gate);
        state.max_limit
    }

    /// Returns the current adaptive concurrency limit for `host`.
    pub fn current_limit(&self, host: &str) -> usize {
        let gate = self.gate_for(host);
        let state = lock_gate(&gate);
        state.limit
    }

    /// Allocates range workers for a file given `concurrent_files` sharing `host`.
    pub fn range_workers_for(&self, host: &str, concurrent_files: usize) -> usize {
        (self.current_limit(host) / concurrent_files.max(1)).max(1)
    }

    /// Waits cancellably for an active connection slot on `host`.
    pub async fn acquire(&self, host: &str, cancel: &CancelSignal) -> Result<HostPermit> {
        let gate = self.gate_for(host);
        loop {
            cancel.throw_if_cancelled()?;
            let notified = gate.notify.notified();
            tokio::pin!(notified);
            notified.as_mut().enable();
            {
                let mut state = lock_gate(&gate);
                if state.active < state.limit {
                    state.active += 1;
                    return Ok(HostPermit {
                        gate: Arc::clone(&gate),
                    });
                }
            }
            cancel
                .run_cancellable(async {
                    notified.await;
                    Ok(())
                })
                .await?;
        }
    }

    /// Reserves the next paced request slot for `host` and waits cancellably if needed.
    pub async fn pace(&self, host: &str, cancel: &CancelSignal) -> Result<()> {
        let gate = self.gate_for(host);
        let sleep_for = {
            let mut state = lock_gate(&gate);
            let now = Instant::now();
            let floor = state.backoff_until.filter(|t| *t > now).unwrap_or(now);
            let slot = state.next_slot.max(floor).max(now);
            state.next_slot = slot + state.current_interval;
            slot.saturating_duration_since(now)
        };
        if !sleep_for.is_zero() {
            cancel.wait(sleep_for).await?;
        }
        Ok(())
    }

    /// Records a clean archive response and restores concurrency additively after a streak.
    pub fn on_success(&self, host: &str) {
        let gate = self.gate_for(host);
        let mut should_notify = false;
        {
            let mut state = lock_gate(&gate);
            state.transient_streak = 0;
            state.clean_streak = state.clean_streak.saturating_add(1);
            let threshold = (state.limit.saturating_mul(2)).max(4);
            if state.clean_streak >= threshold {
                state.clean_streak = 0;
                if state.limit < state.max_limit {
                    state.limit += 1;
                    should_notify = true;
                }
                if state.current_interval > state.base_interval {
                    let excess = state.current_interval.saturating_sub(state.base_interval);
                    if excess <= INTERVAL_SNAP_THRESHOLD {
                        state.current_interval = state.base_interval;
                    } else {
                        state.current_interval = state.base_interval + excess.mul_f32(0.8);
                    }
                }
            }
        }
        if should_notify {
            gate.notify.notify_waiters();
        }
    }

    /// Multiplicatively reduces concurrency and enforces backoff after HTTP 429 or 5xx.
    pub fn on_throttled(&self, host: &str, delay: Duration) {
        let gate = self.gate_for(host);
        let mut state = lock_gate(&gate);
        state.clean_streak = 0;
        state.transient_streak = 0;
        state.limit = (state.limit / 2).max(1);
        state.current_interval = state
            .current_interval
            .mul_f32(1.5)
            .min(MAX_ADAPTIVE_INTERVAL);
        let target = Instant::now() + delay;
        state.backoff_until = Some(state.backoff_until.map_or(target, |cur| cur.max(target)));
        if state.next_slot < target {
            state.next_slot = target;
        }
    }

    /// Steps concurrency down slightly after repeated transient network failures.
    pub fn on_transient_error(&self, host: &str) {
        let gate = self.gate_for(host);
        let mut state = lock_gate(&gate);
        state.clean_streak = 0;
        state.transient_streak = state.transient_streak.saturating_add(1);
        if state.transient_streak >= 2 {
            state.transient_streak = 0;
            state.limit = state.limit.saturating_sub(1).max(1);
        }
    }

    fn gate_for(&self, host: &str) -> Arc<HostGate> {
        let key = host.trim().to_ascii_lowercase();
        let mut map = self.gates.lock().unwrap_or_else(|e| e.into_inner());
        if let Some(existing) = map.get(&key) {
            return Arc::clone(existing);
        }
        let class = HostClass::classify(&key);
        let (safe_max, base_interval) = class.policy(self.has_ncbi_key);
        let max_limit = effective_max(safe_max, self.profile);
        let gate = Arc::new(HostGate {
            state: Mutex::new(HostState {
                active: 0,
                limit: max_limit,
                max_limit,
                base_interval,
                current_interval: base_interval,
                next_slot: Instant::now(),
                backoff_until: None,
                clean_streak: 0,
                transient_streak: 0,
            }),
            notify: Notify::new(),
        });
        map.insert(key, Arc::clone(&gate));
        gate
    }
}

/// Parses an HTTP `Retry-After` header value (integer seconds or RFC 2822 HTTP-date).
pub fn parse_retry_after_value(raw: &str) -> Option<Duration> {
    let trimmed = raw.trim();
    if let Ok(seconds) = trimmed.parse::<u64>() {
        return Some(Duration::from_secs(
            seconds.min(MAX_RETRY_AFTER_SECS as u64),
        ));
    }
    let dt = DateTime::parse_from_rfc2822(trimmed).ok()?;
    let secs = dt
        .with_timezone(&Utc)
        .signed_duration_since(Utc::now())
        .num_seconds()
        .clamp(0, MAX_RETRY_AFTER_SECS);
    Some(Duration::from_secs(secs as u64))
}

/// Computes exponential backoff with +/-20% jitter (or honors `retry_after` with light jitter).
pub fn jittered_backoff(base: Duration, attempt: u32, retry_after: Option<Duration>) -> Duration {
    let raw = match retry_after {
        Some(wait) => wait.min(Duration::from_secs(MAX_RETRY_AFTER_SECS as u64)),
        None => {
            let shift = attempt.saturating_sub(1).min(4);
            base.saturating_mul(1u32 << shift).min(MAX_BACKOFF)
        }
    };
    if raw < Duration::from_millis(50) {
        return raw;
    }
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.subsec_nanos() as u64);
    let seq = JITTER_SEQ.fetch_add(1, Ordering::Relaxed);
    let bucket = (nanos ^ seq.wrapping_mul(0x9E37_79B9_7F4A_7C15)) % 41;
    // Spread across [0.80, 1.20] so concurrent retries do not collide.
    let factor = 0.80 + (bucket as f32) * 0.01;
    raw.mul_f32(factor)
}

fn effective_max(safe_max: usize, profile: ConcurrencyProfile) -> usize {
    match profile {
        ConcurrencyProfile::Auto => safe_max,
        ConcurrencyProfile::Conservative => safe_max.div_ceil(2).max(1),
        ConcurrencyProfile::Max(n) => n.clamp(1, safe_max),
    }
}

fn lock_gate(gate: &HostGate) -> MutexGuard<'_, HostState> {
    gate.state.lock().unwrap_or_else(|e| e.into_inner())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicUsize;

    #[test]
    fn test_hard_clamping_respects_archive_policy() {
        let unkeyed = HostLimiter::new(false, ConcurrencyProfile::Max(100));
        assert_eq!(unkeyed.max_limit("eutils.ncbi.nlm.nih.gov"), 3);
        assert_eq!(unkeyed.max_limit("ftp.ncbi.nlm.nih.gov"), 4);
        assert_eq!(unkeyed.max_limit("www.ebi.ac.uk"), 8);
        assert_eq!(unkeyed.max_limit("ftp.sra.ebi.ac.uk"), 8);

        let keyed = HostLimiter::new(true, ConcurrencyProfile::Max(100));
        assert_eq!(keyed.max_limit("eutils.ncbi.nlm.nih.gov"), 8);
        assert_eq!(keyed.max_limit("ftp.ncbi.nlm.nih.gov"), 4);

        let mobile_auto = HostLimiter::new(true, ConcurrencyProfile::Max(4));
        assert_eq!(mobile_auto.max_limit("eutils.ncbi.nlm.nih.gov"), 4);
        assert_eq!(mobile_auto.max_limit("ftp.sra.ebi.ac.uk"), 4);
        assert_eq!(
            HostLimiter::new(false, ConcurrencyProfile::Max(4))
                .max_limit("eutils.ncbi.nlm.nih.gov"),
            3
        );

        let conservative = HostLimiter::new(true, ConcurrencyProfile::Conservative);
        assert_eq!(conservative.max_limit("eutils.ncbi.nlm.nih.gov"), 4);
        assert_eq!(conservative.max_limit("ftp.ncbi.nlm.nih.gov"), 2);
        assert_eq!(conservative.max_limit("ftp.sra.ebi.ac.uk"), 4);
    }

    #[test]
    fn test_aimd_halves_on_throttle_and_recovers_on_clean_streak() {
        let limiter = HostLimiter::new(true, ConcurrencyProfile::Auto);
        let host = "ftp.sra.ebi.ac.uk";
        assert_eq!(limiter.current_limit(host), 8);

        limiter.on_throttled(host, Duration::from_millis(100));
        assert_eq!(limiter.current_limit(host), 4);
        limiter.on_throttled(host, Duration::from_millis(100));
        assert_eq!(limiter.current_limit(host), 2);
        limiter.on_throttled(host, Duration::from_millis(100));
        assert_eq!(limiter.current_limit(host), 1);

        for _ in 0..4 {
            limiter.on_success(host);
        }
        assert_eq!(limiter.current_limit(host), 2);

        for _ in 0..4 {
            limiter.on_success(host);
        }
        assert_eq!(limiter.current_limit(host), 3);

        limiter.on_transient_error(host);
        assert_eq!(limiter.current_limit(host), 3);
        limiter.on_transient_error(host);
        assert_eq!(limiter.current_limit(host), 2);
    }

    #[test]
    fn test_range_worker_allocation_shares_host_budget() {
        let ena = HostLimiter::new(false, ConcurrencyProfile::Auto);
        assert_eq!(ena.range_workers_for("ftp.sra.ebi.ac.uk", 1), 8);
        assert_eq!(ena.range_workers_for("ftp.sra.ebi.ac.uk", 2), 4);
        assert_eq!(ena.range_workers_for("ftp.sra.ebi.ac.uk", 8), 1);
        assert_eq!(ena.range_workers_for("ftp.sra.ebi.ac.uk", 20), 1);

        let ncbi_ftp = HostLimiter::new(true, ConcurrencyProfile::Auto);
        assert_eq!(ncbi_ftp.range_workers_for("ftp.ncbi.nlm.nih.gov", 1), 4);
        assert_eq!(ncbi_ftp.range_workers_for("ftp.ncbi.nlm.nih.gov", 2), 2);
        assert_eq!(ncbi_ftp.range_workers_for("ftp.ncbi.nlm.nih.gov", 4), 1);
    }

    #[test]
    fn test_parse_retry_after_seconds_and_http_date() {
        assert_eq!(parse_retry_after_value("5"), Some(Duration::from_secs(5)));
        assert_eq!(
            parse_retry_after_value("999"),
            Some(Duration::from_secs(60))
        );
        let future = (Utc::now() + chrono::Duration::seconds(10)).to_rfc2822();
        let parsed = parse_retry_after_value(&future).unwrap();
        assert!(parsed.as_secs() >= 8 && parsed.as_secs() <= 11);
        assert_eq!(parse_retry_after_value("not-a-date"), None);
    }

    #[tokio::test]
    async fn test_acquire_enforces_peak_active_invariant() {
        let limiter = Arc::new(HostLimiter::new(false, ConcurrencyProfile::Max(3)));
        let cancel = CancelSignal::new();
        let active = Arc::new(AtomicUsize::new(0));
        let peak = Arc::new(AtomicUsize::new(0));

        let mut handles = Vec::new();
        for _ in 0..24 {
            let limiter = Arc::clone(&limiter);
            let cancel = cancel.clone();
            let active = Arc::clone(&active);
            let peak = Arc::clone(&peak);
            handles.push(tokio::spawn(async move {
                let _permit = limiter.acquire("ftp.sra.ebi.ac.uk", &cancel).await.unwrap();
                let now = active.fetch_add(1, Ordering::SeqCst) + 1;
                peak.fetch_max(now, Ordering::SeqCst);
                tokio::time::sleep(Duration::from_millis(5)).await;
                active.fetch_sub(1, Ordering::SeqCst);
            }));
        }
        for handle in handles {
            handle.await.unwrap();
        }
        assert!(peak.load(Ordering::SeqCst) <= 3);
        assert_eq!(active.load(Ordering::SeqCst), 0);
    }
}
