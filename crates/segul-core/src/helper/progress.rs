//! Per-file progress and cancellation for long-running tasks.
//!
//! Hosts such as a GUI implement [`ProgressReporter`] and pass it to a task's
//! `with_progress` builder to follow which files are done and to stop a run
//! between files.

use std::fmt;
use std::path::Path;

use anyhow::Result;

/// Receives per-file progress from long-running tasks.
///
/// Files may run in parallel, so methods are called from several threads.
/// All methods do nothing by default.
pub trait ProgressReporter: Send + Sync {
    /// Called once before any file starts.
    fn on_start(&self, _total_files: usize) {}

    fn on_file_started(&self, _path: &Path) {}

    fn on_file_finished(&self, _path: &Path) {}

    /// A step after the file loop, such as writing outputs.
    fn on_stage(&self, _stage: &str) {}

    /// Checked before each file starts. Returning true skips the remaining
    /// files and the task returns [`Cancelled`]. Files already running finish.
    fn is_cancelled(&self) -> bool {
        false
    }
}

/// Error returned when a [`ProgressReporter`] cancels a run.
///
/// Hosts can tell it apart from other errors with
/// `error.downcast_ref::<Cancelled>()`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Cancelled;

impl fmt::Display for Cancelled {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Run cancelled")
    }
}

impl std::error::Error for Cancelled {}

pub(crate) fn report_start(reporter: Option<&dyn ProgressReporter>, total_files: usize) {
    if let Some(reporter) = reporter {
        reporter.on_start(total_files);
    }
}

pub(crate) fn report_stage(reporter: Option<&dyn ProgressReporter>, stage: &str) {
    if let Some(reporter) = reporter {
        reporter.on_stage(stage);
    }
}

/// Returns [`Cancelled`] if the reporter asked to stop.
pub(crate) fn check_cancelled(reporter: Option<&dyn ProgressReporter>) -> Result<()> {
    match reporter {
        Some(reporter) if reporter.is_cancelled() => Err(Cancelled.into()),
        _ => Ok(()),
    }
}

/// Runs `task` for one file, reporting when it starts and when it finishes
/// successfully. Returns [`Cancelled`] without running it if the reporter
/// asked to stop.
pub(crate) fn track<T>(
    reporter: Option<&dyn ProgressReporter>,
    path: &Path,
    task: impl FnOnce() -> Result<T>,
) -> Result<T> {
    let Some(reporter) = reporter else {
        return task();
    };
    check_cancelled(Some(reporter))?;
    reporter.on_file_started(path);
    let result = task();
    if result.is_ok() {
        reporter.on_file_finished(path);
    }
    result
}

/// Reporter for tests that lets the first `allowed_checks` cancel checks
/// pass, then cancels, so a run stops partway through a file.
#[cfg(test)]
pub(crate) struct CancelAfter {
    allowed_checks: usize,
    checks: std::sync::atomic::AtomicUsize,
    pub(crate) finished: std::sync::atomic::AtomicUsize,
}

#[cfg(test)]
impl CancelAfter {
    pub(crate) fn new(allowed_checks: usize) -> Self {
        Self {
            allowed_checks,
            checks: Default::default(),
            finished: Default::default(),
        }
    }
}

#[cfg(test)]
impl ProgressReporter for CancelAfter {
    fn on_file_finished(&self, _path: &Path) {
        self.finished
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    }

    fn is_cancelled(&self) -> bool {
        let checks = self
            .checks
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        checks >= self.allowed_checks
    }
}

#[cfg(test)]
mod test {
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

    use super::*;

    #[derive(Default)]
    struct CountingReporter {
        started: AtomicUsize,
        finished: AtomicUsize,
        cancelled: AtomicBool,
    }

    impl ProgressReporter for CountingReporter {
        fn on_file_started(&self, _path: &Path) {
            self.started.fetch_add(1, Ordering::Relaxed);
        }

        fn on_file_finished(&self, _path: &Path) {
            self.finished.fetch_add(1, Ordering::Relaxed);
        }

        fn is_cancelled(&self) -> bool {
            self.cancelled.load(Ordering::Relaxed)
        }
    }

    #[test]
    fn test_track_reports_each_file() {
        let reporter = CountingReporter::default();
        let paths = [PathBuf::from("a.fastq"), PathBuf::from("b.fastq")];
        for path in &paths {
            track(Some(&reporter), path, || Ok(())).unwrap();
        }
        assert_eq!(reporter.started.load(Ordering::Relaxed), 2);
        assert_eq!(reporter.finished.load(Ordering::Relaxed), 2);
    }

    #[test]
    fn test_track_cancelled() {
        let reporter = CountingReporter::default();
        reporter.cancelled.store(true, Ordering::Relaxed);
        let mut ran = false;
        let result = track(Some(&reporter), Path::new("a.fastq"), || {
            ran = true;
            Ok(())
        });
        assert!(!ran);
        assert!(result.unwrap_err().downcast_ref::<Cancelled>().is_some());
        assert_eq!(reporter.started.load(Ordering::Relaxed), 0);
    }
}
