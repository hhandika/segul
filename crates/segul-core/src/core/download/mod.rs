//! Concurrent GenBank and ENA archive downloader.

pub mod assembly;
pub mod client;
pub mod fasta;
pub mod job;
pub mod limiter;
pub mod types;

pub use client::{ArchiveClient, CancelSignal};
pub use job::DownloadJob;
pub use limiter::{ConcurrencyProfile, HostLimiter};
pub use types::{
    CompletedDownload, DatabaseSource, DownloadCancelled, DownloadConfig, DownloadContent,
    DownloadError, DownloadFailure, DownloadOptions, DownloadProgress, DownloadProgressReporter,
    DownloadRecord, DownloadReport, DownloadStage, RecordDownloadStatus, RemoteFile,
    SequenceGrouping, SequenceNaming,
};
