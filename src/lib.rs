//! # SEGUL command-line interface
//! SEGUL provides fast and memory efficient utilities for alignment manipulation and summary statistics.
//!
//! ## Library users
//! The library code now lives in the [`segul-core`](https://crates.io/crates/segul-core) crate,
//! which lets you enable only the features you need:
//!
//! ```toml
//! segul-core = { version = "0.24", default-features = false, features = ["align"] }
//! ```
//!
//! Replace `segul::` with `segul_core::` in your imports. The modules below are
//! re-exported from `segul-core` for backward compatibility and will be removed
//! in a future release.
pub use segul_core::{core, helper, parser, stats, writer};

pub mod cli;
