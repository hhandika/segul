//! Helper functions for the SEGUL features.
pub mod alphabet;
#[cfg(feature = "align")]
pub mod concat;
pub mod files;
pub mod finder;
pub mod logger;
pub mod partition;
pub mod sequence;
#[cfg(feature = "sequence")]
pub mod translation;
pub mod types;
pub mod utils;
