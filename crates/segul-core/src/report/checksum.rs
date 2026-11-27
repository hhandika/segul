//! Optional SHA-256 checksum configuration and streaming file hashing.

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use std::str::FromStr;

use anyhow::{Context, Result};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

/// Environment variable that sets the default [`ChecksumMode`] when none is
/// supplied on the CLI or API call.
pub const ENV_CHECKSUM_MODE: &str = "SEGUL_CHECKSUM";

const HASH_BUFFER_BYTES: usize = 64 * 1024;

/// Which files receive SHA-256 checksums in `report.json`.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ChecksumMode {
    /// Do not compute SHA-256 checksums (fastest).
    #[default]
    None,
    /// Compute SHA-256 for written output files only.
    Outputs,
    /// Compute SHA-256 for primary and auxiliary input files only.
    Inputs,
    /// Compute SHA-256 for both inputs and outputs.
    All,
}

impl ChecksumMode {
    /// Read the default mode from `SEGUL_CHECKSUM`, falling back to `None`
    /// if unset or unrecognized.
    pub fn from_env() -> Self {
        match std::env::var(ENV_CHECKSUM_MODE) {
            Ok(raw) => raw.trim().to_ascii_lowercase().parse().unwrap_or_default(),
            Err(_) => Self::None,
        }
    }

    /// Resolve an explicit optional mode against `SEGUL_CHECKSUM` and `None`.
    pub fn resolve(explicit: Option<Self>) -> Self {
        explicit.unwrap_or_else(Self::from_env)
    }

    /// Whether input files (primary and auxiliary) should be hashed.
    pub fn includes_inputs(self) -> bool {
        matches!(self, Self::Inputs | Self::All)
    }

    /// Whether output files should be hashed.
    pub fn includes_outputs(self) -> bool {
        matches!(self, Self::Outputs | Self::All)
    }

    /// Name of the hash algorithm recorded in `report.json`, or `None` when
    /// checksums are disabled.
    pub fn algorithm_name(self) -> Option<&'static str> {
        match self {
            Self::None => None,
            Self::Outputs | Self::Inputs | Self::All => Some("sha256"),
        }
    }
}

impl std::fmt::Display for ChecksumMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::None => write!(f, "none"),
            Self::Outputs => write!(f, "outputs"),
            Self::Inputs => write!(f, "inputs"),
            Self::All => write!(f, "all"),
        }
    }
}

impl FromStr for ChecksumMode {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "none" => Ok(Self::None),
            "outputs" => Ok(Self::Outputs),
            "inputs" => Ok(Self::Inputs),
            "all" => Ok(Self::All),
            _ => Err(format!(
                "'{}' is not a valid checksum mode. Valid options: none, outputs, inputs, all",
                s
            )),
        }
    }
}

/// Compute lowercase hex SHA-256 of a file using a 64 KiB streaming buffer.
pub fn sha256_file(path: &Path) -> Result<String> {
    let file = File::open(path).with_context(|| format!("Failed opening {}", path.display()))?;
    let mut reader = BufReader::with_capacity(HASH_BUFFER_BYTES, file);
    let mut hasher = Sha256::new();
    let mut buffer = [0u8; HASH_BUFFER_BYTES];
    loop {
        let read = reader
            .read(&mut buffer)
            .with_context(|| format!("Failed reading {}", path.display()))?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }
    let digest = hasher.finalize();
    Ok(hex_lower(&digest))
}

/// Compute SHA-256 digests for `paths` in parallel, preserving order.
pub fn sha256_files(paths: &[PathBuf]) -> Result<Vec<String>> {
    paths.par_iter().map(|path| sha256_file(path)).collect()
}

fn hex_lower(bytes: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut out = String::with_capacity(bytes.len() * 2);
    for &b in bytes {
        out.push(HEX[(b >> 4) as usize] as char);
        out.push(HEX[(b & 0x0f) as usize] as char);
    }
    out
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use tempfile::NamedTempFile;

    use super::*;

    #[test]
    fn test_checksum_mode_roundtrip() {
        for mode in [
            ChecksumMode::None,
            ChecksumMode::Outputs,
            ChecksumMode::Inputs,
            ChecksumMode::All,
        ] {
            let parsed: ChecksumMode = mode.to_string().parse().unwrap();
            assert_eq!(parsed, mode);
        }
        assert!("invalid".parse::<ChecksumMode>().is_err());
    }

    #[test]
    fn test_sha256_file() {
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(b"abc").unwrap();
        let digest = sha256_file(file.path()).unwrap();
        assert_eq!(
            digest,
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
    }
}
