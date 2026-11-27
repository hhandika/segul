//! Backward-compatible deserialization for v2 `RunReport` envelopes and
//! legacy flat `*-report.json` files.

use anyhow::{Context, Result};
use serde::de::DeserializeOwned;
use serde_json::Value;

/// Deserialize a summary payload (`AlignReport`, `ContigReport`,
/// `DistanceReport`, `ReadReport`) from either a v2 `RunReport` JSON
/// (extracting `"statistics"`) or a legacy flat JSON report.
pub fn parse_report_payload<T: DeserializeOwned>(json: &str) -> Result<T> {
    let value: Value = serde_json::from_str(json).context("Invalid JSON report")?;
    if value.get("schema_version").is_some()
        && let Some(stats) = value.get("statistics")
        && !stats.is_null()
    {
        return serde_json::from_value(stats.clone())
            .context("Failed parsing report statistics payload");
    }
    serde_json::from_value(value).context("Failed parsing legacy report payload")
}
