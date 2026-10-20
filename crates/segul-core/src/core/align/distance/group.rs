//! Assign samples to groups, such as populations or species,
//! and accumulate the mean pairwise distance within and between groups.
use std::collections::{HashMap, HashSet};
use std::ffi::OsStr;
use std::path::Path;

use anyhow::{Context, Result, bail};
use regex::Regex;

/// Delimiters tried by [`SampleGrouping::detect`], in order of preference.
pub const DETECT_DELIMITERS: [char; 5] = ['_', '-', '.', '|', ' '];

/// Minimum proportion of samples that must share a group with another sample
/// for a detected grouping to be accepted.
const MIN_POOLED_PROPORTION: f64 = 0.8;

/// How to assign each sample to a group.
#[derive(Debug, Clone)]
pub enum SampleGrouping {
    /// Split the sample ID on `delimiter` and keep the first `fields` parts.
    /// For example, `Rattus_rattus_FMNH123` with `_` and 2 fields
    /// belongs to `Rattus_rattus`.
    Delimiter { delimiter: char, fields: usize },
    /// Use the first capture group of the pattern,
    /// or the whole match if the pattern has no capture group.
    Regex(Regex),
    /// Look up the group of each sample in a sample-to-group map.
    Map(HashMap<String, String>),
}

impl SampleGrouping {
    /// Group samples by the first `fields` parts of their IDs split on `delimiter`.
    pub fn from_delimiter(delimiter: char, fields: usize) -> Result<Self> {
        if fields == 0 {
            bail!("The number of fields to keep must be at least 1");
        }
        Ok(Self::Delimiter { delimiter, fields })
    }

    /// Group samples by a regular expression.
    pub fn from_regex(pattern: &str) -> Result<Self> {
        let re = Regex::new(pattern)
            .with_context(|| format!("Invalid group regular expression: {}", pattern))?;
        Ok(Self::Regex(re))
    }

    /// Group samples using a CSV or TSV file with two columns: sample and group.
    /// The first row is a header and is skipped.
    pub fn from_map_file(path: &Path) -> Result<Self> {
        let delimiter = match path.extension().and_then(OsStr::to_str) {
            Some(ext) if ext.eq_ignore_ascii_case("csv") => b',',
            Some(ext) if ext.eq_ignore_ascii_case("tsv") => b'\t',
            _ => bail!(
                "The group file {} must be a .csv or .tsv file",
                path.display()
            ),
        };
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(delimiter)
            .trim(csv::Trim::All)
            .from_path(path)
            .with_context(|| format!("Failed to open group file {}", path.display()))?;
        let mut map = HashMap::new();
        for (line, record) in reader.records().enumerate() {
            let record =
                record.with_context(|| format!("Failed to read group file {}", path.display()))?;
            match (record.get(0), record.get(1)) {
                (Some(sample), Some(group)) if !sample.is_empty() && !group.is_empty() => {
                    map.insert(sample.to_string(), group.to_string());
                }
                _ => bail!(
                    "Line {} of {} needs a sample and a group",
                    line + 2,
                    path.display()
                ),
            }
        }
        Ok(Self::Map(map))
    }

    /// Suggest a delimiter and field count from the sample IDs.
    ///
    /// Each delimiter in [`DETECT_DELIMITERS`] is tried with every field count
    /// that leaves at least one part out of every ID. A candidate is accepted
    /// when it gives more than one group, fewer groups than samples, and at least
    /// 80% of samples share their group with another sample. Among accepted
    /// candidates, the one with the most groups wins, so `Genus_species_ID`
    /// names are grouped by species rather than genus. Returns `None` if no
    /// candidate is accepted.
    pub fn detect(ids: &[String]) -> Option<GroupingSuggestion> {
        let mut best: Option<GroupingSuggestion> = None;
        for delimiter in DETECT_DELIMITERS {
            let min_parts = ids.iter().map(|id| id.split(delimiter).count()).min()?;
            for fields in 1..min_parts {
                let grouping = Self::Delimiter { delimiter, fields };
                let counts = grouping.assign(ids).counts();
                let groups = counts.len();
                let pooled: usize = counts.iter().map(|(_, n)| *n).filter(|n| *n > 1).sum();
                let is_accepted = groups > 1
                    && groups < ids.len()
                    && pooled as f64 / ids.len() as f64 >= MIN_POOLED_PROPORTION;
                let is_better = best.as_ref().is_none_or(|b| groups > b.groups);
                if is_accepted && is_better {
                    best = Some(GroupingSuggestion {
                        delimiter,
                        fields,
                        groups,
                    });
                }
            }
        }
        best
    }

    /// The group of a sample, or `None` if the sample cannot be assigned.
    pub fn group_of(&self, id: &str) -> Option<String> {
        let group = match self {
            Self::Delimiter { delimiter, fields } => {
                let parts: Vec<&str> = id.split(*delimiter).collect();
                if parts.len() < *fields {
                    return None;
                }
                parts[..*fields].join(&delimiter.to_string())
            }
            Self::Regex(re) => {
                let caps = re.captures(id)?;
                caps.get(1).or_else(|| caps.get(0))?.as_str().to_string()
            }
            Self::Map(map) => map.get(id)?.to_string(),
        };
        if group.is_empty() { None } else { Some(group) }
    }

    /// Assign each sample to a group, keeping the order of `ids`.
    pub fn assign(&self, ids: &[String]) -> GroupAssignment {
        GroupAssignment {
            ids: ids.to_vec(),
            groups: ids.iter().map(|id| self.group_of(id)).collect(),
        }
    }
}

/// A delimiter and field count suggested by [`SampleGrouping::detect`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GroupingSuggestion {
    pub delimiter: char,
    pub fields: usize,
    /// Number of groups the suggestion creates.
    pub groups: usize,
}

impl GroupingSuggestion {
    pub fn to_grouping(self) -> SampleGrouping {
        SampleGrouping::Delimiter {
            delimiter: self.delimiter,
            fields: self.fields,
        }
    }
}

/// Groups assigned to samples. `groups[i]` is the group of `ids[i]`.
#[derive(Debug, Clone)]
pub struct GroupAssignment {
    pub ids: Vec<String>,
    pub groups: Vec<Option<String>>,
}

impl GroupAssignment {
    /// Each group with its sample count, sorted by group name.
    pub fn counts(&self) -> Vec<(String, usize)> {
        let mut counts: HashMap<&str, usize> = HashMap::new();
        self.groups
            .iter()
            .flatten()
            .for_each(|g| *counts.entry(g.as_str()).or_default() += 1);
        let mut counts: Vec<(String, usize)> = counts
            .into_iter()
            .map(|(g, n)| (g.to_string(), n))
            .collect();
        counts.sort_by(|a, b| alphanumeric_sort::compare_str(&a.0, &b.0));
        counts
    }

    /// Samples that belong to no group.
    pub fn unassigned(&self) -> Vec<String> {
        self.ids
            .iter()
            .zip(&self.groups)
            .filter(|(_, g)| g.is_none())
            .map(|(id, _)| id.to_string())
            .collect()
    }
}

/// Running statistics of the pairwise distances between two groups.
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct GroupStats {
    /// Number of sample pairs with a defined distance.
    pub n_pairs: usize,
    sum: f64,
    sum_sq: f64,
    pub min: Option<f64>,
    pub max: Option<f64>,
}

impl GroupStats {
    pub fn mean(&self) -> Option<f64> {
        if self.n_pairs == 0 {
            return None;
        }
        Some(self.sum / self.n_pairs as f64)
    }

    /// Sample standard deviation. `None` if there are fewer than two pairs.
    pub fn sd(&self) -> Option<f64> {
        if self.n_pairs < 2 {
            return None;
        }
        let n = self.n_pairs as f64;
        let variance = (self.sum_sq - self.sum * self.sum / n) / (n - 1.0);
        Some(variance.max(0.0).sqrt())
    }

    fn add(&mut self, distance: f64) {
        self.n_pairs += 1;
        self.sum += distance;
        self.sum_sq += distance * distance;
        self.min = Some(self.min.map_or(distance, |m| m.min(distance)));
        self.max = Some(self.max.map_or(distance, |m| m.max(distance)));
    }

    fn merge(&mut self, other: &Self) {
        self.n_pairs += other.n_pairs;
        self.sum += other.sum;
        self.sum_sq += other.sum_sq;
        self.min = combine_option(self.min, other.min, f64::min);
        self.max = combine_option(self.max, other.max, f64::max);
    }
}

fn combine_option(a: Option<f64>, b: Option<f64>, pick: fn(f64, f64) -> f64) -> Option<f64> {
    match (a, b) {
        (Some(a), Some(b)) => Some(pick(a, b)),
        (a, b) => a.or(b),
    }
}

/// Distance statistics for every pair of groups.
#[derive(Debug, Default, Clone)]
pub struct GroupDistances {
    /// Groups with at least one sample, including groups without any pair.
    names: HashSet<String>,
    /// Statistics keyed by the two group names in sorted order.
    stats: HashMap<(String, String), GroupStats>,
}

impl GroupDistances {
    /// Group names sorted alphanumerically.
    pub fn names(&self) -> Vec<&str> {
        let mut names: Vec<&str> = self.names.iter().map(|n| n.as_str()).collect();
        names.sort_by(|a, b| alphanumeric_sort::compare_str(a, b));
        names
    }

    /// Statistics between groups `a` and `b`, in either order.
    pub fn get(&self, a: &str, b: &str) -> GroupStats {
        self.stats
            .get(&Self::key(a, b))
            .copied()
            .unwrap_or_default()
    }

    pub(crate) fn add_group(&mut self, group: &str) {
        if !self.names.contains(group) {
            self.names.insert(group.to_string());
        }
    }

    pub(crate) fn add(&mut self, a: &str, b: &str, distance: f64) {
        self.stats.entry(Self::key(a, b)).or_default().add(distance);
    }

    pub(crate) fn merge(&mut self, other: &Self) {
        self.names.extend(other.names.iter().cloned());
        other.stats.iter().for_each(|(key, stats)| {
            self.stats.entry(key.clone()).or_default().merge(stats);
        });
    }

    fn key(a: &str, b: &str) -> (String, String) {
        if a <= b {
            (a.to_string(), b.to_string())
        } else {
            (b.to_string(), a.to_string())
        }
    }
}

#[cfg(test)]
mod tests {
    use std::fs;

    use assert_approx_eq::assert_approx_eq;
    use tempfile::TempDir;

    use super::*;

    fn ids(names: &[&str]) -> Vec<String> {
        names.iter().map(|n| n.to_string()).collect()
    }

    #[test]
    fn test_delimiter_grouping() {
        let grouping = SampleGrouping::from_delimiter('_', 2).unwrap();
        assert_eq!(
            grouping.group_of("Rattus_rattus_FMNH123"),
            Some(String::from("Rattus_rattus"))
        );
        assert_eq!(
            grouping.group_of("Rattus_rattus"),
            Some(String::from("Rattus_rattus"))
        );
        assert_eq!(grouping.group_of("Rattus"), None);
        assert!(SampleGrouping::from_delimiter('_', 0).is_err());
    }

    #[test]
    fn test_regex_grouping() {
        let grouping = SampleGrouping::from_regex(r"^(\w+?)_\d+$").unwrap();
        assert_eq!(grouping.group_of("popA_12"), Some(String::from("popA")));
        assert_eq!(grouping.group_of("popA-12"), None);
        let grouping = SampleGrouping::from_regex(r"^[A-Z]+").unwrap();
        assert_eq!(grouping.group_of("ABC123"), Some(String::from("ABC")));
        assert!(SampleGrouping::from_regex("(").is_err());
    }

    #[test]
    fn test_map_file_grouping() {
        let dir = TempDir::with_prefix("temp").unwrap();
        let path = dir.path().join("groups.tsv");
        fs::write(&path, "sample\tgroup\nA1\tNorth pop\nA2\tSouth\n").unwrap();
        let grouping = SampleGrouping::from_map_file(&path).unwrap();
        assert_eq!(grouping.group_of("A1"), Some(String::from("North pop")));
        assert_eq!(grouping.group_of("A3"), None);

        let path = dir.path().join("groups.csv");
        fs::write(&path, "sample,group\nA1,\n").unwrap();
        assert!(SampleGrouping::from_map_file(&path).is_err());
        assert!(SampleGrouping::from_map_file(&dir.path().join("groups.txt")).is_err());
    }

    #[test]
    fn test_assignment() {
        let grouping = SampleGrouping::from_regex(r"^(sp\d)_").unwrap();
        let assignment = grouping.assign(&ids(&["sp1_a", "sp2_a", "sp1_b", "other"]));
        assert_eq!(
            assignment.counts(),
            vec![(String::from("sp1"), 2), (String::from("sp2"), 1)]
        );
        assert_eq!(assignment.unassigned(), vec![String::from("other")]);
    }

    #[test]
    fn test_detect_species() {
        let names = ids(&[
            "Rattus_rattus_A1",
            "Rattus_rattus_A2",
            "Rattus_exulans_B1",
            "Rattus_exulans_B2",
            "Mus_musculus_C1",
            "Mus_musculus_C2",
        ]);
        let suggestion = SampleGrouping::detect(&names).unwrap();
        assert_eq!(suggestion.delimiter, '_');
        assert_eq!(suggestion.fields, 2);
        assert_eq!(suggestion.groups, 3);
    }

    #[test]
    fn test_detect_other_delimiter() {
        let names = ids(&["popA-1", "popA-2", "popB-1", "popB-2"]);
        let suggestion = SampleGrouping::detect(&names).unwrap();
        assert_eq!((suggestion.delimiter, suggestion.fields), ('-', 1));
    }

    #[test]
    fn test_detect_none() {
        assert_eq!(
            SampleGrouping::detect(&ids(&["ABCD", "ABCE", "ABCF"])),
            None
        );
        assert_eq!(SampleGrouping::detect(&ids(&["a_1", "b_2", "c_3"])), None);
        assert_eq!(SampleGrouping::detect(&[]), None);
    }

    #[test]
    fn test_group_stats() {
        let mut distances = GroupDistances::default();
        distances.add_group("b");
        distances.add_group("a");
        distances.add("b", "a", 0.1);
        distances.add("a", "b", 0.3);
        let stats = distances.get("a", "b");
        assert_eq!(stats.n_pairs, 2);
        assert_approx_eq!(stats.mean().unwrap(), 0.2);
        assert_approx_eq!(stats.sd().unwrap(), 0.02_f64.sqrt());
        assert_eq!(stats.min, Some(0.1));
        assert_eq!(distances.names(), vec!["a", "b"]);
        assert_eq!(distances.get("a", "a").mean(), None);

        let mut other = GroupDistances::default();
        other.add_group("c");
        other.add("a", "b", 0.5);
        distances.merge(&other);
        let stats = distances.get("b", "a");
        assert_eq!(stats.n_pairs, 3);
        assert_approx_eq!(stats.mean().unwrap(), 0.3);
        assert_eq!(stats.max, Some(0.5));
        assert_eq!(distances.names(), vec!["a", "b", "c"]);
    }
}
