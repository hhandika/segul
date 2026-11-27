//! Multi-step workflow recipes (`segul-workflow.json`) built from chained
//! `parent_reports` lineage or explicit report sequences.

use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use chrono::{SecondsFormat, Utc};
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::helper::utils;

use super::RunReport;
use super::lineage::canonical_or_normalized;
use super::paths::resolve_portable_path;
use super::spec::TaskKind;

/// Schema version for `segul-workflow.json`.
pub const WORKFLOW_SCHEMA_VERSION: u32 = 1;

/// Declarative multi-step pipeline recipe generated from chained runs.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WorkflowRecipe {
    pub schema_version: u32,
    pub name: String,
    pub created_utc: String,
    pub segul_core_version: String,
    pub steps: Vec<WorkflowStep>,
}

/// Source of primary inputs for a [`WorkflowStep`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "source", rename_all = "snake_case")]
pub enum WorkflowInputSource {
    /// External files resolved from `base_dir`.
    Files {
        base_dir: String,
        files: Vec<String>,
    },
    /// Outputs produced by an earlier step in the same workflow.
    FromStep {
        step_id: String,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        roles: Vec<String>,
    },
}

/// Auxiliary input file referenced by a [`WorkflowStep`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct WorkflowAuxInput {
    pub role: String,
    pub base_dir: String,
    pub path: String,
}

/// Single step in a [`WorkflowRecipe`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WorkflowStep {
    pub id: String,
    pub kind: TaskKind,
    pub inputs: WorkflowInputSource,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub auxiliary: Vec<WorkflowAuxInput>,
    pub parameters: Value,
    pub output_subdir: String,
}

impl WorkflowRecipe {
    /// Load a workflow recipe from `path`.
    pub fn from_path(path: &Path) -> Result<Self> {
        let text = fs::read_to_string(path)
            .with_context(|| format!("Failed reading {}", path.display()))?;
        serde_json::from_str(&text).with_context(|| format!("Failed parsing {}", path.display()))
    }

    /// Write formatted workflow JSON to `path`.
    pub fn write(&self, path: &Path) -> Result<()> {
        if let Some(parent) = path.parent()
            && !parent.as_os_str().is_empty()
        {
            fs::create_dir_all(parent)
                .with_context(|| format!("Failed creating {}", parent.display()))?;
        }
        let json = serde_json::to_string_pretty(self)?;
        fs::write(path, json).with_context(|| format!("Failed writing {}", path.display()))
    }

    /// Walk `parent_reports` transitively starting at `leaf_report` and
    /// construct an executable topological [`WorkflowRecipe`].
    pub fn from_report_chain(leaf_report: &Path, name: Option<&str>) -> Result<Self> {
        let reports = collect_chain_with_paths(leaf_report)?;
        let default_name = name.unwrap_or("segul-workflow");
        Ok(build_recipe_from_pairs(&reports, default_name))
    }

    /// Build a [`WorkflowRecipe`] from an ordered list of report files
    /// (for example, runs recorded in a SEGUI Workbench session).
    pub fn from_report_paths(report_paths: &[PathBuf], name: Option<&str>) -> Result<Self> {
        let mut pairs = Vec::with_capacity(report_paths.len());
        for path in report_paths {
            let report = RunReport::from_path(path)?;
            pairs.push((canonical_or_normalized(path), report));
        }
        Ok(build_recipe_from_pairs(
            &pairs,
            name.unwrap_or("segul-workflow"),
        ))
    }
}

/// Collect a leaf report and all transitive parent reports in chronological
/// topological order (roots first, leaf last).
pub(super) fn collect_chain_with_paths(leaf: &Path) -> Result<Vec<(PathBuf, RunReport)>> {
    let mut visited: HashSet<PathBuf> = HashSet::new();
    let mut ordered: Vec<(PathBuf, RunReport)> = Vec::new();
    visit_postorder(leaf, &mut visited, &mut ordered)?;
    Ok(ordered)
}

fn visit_postorder(
    report_path: &Path,
    visited: &mut HashSet<PathBuf>,
    ordered: &mut Vec<(PathBuf, RunReport)>,
) -> Result<()> {
    let canon = canonical_or_normalized(report_path);
    if !visited.insert(canon.clone()) {
        return Ok(());
    }
    let report = RunReport::from_path(report_path)?;
    let report_dir = report_path.parent().unwrap_or_else(|| Path::new("."));
    for parent in &report.run.parent_reports {
        let parent_path = resolve_portable_path(report_dir, &parent.path);
        if parent_path.is_file() {
            visit_postorder(&parent_path, visited, ordered)?;
        }
    }
    ordered.push((canon, report));
    Ok(())
}

fn build_recipe_from_pairs(pairs: &[(PathBuf, RunReport)], name: &str) -> WorkflowRecipe {
    let mut path_to_step: HashMap<PathBuf, String> = HashMap::new();
    let mut steps = Vec::with_capacity(pairs.len());

    for (idx, (report_path, report)) in pairs.iter().enumerate() {
        let step_id = format!(
            "step_{}_{}",
            idx + 1,
            report.kind.as_str().replace('.', "_")
        );
        let report_dir = report_path.parent().unwrap_or_else(|| Path::new("."));

        // Link to upstream step when all primary inputs match a single parent
        // report in this recipe.
        let mut linked_source = None;
        if let Some(parent_ref) = report.run.parent_reports.first()
            && report.run.parent_reports.len() == 1
        {
            let parent_abs =
                canonical_or_normalized(&resolve_portable_path(report_dir, &parent_ref.path));
            if let Some(upstream_id) = path_to_step.get(&parent_abs)
                && parent_ref.matched_files.len() == report.inputs.primary.len()
            {
                let parent_report = pairs.iter().find(|(p, _)| p == &parent_abs).map(|(_, r)| r);
                let roles = infer_output_roles(parent_report, parent_abs.parent(), report);
                linked_source = Some(WorkflowInputSource::FromStep {
                    step_id: upstream_id.clone(),
                    roles,
                });
            }
        }

        let inputs = linked_source.unwrap_or_else(|| WorkflowInputSource::Files {
            base_dir: report.inputs.base_dir.clone(),
            files: report
                .inputs
                .primary
                .iter()
                .map(|f| f.path.clone())
                .collect(),
        });

        let auxiliary = report
            .inputs
            .auxiliary
            .iter()
            .map(|a| WorkflowAuxInput {
                role: a.role.clone(),
                base_dir: report.inputs.base_dir.clone(),
                path: a.path.clone(),
            })
            .collect();

        path_to_step.insert(report_path.clone(), step_id.clone());
        steps.push(WorkflowStep {
            id: step_id,
            kind: report.kind,
            inputs,
            auxiliary,
            parameters: report.parameters.clone(),
            output_subdir: format!("{:02}_{}", idx + 1, report.kind.as_str().replace('.', "_")),
        });
    }

    WorkflowRecipe {
        schema_version: WORKFLOW_SCHEMA_VERSION,
        name: name.to_string(),
        created_utc: Utc::now().to_rfc3339_opts(SecondsFormat::Secs, true),
        segul_core_version: utils::get_crate_version(),
        steps,
    }
}

fn infer_output_roles(
    parent_report: Option<&RunReport>,
    parent_dir: Option<&Path>,
    child_report: &RunReport,
) -> Vec<String> {
    let (Some(parent), Some(pdir)) = (parent_report, parent_dir) else {
        return Vec::new();
    };
    let child_base = resolve_portable_path(Path::new("."), &child_report.inputs.base_dir);
    let child_set: HashSet<PathBuf> = child_report
        .inputs
        .primary
        .iter()
        .map(|f| canonical_or_normalized(&resolve_portable_path(&child_base, &f.path)))
        .collect();
    let mut roles: Vec<String> = parent
        .outputs
        .files
        .iter()
        .filter(|o| {
            let p = canonical_or_normalized(&resolve_portable_path(pdir, &o.path));
            child_set.contains(&p)
        })
        .map(|o| o.role.clone())
        .collect();
    roles.sort();
    roles.dedup();
    roles
}
