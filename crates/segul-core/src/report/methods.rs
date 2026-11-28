//! Publication-ready Methods text and summary tables rendered from
//! `report.json`, chained `parent_reports`, or `segul-workflow.json`.

use std::fmt::Write;
use std::path::Path;
use std::str::FromStr;

use anyhow::Result;

use super::RunReport;
use super::spec::TaskParameters;
use super::workflow::{WorkflowRecipe, collect_chain_with_paths};

/// Output format for generated Methods text.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum MethodsFormat {
    #[default]
    Markdown,
    PlainText,
}

impl std::fmt::Display for MethodsFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Markdown => write!(f, "markdown"),
            Self::PlainText => write!(f, "text"),
        }
    }
}

impl FromStr for MethodsFormat {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s.trim().to_ascii_lowercase().as_str() {
            "markdown" | "md" => Ok(Self::Markdown),
            "text" | "txt" | "plain" => Ok(Self::PlainText),
            _ => Err(format!(
                "'{}' is not a valid methods format. Valid options: markdown, text",
                s
            )),
        }
    }
}

/// Collect a leaf `report.json` and all ancestor reports in execution order.
pub fn collect_report_chain(leaf_report: &Path) -> Result<Vec<RunReport>> {
    Ok(collect_chain_with_paths(leaf_report)?
        .into_iter()
        .map(|(_, report)| report)
        .collect())
}

/// Render a publication-ready Methods section from a `report.json` (optionally
/// including its full `parent_reports` chain) or a `segul-workflow.json`.
pub fn render_methods_from_path(path: &Path, chain: bool, format: MethodsFormat) -> Result<String> {
    if let Ok(report) = RunReport::from_path(path) {
        let reports = if chain {
            collect_report_chain(path)?
        } else {
            vec![report]
        };
        return Ok(render_methods(&reports, format));
    }
    let workflow = WorkflowRecipe::from_path(path)?;
    Ok(render_workflow_methods(&workflow, format))
}

/// Render Methods prose and a summary parameter table from an ordered slice of
/// [`RunReport`]s.
pub fn render_methods(reports: &[RunReport], format: MethodsFormat) -> String {
    if reports.is_empty() {
        return String::new();
    }
    let mut out = String::new();
    let version = &reports
        .last()
        .map(|r| r.run.segul_core_version.as_str())
        .unwrap_or("unknown");

    let frontends: Vec<String> = {
        let mut names = Vec::new();
        for r in reports {
            let label = format!("{} v{}", r.run.frontend.name, r.run.frontend.version);
            if !names.contains(&label) {
                names.push(label);
            }
        }
        names
    };

    match format {
        MethodsFormat::Markdown => {
            let _ = writeln!(out, "## Methods\n");
            let _ = writeln!(
                out,
                "Analyses were performed using SEGUL core v{} (via {}).\n",
                version,
                frontends.join(", ")
            );
            for (idx, report) in reports.iter().enumerate() {
                let sentence = describe_report_step(report);
                let _ = writeln!(
                    out,
                    "{}. **{}**: {}",
                    idx + 1,
                    report.kind.label(),
                    sentence
                );
            }
            let _ = writeln!(out, "\n### Reproducibility Parameters\n");
            let _ = writeln!(
                out,
                "| Step | Task | Parameters | Inputs | Outputs | SHA-256 |"
            );
            let _ = writeln!(out, "| --- | --- | --- | --- | --- | --- |");
            for (idx, report) in reports.iter().enumerate() {
                let _ = writeln!(
                    out,
                    "| {} | `{}` | {} | {} | {} | {} |",
                    idx + 1,
                    report.kind,
                    summarize_params(report),
                    report.inputs.primary.len(),
                    report.outputs.files.len(),
                    report.run.checksum_mode
                );
            }
        }
        MethodsFormat::PlainText => {
            let _ = writeln!(out, "Methods");
            let _ = writeln!(
                out,
                "Analyses were performed using SEGUL core v{} (via {}).\n",
                version,
                frontends.join(", ")
            );
            for (idx, report) in reports.iter().enumerate() {
                let _ = writeln!(
                    out,
                    "{}. {}: {}",
                    idx + 1,
                    report.kind.label(),
                    describe_report_step(report)
                );
            }
            let _ = writeln!(out, "\nParameters Summary:");
            for (idx, report) in reports.iter().enumerate() {
                let _ = writeln!(
                    out,
                    "- Step {} [{}]: {} (inputs: {}, outputs: {}, checksums: {})",
                    idx + 1,
                    report.kind,
                    summarize_params(report),
                    report.inputs.primary.len(),
                    report.outputs.files.len(),
                    report.run.checksum_mode
                );
            }
        }
    }
    out
}

fn render_workflow_methods(workflow: &WorkflowRecipe, format: MethodsFormat) -> String {
    let mut out = String::new();
    match format {
        MethodsFormat::Markdown => {
            let _ = writeln!(out, "## Methods ({})\n", workflow.name);
            let _ = writeln!(
                out,
                "A {}-step reproducible SEGUL workflow (segul-core v{}) was configured:\n",
                workflow.steps.len(),
                workflow.segul_core_version
            );
            for (idx, step) in workflow.steps.iter().enumerate() {
                let desc = describe_value_params(step.kind, &step.parameters);
                let _ = writeln!(out, "{}. **{}** ({})", idx + 1, step.kind.label(), desc);
            }
        }
        MethodsFormat::PlainText => {
            let _ = writeln!(out, "Methods ({})", workflow.name);
            let _ = writeln!(
                out,
                "A {}-step SEGUL workflow (segul-core v{}):",
                workflow.steps.len(),
                workflow.segul_core_version
            );
            for (idx, step) in workflow.steps.iter().enumerate() {
                let desc = describe_value_params(step.kind, &step.parameters);
                let _ = writeln!(out, "{}. {} ({})", idx + 1, step.kind.label(), desc);
            }
        }
    }
    out
}

fn describe_report_step(report: &RunReport) -> String {
    let param_summary = summarize_params(report);
    format!(
        "Processed {} primary input file(s) producing {} output file(s) ({}; checksum mode: {}).",
        report.inputs.primary.len(),
        report.outputs.files.len(),
        param_summary,
        report.run.checksum_mode
    )
}

fn summarize_params(report: &RunReport) -> String {
    describe_value_params(report.kind, &report.parameters)
}

fn describe_value_params(kind: super::spec::TaskKind, val: &serde_json::Value) -> String {
    let Ok(parsed) = TaskParameters::from_json_value(kind, val) else {
        return val.to_string();
    };
    match parsed {
        TaskParameters::AlignSummary(p) => format!(
            "format={}, datatype={}, interval={}%, per_locus={}",
            p.input_format, p.datatype, p.interval, p.per_locus
        ),
        TaskParameters::AlignDistance(p) => format!(
            "format={}, datatype={}, model={}, triangle={}, long={}",
            p.input_format, p.datatype, p.model, p.triangle, p.long_format
        ),
        TaskParameters::AlignConcat(p) => format!(
            "in={}, out={}, partition={}, datatype={}, sort={}",
            p.input_format, p.output_format, p.partition_format, p.datatype, p.sort
        ),
        TaskParameters::AlignConvert(p) => format!(
            "in={}, out={}, datatype={}, sort={}",
            p.input_format, p.output_format, p.datatype, p.sort
        ),
        TaskParameters::AlignFilter(p) => format!(
            "in={}, datatype={}, criterion={:?}, concat={}",
            p.input_format, p.datatype, p.criterion, p.concat
        ),
        TaskParameters::AlignSplit(p) => format!(
            "in={}, out={}, partition={}, datatype={}",
            p.input_format, p.output_format, p.partition_format, p.datatype
        ),
        TaskParameters::AlignTrim(p) => format!(
            "in={}, out={}, datatype={}, criterion={:?}",
            p.input_format, p.output_format, p.datatype, p.criterion
        ),
        TaskParameters::AlignUnalign(p) => format!(
            "in={}, out={}, datatype={}",
            p.input_format, p.output_format, p.datatype
        ),
        TaskParameters::PartitionConvert(p) => format!(
            "in={}, out={}, datatype={}",
            p.input_partition_format, p.output_partition_format, p.datatype
        ),
        TaskParameters::ContigSummary(p) => format!("format={}", p.input_format),
        TaskParameters::ReadSummary(p) => format!(
            "format={}, mode={}, phred={}, qc_profile={}, dup_sketch_size={}",
            p.input_format, p.mode, p.phred, p.qc_profile, p.dup_sketch_size
        ),
        TaskParameters::MafSummary(p) => format!(
            "naming={}, ref_species={}",
            p.naming,
            p.resolved_ref_species
                .as_deref()
                .or(p.ref_species.as_deref())
                .unwrap_or("auto")
        ),
        TaskParameters::MafConvert(p) => format!(
            "out={}, group_by={}, clip={}, strand_aware={}, ref_species={}",
            p.output_format,
            p.group_by,
            p.clip,
            p.strand_aware,
            p.resolved_ref_species
                .as_deref()
                .or(p.ref_species.as_deref())
                .unwrap_or("auto")
        ),
        TaskParameters::SequenceAdd(p) => format!(
            "in={}, dest={}, out={}, added_only={}, single_file={}",
            p.input_format, p.dest_format, p.output_format, p.added_only, p.single_file
        ),
        TaskParameters::SequenceExtract(p) => format!(
            "in={}, out={}, datatype={}, mode={}",
            p.input_format, p.output_format, p.datatype, p.mode
        ),
        TaskParameters::SequenceFilter(p) => format!(
            "in={}, out={}, datatype={}, criterion={:?}",
            p.input_format, p.output_format, p.datatype, p.criterion
        ),
        TaskParameters::SequenceId(p) => {
            format!(
                "in={}, datatype={}, map={}",
                p.input_format, p.datatype, p.map_id
            )
        }
        TaskParameters::SequenceRemove(p) => format!(
            "in={}, out={}, datatype={}, mode={}",
            p.input_format, p.output_format, p.datatype, p.mode
        ),
        TaskParameters::SequenceRename(p) => format!(
            "in={}, out={}, mode={}, dry_run={}",
            p.input_format, p.output_format, p.mode, p.dry_run
        ),
        TaskParameters::SequenceTranslate(p) => format!(
            "in={}, out={}, table={}, frame={}, auto={}",
            p.input_format, p.output_format, p.table, p.reading_frame, p.auto_frame
        ),
        TaskParameters::SequenceTranspose(p) => {
            format!(
                "in={}, out={}, datatype={}",
                p.input_format, p.output_format, p.datatype
            )
        }
        TaskParameters::TreeConvert(p) => format!(
            "in={}, out={}, translate={}, strip_annotations={}, gzip={}",
            p.input_format, p.output_format, p.translate, p.strip_annotations, p.gzip
        ),
    }
}
