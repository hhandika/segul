//! Cross-frontend replay for single `report.json` files and multi-step
//! `segul-workflow.json` recipes.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};

use super::checksum::ChecksumMode;
use super::paths::{resolve_base_dir, resolve_portable_path};
use super::spec::{TaskKind, TaskParameters};
use super::verify::{VerifyScope, verify_report};
use super::workflow::{WorkflowInputSource, WorkflowRecipe};
use super::{FrontendInfo, ReportContext, RunReport, SourceMode};

/// Options controlling how a `report.json` or `segul-workflow.json` is replayed.
#[derive(Debug, Clone, Default)]
pub struct ReplayOptions {
    /// Override the base directory used to resolve recorded input paths.
    pub input_base: Option<PathBuf>,
    /// Override the output directory.
    pub output_dir: Option<PathBuf>,
    /// Verify recorded input SHA-256 / sizes before executing.
    pub verify_inputs: bool,
    /// Override checksum mode for the newly generated `report.json`.
    pub checksum_mode: Option<ChecksumMode>,
    /// Frontend identity recorded in the new `report.json`.
    pub frontend: FrontendInfo,
}

/// Replay a single `report.json`, returning the newly written [`RunReport`].
pub fn replay_report(report_path: &Path, opts: &ReplayOptions) -> Result<RunReport> {
    let report = RunReport::from_path(report_path)?;
    if opts.verify_inputs {
        let verification = verify_report(
            report_path,
            opts.input_base.as_deref(),
            None,
            VerifyScope::Inputs,
        )?;
        if !verification.passed {
            bail!(
                "Input verification failed for {} ({} missing, {} mismatched)",
                report_path.display(),
                verification.missing,
                verification.mismatched
            );
        }
    }

    let base_dir = match &opts.input_base {
        Some(dir) => dir.clone(),
        None => resolve_base_dir(&report.inputs.base_dir, None),
    };
    let output_dir = match &opts.output_dir {
        Some(dir) => dir.clone(),
        None => resolve_portable_path(Path::new("."), &report.outputs.output_dir),
    };
    let primary: Vec<PathBuf> = report
        .inputs
        .primary
        .iter()
        .map(|f| resolve_portable_path(&base_dir, &f.path))
        .collect();
    let auxiliary: Vec<(String, PathBuf)> = report
        .inputs
        .auxiliary
        .iter()
        .map(|a| (a.role.clone(), resolve_portable_path(&base_dir, &a.path)))
        .collect();
    let params = report.typed_parameters()?;
    let checksum_mode = opts.checksum_mode.unwrap_or(report.run.checksum_mode);
    let ctx = ReportContext {
        frontend: opts.frontend.clone(),
        checksum_mode,
        source_mode: SourceMode::Files,
        input_dir: None,
    };

    execute_task(
        report.kind,
        &params,
        &primary,
        &auxiliary,
        &output_dir,
        &ctx,
    )
}

/// Replay every step of a `segul-workflow.json` sequentially, chaining outputs
/// from earlier steps into later steps.
pub fn replay_workflow(workflow_path: &Path, opts: &ReplayOptions) -> Result<Vec<RunReport>> {
    let recipe = WorkflowRecipe::from_path(workflow_path)?;
    let root_out = opts
        .output_dir
        .clone()
        .unwrap_or_else(|| PathBuf::from(&recipe.name));
    fs::create_dir_all(&root_out)
        .with_context(|| format!("Failed creating {}", root_out.display()))?;

    let checksum_mode = ChecksumMode::resolve(opts.checksum_mode);
    let ctx = ReportContext {
        frontend: opts.frontend.clone(),
        checksum_mode,
        source_mode: SourceMode::Files,
        input_dir: None,
    };

    let mut step_outputs: HashMap<String, (PathBuf, RunReport)> = HashMap::new();
    let mut reports = Vec::with_capacity(recipe.steps.len());

    for step in &recipe.steps {
        let step_out_dir = root_out.join(&step.output_subdir);
        let primary = match &step.inputs {
            WorkflowInputSource::Files { base_dir, files } => {
                let base = match &opts.input_base {
                    Some(dir) => dir.clone(),
                    None => resolve_base_dir(base_dir, None),
                };
                files
                    .iter()
                    .map(|f| resolve_portable_path(&base, f))
                    .collect::<Vec<_>>()
            }
            WorkflowInputSource::FromStep { step_id, roles } => {
                let (prev_dir, prev_report) = step_outputs.get(step_id).with_context(|| {
                    format!(
                        "Workflow step '{}' references missing step '{}'",
                        step.id, step_id
                    )
                })?;
                let mut selected: Vec<PathBuf> = prev_report
                    .outputs
                    .files
                    .iter()
                    .filter(|o| roles.is_empty() || roles.iter().any(|r| r == &o.role))
                    .map(|o| resolve_portable_path(prev_dir, &o.path))
                    .collect();
                alphanumeric_sort::sort_path_slice(&mut selected);
                if selected.is_empty() {
                    bail!(
                        "Workflow step '{}' produced no matching output files for '{}'",
                        step_id,
                        step.id
                    );
                }
                selected
            }
        };

        let auxiliary: Vec<(String, PathBuf)> = step
            .auxiliary
            .iter()
            .map(|a| {
                let base = match &opts.input_base {
                    Some(dir) => dir.clone(),
                    None => resolve_base_dir(&a.base_dir, None),
                };
                (a.role.clone(), resolve_portable_path(&base, &a.path))
            })
            .collect();

        let params = TaskParameters::from_json_value(step.kind, &step.parameters)?;
        let report = execute_task(
            step.kind,
            &params,
            &primary,
            &auxiliary,
            &step_out_dir,
            &ctx,
        )?;
        step_outputs.insert(step.id.clone(), (step_out_dir, report.clone()));
        reports.push(report);
    }

    Ok(reports)
}

#[allow(unused_variables)]
fn execute_task(
    kind: TaskKind,
    params: &TaskParameters,
    primary: &[PathBuf],
    auxiliary: &[(String, PathBuf)],
    output_dir: &Path,
    ctx: &ReportContext,
) -> Result<RunReport> {
    match (kind, params) {
        #[cfg(feature = "align")]
        (TaskKind::AlignSummary, TaskParameters::AlignSummary(p)) => {
            use crate::core::align::summarize::AlignmentSummary;
            use crate::helper::types::{DataType, InputFmt};
            let input_fmt: InputFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            let mut summary = AlignmentSummary::new(&input_fmt, output_dir, p.interval, &datatype)
                .with_report_context(ctx.clone());
            if p.per_locus {
                let files = summary.summarize_locus_with_report(primary, p.prefix.as_deref())?;
                Ok(files)
            } else {
                let files = summary.summarize_all(primary, p.prefix.as_deref());
                RunReport::from_path(&files.report)
            }
        }
        #[cfg(feature = "align")]
        (TaskKind::AlignDistance, TaskParameters::AlignDistance(p)) => {
            use super::spec::DistanceGroupingSpec;
            use crate::core::align::distance::group::SampleGrouping;
            use crate::core::align::distance::{AlignmentDistance, MatrixDiagonal, MatrixTriangle};
            use crate::helper::finder::IDs;
            use crate::helper::types::{DataType, InputFmt};
            use crate::stats::distance::DistanceModel;

            let input_fmt: InputFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            let model: DistanceModel = p.model.parse().map_err(anyhow::Error::msg)?;
            let triangle: MatrixTriangle = p.triangle.parse().map_err(anyhow::Error::msg)?;
            let diagonal = p
                .diagonal
                .as_deref()
                .map(|d| d.parse::<MatrixDiagonal>().map_err(anyhow::Error::msg))
                .transpose()?;
            let mut handle = AlignmentDistance::new(&input_fmt, &datatype, &model, output_dir)
                .with_triangle(triangle)
                .with_long_format(p.long_format)
                .with_report_context(ctx.clone());
            if let Some(d) = diagonal {
                handle = handle.with_diagonal(d);
            }
            if let Some(group_spec) = p.resolved_grouping.as_ref().or(p.grouping.as_ref()) {
                let grouping = match group_spec {
                    DistanceGroupingSpec::Delimiter { delimiter, fields } => {
                        let ch = delimiter
                            .chars()
                            .next()
                            .context("Empty group delimiter in report")?;
                        SampleGrouping::from_delimiter(ch, *fields)?
                    }
                    DistanceGroupingSpec::Regex { pattern } => SampleGrouping::from_regex(pattern)?,
                    DistanceGroupingSpec::MapFile { .. } => {
                        let map_path = aux_path(auxiliary, "group_map")
                            .context("Missing group_map auxiliary input for distance replay")?;
                        handle = handle.with_group_map_path(map_path.to_path_buf());
                        SampleGrouping::from_map_file(map_path)?
                    }
                    DistanceGroupingSpec::Auto => {
                        let ids: Vec<String> = IDs::new(primary, &input_fmt, &datatype)
                            .id_unique()
                            .into_iter()
                            .collect();
                        let suggestion = SampleGrouping::detect(&ids)
                            .context("Could not auto-detect sample groups during replay")?;
                        suggestion.to_grouping()
                    }
                };
                handle = handle.with_grouping(grouping);
            }
            let files = handle.compute(primary, p.prefix.as_deref());
            RunReport::from_path(&files.report)
        }
        #[cfg(feature = "align")]
        (TaskKind::AlignConcat, TaskParameters::AlignConcat(p)) => {
            use crate::core::align::concat::AlignmentConcatenation;
            use crate::helper::types::{DataType, InputFmt, OutputFmt, PartitionFmt};
            let input_fmt: InputFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            let output_fmt: OutputFmt = p.output_format.parse().map_err(anyhow::Error::msg)?;
            let part_fmt: PartitionFmt = p.partition_format.parse().map_err(anyhow::Error::msg)?;
            let prefix_buf = p
                .prefix
                .as_deref()
                .map(PathBuf::from)
                .unwrap_or_else(|| output_dir.to_path_buf());
            let mut files = primary.to_vec();
            let mut concat = AlignmentConcatenation::new(
                &input_fmt,
                output_dir,
                &output_fmt,
                &part_fmt,
                &prefix_buf,
            )
            .with_sort(p.sort)
            .with_report_context(ctx.clone());
            concat.concat_with_report(&mut files, &datatype, p.prefix.as_deref())
        }
        #[cfg(feature = "align")]
        (TaskKind::AlignConvert, TaskParameters::AlignConvert(p)) => {
            use crate::core::align::convert::AlignmentConversion;
            use crate::helper::types::{DataType, InputFmt, OutputFmt};
            let input_fmt: InputFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            let output_fmt: OutputFmt = p.output_format.parse().map_err(anyhow::Error::msg)?;
            AlignmentConversion::new(&input_fmt, &output_fmt, &datatype, p.sort)
                .with_report_context(ctx.clone())
                .convert_with_report(primary, output_dir)
        }
        #[cfg(feature = "align")]
        (TaskKind::AlignFilter, TaskParameters::AlignFilter(p)) => {
            use super::spec::AlignFilterCriterion;
            use crate::core::align::filter::{AlignmentFiltering, FilteringParameters};
            use crate::helper::finder::IDs;
            use crate::helper::types::{DataType, InputFmt, OutputFmt, PartitionFmt};

            let input_fmt: InputFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            let params = match &p.criterion {
                AlignFilterCriterion::MinTaxPercent { percent, .. } => {
                    let count = IDs::new(primary, &input_fmt, &datatype).id_unique().len();
                    FilteringParameters::MinTax((*percent * count as f64).floor() as usize)
                }
                AlignFilterCriterion::MinTaxCount { min_taxa } => {
                    FilteringParameters::MinTax(*min_taxa)
                }
                AlignFilterCriterion::MinLen { min_len } => FilteringParameters::MinLen(*min_len),
                AlignFilterCriterion::MaxLen { max_len } => FilteringParameters::MaxLen(*max_len),
                AlignFilterCriterion::MinParsInf { min_pars_inf } => {
                    FilteringParameters::MinParsInf(*min_pars_inf)
                }
                AlignFilterCriterion::MaxParsInf { max_pars_inf } => {
                    FilteringParameters::MaxParsInf(*max_pars_inf)
                }
                AlignFilterCriterion::PercentInf { percent_inf, .. } => {
                    FilteringParameters::PercInf(*percent_inf)
                }
                AlignFilterCriterion::MissingData { max_missing } => {
                    FilteringParameters::MissingData(*max_missing)
                }
                AlignFilterCriterion::TaxonAll { taxa } => {
                    FilteringParameters::TaxonAll(taxa.clone())
                }
            };
            let mut filter =
                AlignmentFiltering::new(primary, &input_fmt, &datatype, output_dir, &params)
                    .with_report_context(ctx.clone());
            let out_fmt_parsed: Option<OutputFmt> = p
                .output_format
                .as_deref()
                .map(|s| s.parse().map_err(anyhow::Error::msg))
                .transpose()?;
            let part_fmt_parsed: Option<PartitionFmt> = p
                .partition_format
                .as_deref()
                .map(|s| s.parse().map_err(anyhow::Error::msg))
                .transpose()?;
            let prefix_path = p.prefix.as_deref().map(PathBuf::from);
            if p.concat
                && let (Some(of), Some(pf), Some(px)) =
                    (&out_fmt_parsed, &part_fmt_parsed, &prefix_path)
            {
                filter.set_concat(of, pf, px);
            }
            filter.filter_with_report(p.clone(), auxiliary)
        }
        #[cfg(feature = "align")]
        (TaskKind::AlignSplit, TaskParameters::AlignSplit(p)) => {
            use crate::core::align::split::AlignmentSplitting;
            use crate::helper::types::{DataType, InputFmt, OutputFmt, PartitionFmt};
            let input_path = primary
                .first()
                .context("Missing primary input for align.split")?;
            let input_fmt: InputFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            let output_fmt: OutputFmt = p.output_format.parse().map_err(anyhow::Error::msg)?;
            let part_fmt: PartitionFmt = p.partition_format.parse().map_err(anyhow::Error::msg)?;
            let part_path = aux_path(auxiliary, "partition")
                .map(Path::to_path_buf)
                .unwrap_or_else(|| input_path.clone());
            AlignmentSplitting::new(input_path, &datatype, &input_fmt, output_dir, &output_fmt)
                .with_report_context(ctx.clone())
                .try_split_with_report(&part_path, &part_fmt, &p.prefix, p.skip_checking)
        }
        #[cfg(feature = "align")]
        (TaskKind::AlignTrim, TaskParameters::AlignTrim(p)) => {
            use super::spec::AlignTrimCriterion;
            use crate::core::align::trim::{AlignmentTrimming, TrimmingParameters};
            use crate::helper::types::{DataType, InputFmt, OutputFmt};
            let input_fmt: InputFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            let output_fmt: OutputFmt = p.output_format.parse().map_err(anyhow::Error::msg)?;
            let trim_params = match p.criterion {
                AlignTrimCriterion::MissingData { max_missing } => {
                    TrimmingParameters::MissingData(max_missing)
                }
                AlignTrimCriterion::ParsInf { min_pars_inf } => {
                    TrimmingParameters::ParsInf(min_pars_inf)
                }
            };
            AlignmentTrimming::new(
                primary,
                &input_fmt,
                &datatype,
                output_dir,
                &output_fmt,
                &trim_params,
            )
            .with_report_context(ctx.clone())
            .trim_with_report()
        }
        #[cfg(feature = "align")]
        (TaskKind::AlignUnalign, TaskParameters::AlignUnalign(p)) => {
            use crate::core::align::unalign::UnalignAlignment;
            use crate::helper::types::{DataType, InputFmt, OutputFmt};
            let input_fmt: InputFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            let output_fmt: OutputFmt = p.output_format.parse().map_err(anyhow::Error::msg)?;
            UnalignAlignment::new(primary, &input_fmt, &datatype, output_dir, &output_fmt)
                .with_report_context(ctx.clone())
                .unalign_with_report()
        }
        #[cfg(feature = "align")]
        (TaskKind::PartitionConvert, TaskParameters::PartitionConvert(p)) => {
            use crate::core::align::partition::PartitionConverter;
            use crate::helper::partition::construct_partition_path;
            use crate::helper::types::{DataType, PartitionFmt};
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            let in_fmt: PartitionFmt = p
                .input_partition_format
                .parse()
                .map_err(anyhow::Error::msg)?;
            let out_fmt: PartitionFmt = p
                .output_partition_format
                .parse()
                .map_err(anyhow::Error::msg)?;
            fs::create_dir_all(output_dir)?;
            let mut outputs = Vec::new();
            for input in primary {
                let stem = input
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .context("Invalid partition filename")?;
                let out_path = construct_partition_path(&output_dir.join(stem), &out_fmt);
                PartitionConverter::new(input, &in_fmt, &out_path, &out_fmt)
                    .try_convert(&datatype, p.skip_checking)?;
                outputs.push(("partition", out_path));
            }
            let counts = super::stats::TaskCounts {
                files_read: primary.len(),
                files_written: outputs.len(),
                items_processed: Some(outputs.len()),
                items_skipped: None,
            };
            let (_, report) = super::write_task_report(
                ctx,
                TaskKind::PartitionConvert,
                primary,
                Vec::new(),
                p,
                output_dir,
                outputs,
                Some(&counts),
                None,
            )?;
            Ok(report)
        }
        #[cfg(feature = "contig")]
        (TaskKind::ContigSummary, TaskParameters::ContigSummary(p)) => {
            use crate::core::contig::summarize::ContigSummaryHandler;
            use crate::helper::types::ContigFmt;
            let input_fmt: ContigFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let files =
                ContigSummaryHandler::new(primary, &input_fmt, output_dir, p.prefix.as_deref())
                    .with_report_context(ctx.clone())
                    .try_summarize()?;
            RunReport::from_path(&files.report)
        }
        #[cfg(feature = "read")]
        (TaskKind::ReadSummary, TaskParameters::ReadSummary(p)) => {
            use crate::core::read::summarize::GenomicReadSummary;
            use crate::helper::types::{PhredOffset, SeqReadFmt, SummaryMode};
            use crate::stats::fastq::qc::QcProfile;
            let input_fmt: SeqReadFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let mode: SummaryMode = p.mode.parse().map_err(anyhow::Error::msg)?;
            let phred: PhredOffset = p.phred.parse().map_err(anyhow::Error::msg)?;
            let qc_profile: QcProfile = p.qc_profile.parse().map_err(anyhow::Error::msg)?;
            let mut inputs = primary.to_vec();
            let files = GenomicReadSummary::new(
                &mut inputs,
                &input_fmt,
                &mode,
                output_dir,
                p.prefix.as_deref(),
            )
            .with_phred(phred)
            .with_dup_sketch_size(p.dup_sketch_size)
            .with_qc_profile(qc_profile)
            .with_report_context(ctx.clone())
            .try_summarize()?;
            let report_path = files.report.context("Missing read report path")?;
            RunReport::from_path(&report_path)
        }
        #[cfg(feature = "maf")]
        (TaskKind::MafSummary, TaskParameters::MafSummary(p)) => {
            use crate::core::maf::SampleNaming;
            use crate::core::maf::summary::MafSummary;
            let naming = match p.naming.as_str() {
                "full_source" => SampleNaming::FullSource,
                _ => SampleNaming::AssemblyPrefix,
            };
            let files = MafSummary::new(primary, output_dir, p.prefix.as_deref())
                .with_naming(naming)
                .with_ref_species(p.ref_species.clone())
                .with_species(p.species.clone())
                .with_report_context(ctx.clone())
                .try_summarize_with_report()?;
            Ok(files.run_report)
        }
        #[cfg(feature = "maf")]
        (TaskKind::MafConvert, TaskParameters::MafConvert(p)) => {
            use crate::core::maf::BlockFilter;
            use crate::core::maf::SampleNaming;
            use crate::core::maf::convert::{MafConverter, MafGrouping};
            use crate::helper::types::OutputFmt;
            let output_fmt: OutputFmt = p.output_format.parse().map_err(anyhow::Error::msg)?;
            let grouping: MafGrouping = p.group_by.parse().map_err(anyhow::Error::msg)?;
            let naming = match p.naming.as_str() {
                "full_source" => SampleNaming::FullSource,
                _ => SampleNaming::AssemblyPrefix,
            };
            let filter = BlockFilter {
                min_ref_len: p.min_block_len,
                min_species: p.min_species,
                min_score: p.min_score,
                species: p.species.clone(),
            };
            let mut converter = MafConverter::new(primary, output_dir, &output_fmt)
                .with_grouping(grouping)
                .with_clip(p.clip)
                .with_strand_aware(p.strand_aware)
                .with_naming(naming)
                .with_ref_species(p.ref_species.clone())
                .with_filter(filter)
                .with_unmatched(p.write_unmatched)
                .with_prefix(p.prefix.clone())
                .with_report_context(ctx.clone());
            if let Some(bed) = aux_path(auxiliary, "bed") {
                converter = converter.with_bed(bed);
            }
            converter.convert_with_report(auxiliary)
        }
        #[cfg(feature = "sequence")]
        (TaskKind::SequenceAdd, TaskParameters::SequenceAdd(p)) => {
            use crate::core::sequence::addition::SequenceAddition;
            use crate::helper::types::{DataType, InputFmt, OutputFmt};
            let input_fmt: InputFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let dest_fmt: InputFmt = p.dest_format.parse().map_err(anyhow::Error::msg)?;
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            let output_fmt: OutputFmt = p.output_format.parse().map_err(anyhow::Error::msg)?;
            let dests: Vec<PathBuf> = auxiliary
                .iter()
                .filter(|(r, _)| r == "destination")
                .map(|(_, path)| path.clone())
                .collect();
            SequenceAddition::new(
                primary,
                &input_fmt,
                &datatype,
                output_dir,
                &output_fmt,
                p.added_only,
            )
            .with_report_context(ctx.clone())
            .add_with_report(&dests, &dest_fmt, p.single_file, p.include_filename)
        }
        #[cfg(feature = "sequence")]
        (TaskKind::SequenceExtract, TaskParameters::SequenceExtract(p)) => {
            use crate::core::sequence::extract::{SeqExtractionParameters, SequenceExtraction};
            use crate::helper::types::{DataType, InputFmt, OutputFmt};
            let input_fmt: InputFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            let output_fmt: OutputFmt = p.output_format.parse().map_err(anyhow::Error::msg)?;
            let opts = if let Some(re) = &p.regex {
                SeqExtractionParameters::Regex(re.clone())
            } else if let Some(ids) = &p.ids {
                SeqExtractionParameters::Id(ids.clone())
            } else {
                SeqExtractionParameters::None
            };
            SequenceExtraction::new(&input_fmt, &datatype, &opts, output_dir, &output_fmt)
                .with_report_context(ctx.clone())
                .extract_with_report(primary, p.clone(), auxiliary)
        }
        #[cfg(feature = "sequence")]
        (TaskKind::SequenceFilter, TaskParameters::SequenceFilter(p)) => {
            use super::spec::SequenceFilterCriterion;
            use crate::core::sequence::filter::{SeqFilteringParameters, SequenceFiltering};
            use crate::helper::types::{DataType, InputFmt, OutputFmt};
            let input_fmt: InputFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            let output_fmt: OutputFmt = p.output_format.parse().map_err(anyhow::Error::msg)?;
            let params = match p.criterion {
                SequenceFilterCriterion::MinLength { min_length } => {
                    SeqFilteringParameters::MinSequenceLength(min_length)
                }
                SequenceFilterCriterion::MaxLength { max_length } => {
                    SeqFilteringParameters::MaxSequenceLength(max_length)
                }
                SequenceFilterCriterion::MaxGap { max_gap } => {
                    SeqFilteringParameters::PercentMaxGap(max_gap)
                }
            };
            SequenceFiltering::new(
                primary,
                &input_fmt,
                &datatype,
                output_dir,
                &output_fmt,
                &params,
            )
            .with_report_context(ctx.clone())
            .filter_with_report()
        }
        #[cfg(feature = "sequence")]
        (TaskKind::SequenceId, TaskParameters::SequenceId(p)) => {
            use crate::core::sequence::id::SequenceID;
            use crate::helper::types::{DataType, InputFmt};
            let input_fmt: InputFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            SequenceID::new(
                primary,
                &input_fmt,
                &datatype,
                output_dir,
                p.prefix.as_deref(),
            )
            .with_report_context(ctx.clone())
            .run_with_report(p.map_id)
        }
        #[cfg(feature = "sequence")]
        (TaskKind::SequenceRemove, TaskParameters::SequenceRemove(p)) => {
            use crate::core::sequence::remove::{SeqRemovalParameters, SequenceRemoval};
            use crate::helper::types::{DataType, InputFmt, OutputFmt};
            let input_fmt: InputFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            let output_fmt: OutputFmt = p.output_format.parse().map_err(anyhow::Error::msg)?;
            let opts = if let Some(re) = &p.regex {
                SeqRemovalParameters::Regex(re.clone())
            } else {
                SeqRemovalParameters::Id(p.ids.clone().unwrap_or_default())
            };
            SequenceRemoval::new(&input_fmt, &datatype, output_dir, &output_fmt, &opts)
                .with_report_context(ctx.clone())
                .remove_with_report(primary)
        }
        #[cfg(feature = "sequence")]
        (TaskKind::SequenceRename, TaskParameters::SequenceRename(p)) => {
            use crate::core::sequence::rename::{SeqRenamingParameters, SequenceRenaming};
            use crate::helper::types::{DataType, InputFmt, OutputFmt};
            use crate::parser::delimited;
            let input_fmt: InputFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            let output_fmt: OutputFmt = p.output_format.parse().map_err(anyhow::Error::msg)?;
            let opts = match p.mode.as_str() {
                "id_map" => {
                    let map_path = aux_path(auxiliary, "rename_map")
                        .context("Missing rename_map auxiliary input for sequence.rename")?;
                    SeqRenamingParameters::RnId(delimited::parse_delimited_text(map_path))
                }
                "remove_str" => SeqRenamingParameters::RmStr(p.from.clone().unwrap_or_default()),
                "remove_regex" => SeqRenamingParameters::RmRegex(
                    p.from.clone().unwrap_or_default(),
                    p.all_matches,
                ),
                "replace_str" => SeqRenamingParameters::RpStr(
                    p.from.clone().unwrap_or_default(),
                    p.to.clone().unwrap_or_default(),
                ),
                "replace_regex" => SeqRenamingParameters::RpRegex(
                    p.from.clone().unwrap_or_default(),
                    p.to.clone().unwrap_or_default(),
                    p.all_matches,
                ),
                _ => SeqRenamingParameters::None,
            };
            SequenceRenaming::new(&input_fmt, &datatype, output_dir, &output_fmt, &opts)
                .with_report_context(ctx.clone())
                .rename_with_report(primary, p.clone(), auxiliary)
        }
        #[cfg(feature = "sequence")]
        (TaskKind::SequenceTranslate, TaskParameters::SequenceTranslate(p)) => {
            use crate::core::sequence::translate::SequenceTranslation;
            use crate::helper::types::{DataType, GeneticCodes, InputFmt, OutputFmt};
            let input_fmt: InputFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            let output_fmt: OutputFmt = p.output_format.parse().map_err(anyhow::Error::msg)?;
            let table: GeneticCodes = p.table.parse().map_err(anyhow::Error::msg)?;
            let trans = SequenceTranslation::new(&input_fmt, &table, &datatype, &output_fmt)
                .with_report_context(ctx.clone());
            if p.auto_frame {
                trans.translate_all_autoframe_with_report(primary, output_dir)
            } else {
                trans.translate_with_report(primary, p.reading_frame, output_dir)
            }
        }
        #[cfg(feature = "sequence")]
        (TaskKind::SequenceTranspose, TaskParameters::SequenceTranspose(p)) => {
            use crate::core::sequence::transpose::SequenceTranspose;
            use crate::helper::types::{DataType, InputFmt, OutputFmt};
            let input_fmt: InputFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            let output_fmt: OutputFmt = p.output_format.parse().map_err(anyhow::Error::msg)?;
            SequenceTranspose::new(&input_fmt, &datatype, output_dir, &output_fmt)
                .with_report_context(ctx.clone())
                .transpose_with_report(primary)
        }
        #[cfg(feature = "tree")]
        (TaskKind::TreeConvert, TaskParameters::TreeConvert(p)) => {
            use crate::core::tree::convert::TreeConversion;
            use crate::parser::tree::{ParseOptions, TreeFmt, TreeOutputFmt};
            use crate::writer::tree::TreeWriteOptions;
            let input_fmt: TreeFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let output_fmt: TreeOutputFmt = p.output_format.parse().map_err(anyhow::Error::msg)?;
            let parse_opts = ParseOptions {
                parse_annotations: !p.strip_annotations,
                ..Default::default()
            };
            let write_opts = TreeWriteOptions {
                use_translate: p.translate,
                write_annotations: !p.strip_annotations,
                gzip: p.gzip,
                ..Default::default()
            };
            TreeConversion::new(&input_fmt, &output_fmt)
                .with_parse_options(parse_opts)
                .with_write_options(write_opts)
                .with_report_context(ctx.clone())
                .convert_with_report(primary, output_dir)
        }
        #[cfg(feature = "tree")]
        (TaskKind::TreeInfer, TaskParameters::TreeInfer(p)) => {
            use crate::core::tree::estimate::{
                BranchSwap, GapMode, ParsimonyBranchLengths, ParsimonyOptions, ParsimonySearch,
                RootingMethod, TreeEstimateOptions, TreeEstimation, TreeMethod,
            };
            use crate::helper::types::{DataType, InputFmt};
            use crate::parser::tree::TreeOutputFmt;
            use crate::stats::distance::DistanceModel;
            use crate::writer::tree::TreeWriteOptions;

            let input_fmt: InputFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let datatype: DataType = p.datatype.parse().map_err(anyhow::Error::msg)?;
            let method: TreeMethod = p.method.parse().map_err(anyhow::Error::msg)?;
            let model: DistanceModel = p
                .model
                .as_deref()
                .unwrap_or("jc69")
                .parse()
                .map_err(anyhow::Error::msg)?;
            let search: ParsimonySearch = p
                .parsimony_search
                .as_deref()
                .unwrap_or("auto")
                .parse()
                .map_err(anyhow::Error::msg)?;
            let swap: BranchSwap = p
                .swap
                .as_deref()
                .unwrap_or("tbr")
                .parse()
                .map_err(anyhow::Error::msg)?;
            let gap_mode: GapMode = p
                .gap_mode
                .as_deref()
                .unwrap_or("missing")
                .parse()
                .map_err(anyhow::Error::msg)?;
            let branch_lengths: ParsimonyBranchLengths = p
                .mp_branch_lengths
                .as_deref()
                .unwrap_or("acctran")
                .parse()
                .map_err(anyhow::Error::msg)?;
            let rooting = match p.rooting.as_str() {
                "midpoint" => RootingMethod::Midpoint,
                "outgroup" => RootingMethod::Outgroup(p.outgroup.clone().unwrap_or_default()),
                _ => RootingMethod::None,
            };
            let output_fmt: TreeOutputFmt = p.output_format.parse().map_err(anyhow::Error::msg)?;
            let write_opts = TreeWriteOptions {
                use_translate: p.translate,
                gzip: p.gzip,
                ..Default::default()
            };
            let options = TreeEstimateOptions {
                method,
                model,
                parsimony: ParsimonyOptions {
                    search,
                    swap,
                    replicates: p.replicates.unwrap_or(10),
                    gap_mode,
                    branch_lengths,
                    seed: p.seed,
                },
                bootstrap: p.bootstrap,
                seed: p.seed,
                rooting,
                output_fmt,
                single_file: p.single_file,
                write_opts,
            };
            let (_, report) = TreeEstimation::new(&input_fmt, &datatype, options)
                .with_report_context(ctx.clone())
                .estimate_with_report(primary, output_dir, p.prefix.as_deref())?;
            Ok(report)
        }
        #[cfg(feature = "tree")]
        (TaskKind::TreeSummary, TaskParameters::TreeSummary(p)) => {
            use crate::core::tree::estimate::RootingMethod;
            use crate::core::tree::summary::{
                BranchLengthStrategy, BurninSpec, MccScoring, RootedMode, SummaryMethod,
                SupportFormat, SupportPlacement, TreeSummarization, TreeSummaryOptions,
            };
            use crate::parser::tree::{TreeFmt, TreeOutputFmt};
            use crate::writer::tree::TreeWriteOptions;

            let input_fmt: TreeFmt = p.input_format.parse().map_err(anyhow::Error::msg)?;
            let method: SummaryMethod = p.method.parse().map_err(anyhow::Error::msg)?;
            let burnin: BurninSpec = p.burnin.parse().map_err(anyhow::Error::msg)?;
            let rooted_mode: RootedMode = p.rooted_mode.parse().map_err(anyhow::Error::msg)?;
            let mcc_scoring: MccScoring = p.mcc_scoring.parse().map_err(anyhow::Error::msg)?;
            let branch_lengths: BranchLengthStrategy =
                p.branch_lengths.parse().map_err(anyhow::Error::msg)?;
            let support_format: SupportFormat =
                p.support_format.parse().map_err(anyhow::Error::msg)?;
            let support_placement: SupportPlacement =
                p.support_placement.parse().map_err(anyhow::Error::msg)?;
            let output_fmt: TreeOutputFmt = p.output_format.parse().map_err(anyhow::Error::msg)?;
            let rooting = match p.rooting.as_str() {
                "midpoint" => RootingMethod::Midpoint,
                "outgroup" => RootingMethod::Outgroup(p.outgroup.clone().unwrap_or_default()),
                _ => RootingMethod::None,
            };
            let target_tree = aux_path(auxiliary, "target").map(Path::to_path_buf);
            let options = TreeSummaryOptions {
                method,
                target_tree,
                burnin,
                min_clade_prob: Some(p.min_clade_prob),
                rooted_mode,
                mcc_scoring,
                branch_lengths,
                support_format,
                support_placement,
                hpd_level: p.hpd_level,
                summarize_annotations: p.summarize_annotations,
                allow_partial_taxa: p.allow_partial_taxa,
                rooting,
                output_fmt,
                write_opts: TreeWriteOptions {
                    use_translate: p.translate,
                    write_annotations: !p.strip_annotations,
                    gzip: p.gzip,
                    ..Default::default()
                },
            };
            let (_, report) = TreeSummarization::new(&input_fmt, options)
                .with_report_context(ctx.clone())
                .summarize_with_report(primary, output_dir, p.prefix.as_deref())?;
            Ok(report)
        }
        _ => bail!("Task '{}' is not enabled in this build of segul-core", kind),
    }
}

#[cfg_attr(not(any(feature = "align", feature = "sequence")), allow(dead_code))]
fn aux_path<'a>(auxiliary: &'a [(String, PathBuf)], role: &str) -> Option<&'a Path> {
    auxiliary
        .iter()
        .find(|(r, _)| r == role)
        .map(|(_, p)| p.as_path())
}
