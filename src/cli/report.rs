//! CLI handlers for `segul run` and `segul report`.

use std::fs;
use std::path::Path;

use colored::Colorize;
use segul_core::report::methods::render_methods_from_path;
use segul_core::report::replay::{ReplayOptions, replay_report, replay_workflow};
use segul_core::report::verify::{FileVerifyStatus, verify_report};
use segul_core::report::workflow::WorkflowRecipe;
use segul_core::report::{FrontendInfo, ReportContext, RunReport};

use crate::cli::OutputCli;
use crate::cli::args::{ReportMethodsArgs, ReportVerifyArgs, ReportWorkflowArgs, RunArgs};

pub(in crate::cli) struct RunCliParser<'a> {
    args: &'a RunArgs,
}

impl OutputCli for RunCliParser<'_> {}

impl<'a> RunCliParser<'a> {
    pub(in crate::cli) fn new(args: &'a RunArgs) -> Self {
        Self { args }
    }

    pub(in crate::cli) fn run(&self, ctx: &ReportContext) {
        if let Some(output) = &self.args.output {
            self.check_output_dir_exist(output, self.args.force);
        }
        let opts = ReplayOptions {
            input_base: self.args.input_dir.clone(),
            output_dir: self.args.output.clone(),
            verify_inputs: self.args.verify_inputs,
            checksum_mode: Some(ctx.checksum_mode),
            frontend: FrontendInfo::cli(clap::crate_version!()),
        };
        if RunReport::is_run_report_file(&self.args.file) {
            match replay_report(&self.args.file, &opts) {
                Ok(report) => {
                    log::info!("{}", "Replay Complete".yellow());
                    log::info!("{:18}: {}", "Task", report.kind);
                    log::info!("{:18}: {}", "Output dir", report.outputs.output_dir);
                }
                Err(err) => {
                    log::error!("{:#}", err);
                    std::process::exit(1);
                }
            }
        } else {
            match replay_workflow(&self.args.file, &opts) {
                Ok(reports) => {
                    log::info!("{}", "Workflow Complete".yellow());
                    log::info!("{:18}: {}", "Steps executed", reports.len());
                    if let Some(last) = reports.last() {
                        log::info!("{:18}: {}", "Final output", last.outputs.output_dir);
                    }
                }
                Err(err) => {
                    log::error!("{:#}", err);
                    std::process::exit(1);
                }
            }
        }
    }
}

pub(in crate::cli) struct ReportCliParser;

impl ReportCliParser {
    pub(in crate::cli) fn verify(args: &ReportVerifyArgs) {
        let result = match verify_report(
            &args.report,
            args.input_dir.as_deref(),
            args.output_dir.as_deref(),
            args.scope,
        ) {
            Ok(res) => res,
            Err(err) => {
                eprintln!("{err:#}");
                std::process::exit(1);
            }
        };

        if args.json {
            match result.to_json_pretty() {
                Ok(text) => println!("{text}"),
                Err(err) => {
                    eprintln!("Failed serializing verification result: {err}");
                    std::process::exit(1);
                }
            }
        } else {
            let status_label = if result.passed {
                "PASSED".green()
            } else {
                "FAILED".red()
            };
            println!(
                "Verification {} ({}: checked {}, ok {}, mismatched {}, missing {}, unchecked {})",
                status_label,
                result.scope,
                result.checked,
                result.ok,
                result.mismatched,
                result.missing,
                result.unchecked
            );
            for (group, entries) in [
                ("input", &result.inputs),
                ("auxiliary", &result.auxiliary),
                ("output", &result.outputs),
                ("parent", &result.parents),
            ] {
                for entry in entries {
                    match &entry.status {
                        FileVerifyStatus::Ok | FileVerifyStatus::Unchecked => {}
                        FileVerifyStatus::Missing => {
                            eprintln!(
                                "  [MISSING] {group}: {} ({})",
                                entry.path, entry.resolved_path
                            );
                        }
                        FileVerifyStatus::Mismatch { expected, actual } => {
                            eprintln!(
                                "  [MISMATCH] {group}: {} (expected {expected}, got {actual})",
                                entry.path
                            );
                        }
                        FileVerifyStatus::SizeMismatch { expected, actual } => {
                            eprintln!(
                                "  [SIZE MISMATCH] {group}: {} (expected {expected} bytes, got {actual})",
                                entry.path
                            );
                        }
                    }
                }
            }
        }

        if !result.passed {
            std::process::exit(1);
        }
    }

    pub(in crate::cli) fn methods(args: &ReportMethodsArgs) {
        let text = match render_methods_from_path(&args.report, !args.no_chain, args.format) {
            Ok(rendered) => rendered,
            Err(err) => {
                eprintln!("{err:#}");
                std::process::exit(1);
            }
        };
        if let Some(out_path) = &args.output {
            if let Err(err) = write_text_file(out_path, &text) {
                eprintln!("{err:#}");
                std::process::exit(1);
            }
        } else {
            print!("{text}");
        }
    }

    pub(in crate::cli) fn workflow(args: &ReportWorkflowArgs) {
        let recipe = if args.reports.len() == 1 {
            WorkflowRecipe::from_report_chain(&args.reports[0], args.name.as_deref())
        } else {
            WorkflowRecipe::from_report_paths(&args.reports, args.name.as_deref())
        };
        let recipe = match recipe {
            Ok(r) => r,
            Err(err) => {
                log::error!("{:#}", err);
                std::process::exit(1);
            }
        };
        if let Err(err) = recipe.write(&args.output) {
            log::error!("{:#}", err);
            std::process::exit(1);
        }
        log::info!("{}", "Workflow Recipe".yellow());
        log::info!("{:18}: {}", "Steps", recipe.steps.len());
        log::info!("{:18}: {}", "Output path", args.output.display());
    }
}

fn write_text_file(path: &Path, text: &str) -> anyhow::Result<()> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    fs::write(path, text)?;
    Ok(())
}
