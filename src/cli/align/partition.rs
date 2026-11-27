use std::path::{Path, PathBuf};

use segul_core::report::spec::{PartitionConvertParams, canonical};
use segul_core::report::stats::TaskCounts;
use segul_core::report::{ReportContext, TaskKind, write_task_report};

use crate::core::align::partition::PartitionConverter;
use crate::helper::partition::{construct_partition_path, infer_partition_fmt};
use crate::helper::{logger, utils};

use crate::cli::args::align::PartitionArgs;
use crate::cli::{AlignSeqInput, ConcatCli, InputCli, OutputCli};

impl InputCli for PartParser<'_> {}
impl ConcatCli for PartParser<'_> {}
impl OutputCli for PartParser<'_> {}
impl AlignSeqInput for PartParser<'_> {}

pub(in crate::cli) struct PartParser<'a> {
    args: &'a PartitionArgs,
}

impl<'a> PartParser<'a> {
    pub(in crate::cli) fn new(args: &'a PartitionArgs) -> Self {
        Self { args }
    }

    pub(in crate::cli) fn convert(&self, ctx: &ReportContext) {
        let inputs = self.collect_paths(&self.args.input);
        let input_counts = inputs.len();
        let explicit_fmt = self
            .args
            .part_fmt
            .as_deref()
            .map(|fmt| self.parse_partition_fmt_std(fmt));
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let out_part_fmt = self.parse_partition_fmt(&self.args.out_part, self.args.codon);

        let mut converted_outputs: Vec<(&str, PathBuf)> = Vec::with_capacity(input_counts);
        let mut first_in_fmt = None;
        inputs.iter().for_each(|input| {
            logger::log_input_partition(Some(input), input_counts);
            let in_part_fmt = explicit_fmt.unwrap_or_else(|| {
                infer_partition_fmt(input).unwrap_or_else(|err| panic!("{err}"))
            });
            if first_in_fmt.is_none() {
                first_in_fmt = Some(in_part_fmt);
            }
            let output_base = match &self.args.output {
                Some(out_dir) => {
                    let fname = input
                        .file_name()
                        .expect("Failed to parse input partition file name");
                    out_dir.join(fname)
                }
                None => input.to_path_buf(),
            };
            let output = construct_partition_path(&output_base, &out_part_fmt);
            self.check_output_file_exist(&output, self.args.force);
            let converter = PartitionConverter::new(input, &in_part_fmt, &output, &out_part_fmt);
            converter
                .try_convert(&datatype, self.args.skip_checking)
                .unwrap_or_else(|err| panic!("{err}"));
            converted_outputs.push(("partition", output));
            if input_counts > 1 {
                utils::print_divider();
            }
        });
        if let Some(in_part_fmt) = first_in_fmt {
            let report_dir = self.args.output.clone().unwrap_or_else(|| {
                converted_outputs
                    .first()
                    .and_then(|(_, p)| p.parent())
                    .filter(|p| !p.as_os_str().is_empty())
                    .unwrap_or_else(|| Path::new("."))
                    .to_path_buf()
            });
            let params = PartitionConvertParams {
                datatype: canonical::datatype(&datatype).to_string(),
                input_partition_format: canonical::partition_fmt(&in_part_fmt).to_string(),
                output_partition_format: canonical::partition_fmt(&out_part_fmt).to_string(),
                skip_checking: self.args.skip_checking,
            };
            let counts = TaskCounts {
                files_read: inputs.len(),
                files_written: converted_outputs.len(),
                items_processed: Some(converted_outputs.len()),
                items_skipped: None,
            };
            write_task_report(
                ctx,
                TaskKind::PartitionConvert,
                &inputs,
                Vec::new(),
                &params,
                &report_dir,
                converted_outputs,
                Some(&counts),
                None,
            )
            .unwrap_or_else(|err| panic!("{err}"));
        }
    }
}
