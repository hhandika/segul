//! Multi-Sequence Alignment Module.
//!
//! Parse user input and call the appropriate functions
//! to concatenate sequence alignments.

use std::path::{Path, PathBuf};

use colored::Colorize;

use crate::helper::concat::Concat;
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, PartitionFmt};
use crate::helper::{files, utils};
use crate::report::spec::{AlignConcatParams, TaskKind, canonical};
use crate::report::stats::ConcatStats;
use crate::report::{ReportContext, RunReport, RunReportBuilder, report_path};
use crate::writer::{partition::PartWriter, sequences::SeqWriter};

pub struct AlignmentConcatenation<'a> {
    input_fmt: &'a InputFmt,
    output_dir: &'a Path,
    output_fmt: &'a OutputFmt,
    part_fmt: &'a PartitionFmt,
    output_prefix: &'a Path,
    sort: bool,
    report_context: ReportContext,
}

impl<'a> AlignmentConcatenation<'a> {
    pub fn new(
        input_fmt: &'a InputFmt,
        output_dir: &'a Path,
        output_fmt: &'a OutputFmt,
        part_fmt: &'a PartitionFmt,
        output_prefix: &'a Path,
    ) -> Self {
        Self {
            input_fmt,
            output_dir,
            output_fmt,
            part_fmt,
            output_prefix,
            sort: false,
            report_context: ReportContext::default(),
        }
    }

    /// Sort sequences alphabetically by ID before writing the concatenated matrix.
    pub fn with_sort(mut self, sort: bool) -> Self {
        self.sort = sort;
        self
    }

    /// Attach frontend provenance and SHA-256 checksum settings for `report.json`.
    pub fn with_report_context(mut self, context: ReportContext) -> Self {
        self.report_context = context;
        self
    }

    pub fn concat(&mut self, files: &mut [PathBuf], datatype: &DataType) {
        self.concat_with_report(files, datatype, None)
            .expect("Failed concatenating alignments");
    }

    /// Concatenate alignments and write a unified `report.json`.
    pub fn concat_with_report(
        &mut self,
        files: &mut [PathBuf],
        datatype: &DataType,
        prefix_name: Option<&str>,
    ) -> anyhow::Result<RunReport> {
        alphanumeric_sort::sort_path_slice(files);
        let prefix_label = prefix_name.map(str::to_string).or_else(|| {
            if self.output_prefix != self.output_dir {
                self.output_prefix
                    .file_name()
                    .and_then(|n| n.to_str())
                    .map(str::to_string)
            } else {
                None
            }
        });
        let params = AlignConcatParams {
            input_format: canonical::input_fmt(self.input_fmt).to_string(),
            datatype: canonical::datatype(datatype).to_string(),
            output_format: canonical::output_fmt(self.output_fmt).to_string(),
            partition_format: canonical::partition_fmt(self.part_fmt).to_string(),
            sort: self.sort,
            prefix: prefix_label.clone(),
        };
        let builder = RunReportBuilder::new(TaskKind::AlignConcat, &self.report_context)
            .primary_inputs(files)
            .parameters(&params)?;
        let mut concat = Concat::new(files, self.input_fmt, datatype);
        let output_path =
            files::create_output_fname(self.output_dir, self.output_prefix, self.output_fmt);
        let spin = utils::set_spinner();
        concat.concat_alignment(&spin);
        if self.sort {
            concat.alignment.sort_keys();
        }
        let mut seq_writer = SeqWriter::new(&output_path, &concat.alignment, &concat.header);
        let part_fname = self.construct_part_fpath(&output_path);
        let part_writer = PartWriter::new(&part_fname, &concat.partition, self.part_fmt, datatype);
        spin.set_message("Writing output files...");
        seq_writer
            .write_sequence(self.output_fmt)
            .expect("Failed writing the output file");
        part_writer.write_partition();
        spin.finish_with_message("Finished concatenating alignments!\n");
        self.print_output_info(concat.partition.len(), &concat.header);
        let mut out_files = vec![("alignment", output_path)];
        if !matches!(
            self.part_fmt,
            PartitionFmt::Charset | PartitionFmt::CharsetCodon
        ) {
            out_files.push(("partition", part_fname));
        }
        let stats = ConcatStats {
            taxa: concat.header.ntax,
            loci: concat.partition.len(),
            sites: concat.header.nchar,
            partitions: concat.partition.len(),
        };
        let report_file = report_path(self.output_dir, prefix_label.as_deref());
        builder
            .outputs(self.output_dir, out_files)
            .statistics(&stats)?
            .write_to(&report_file)
    }

    fn construct_part_fpath(&mut self, output_path: &Path) -> PathBuf {
        match self.part_fmt {
            PartitionFmt::Charset | PartitionFmt::CharsetCodon => PathBuf::from(output_path),
            PartitionFmt::Nexus | PartitionFmt::NexusCodon => {
                self.get_part_fname(output_path, "nex")
            }
            PartitionFmt::Raxml | PartitionFmt::RaxmlCodon => {
                self.get_part_fname(output_path, "txt")
            }
        }
    }

    fn get_part_fname(&self, output_path: &Path, ext: &str) -> PathBuf {
        let fname = format!(
            "{}_partition.{}",
            output_path
                .file_stem()
                .expect("Failed getting file name for partition")
                .to_string_lossy(),
            ext
        );

        output_path
            .parent()
            .expect("Failed getting output parent directory")
            .join(Path::new(&fname))
    }

    fn print_output_info(&self, count: usize, header: &Header) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Taxa", utils::fmt_num(&header.ntax));
        log::info!("{:18}: {}", "Loci", utils::fmt_num(&count));
        log::info!(
            "{:18}: {}",
            "Alignment length",
            utils::fmt_num(&header.nchar)
        );
        log::info!("{:18}: {}", "Directory", self.output_dir.display());
    }
}
