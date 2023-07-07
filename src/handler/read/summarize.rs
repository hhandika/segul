//! A handler for summarizing raw sequence data.

use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
    sync::mpsc::channel,
};

use colored::Colorize;

use indicatif::ProgressBar;
use rayon::prelude::*;

use crate::{
    helper::{
        files,
        types::{infer_raw_input_auto, SeqReadFmt, SummaryMode},
        utils::set_spinner,
    },
    stats::fastq::{self, FastqSummary},
    writer::read::ReadSummaryWriter,
};

/// Include support for any compressed or uncompressed fastq files.
pub struct ReadSummaryHandler<'a> {
    pub inputs: &'a mut [PathBuf],
    pub input_fmt: &'a SeqReadFmt,
    /// Summary mode
    /// * `Minimal` - Only write the number of reads in each file
    /// * `Default` - Write all the essential summary statistics
    /// * `Complete` - Write all the summary statistics
    pub mode: &'a SummaryMode,
    pub output: &'a Path,
}

impl<'a> ReadSummaryHandler<'a> {
    pub fn new(
        inputs: &'a mut [PathBuf],
        input_fmt: &'a SeqReadFmt,
        mode: &'a SummaryMode,
        output: &'a Path,
    ) -> Self {
        Self {
            inputs,
            input_fmt,
            mode,
            output,
        }
    }

    pub fn summarize(&mut self) {
        let spin = set_spinner();
        spin.set_message("Calculating summary of fastq files");
        match self.mode {
            SummaryMode::Minimal => {
                let counts = self.par_summarize_minimal();
                let writer = ReadSummaryWriter::new(self.output);
                spin.set_message("Writing records\n");
                writer
                    .write_read_count_only(&counts)
                    .expect("Failed writing to file");
            }
            SummaryMode::Default => {
                let mut records = self.par_summarize_default();
                self.write_records(&spin, &mut records);
            }
            SummaryMode::Complete => {
                let mut records = self.par_summarize_complete();
                self.write_records(&spin, &mut records);
            }
        }
        spin.finish_with_message("Finished processing fastq files\n");
        self.print_output_info();
    }

    fn write_records(&mut self, spin: &ProgressBar, records: &mut [FastqSummary]) {
        // Sort records by file name
        records.sort_by(|a, b| a.path.cmp(&b.path));
        spin.set_message("Writing records\n");
        let writer = ReadSummaryWriter::new(self.output);
        writer.write(records).expect("Failed writing to file");
    }

    fn par_summarize_default(&self) -> Vec<FastqSummary> {
        let (sender, receiver) = channel();

        self.inputs.par_iter().for_each_with(sender, |s, p| {
            let record = self.summarize_default(p);
            s.send(record)
                .expect("Failed parallel processing fastq files");
        });

        receiver.iter().collect()
    }

    fn par_summarize_complete(&self) -> Vec<FastqSummary> {
        let (sender, receiver) = channel();

        self.inputs.par_iter().for_each_with(sender, |s, p| {
            let record = self.summarize_complete(p);
            s.send(record)
                .expect("Failed parallel processing fastq files");
        });

        receiver.iter().collect()
    }

    fn par_summarize_minimal(&self) -> BTreeMap<String, usize> {
        let (sender, receiver) = channel();

        self.inputs.par_iter().for_each_with(sender, |s, p| {
            let count = self.summarize_minimal(p, self.input_fmt);
            s.send((p.display().to_string(), count))
                .expect("Failed parallel processing fastq files");
        });

        receiver.iter().collect()
    }

    fn summarize_minimal(&self, p: &Path, input_fmt: &SeqReadFmt) -> usize {
        match input_fmt {
            SeqReadFmt::Fastq => {
                let mut buff = files::open_file(p);
                fastq::summarize_minimal(&mut buff)
            }
            SeqReadFmt::Gzip => {
                let mut decoder = files::decode_gzip(p);
                fastq::summarize_minimal(&mut decoder)
            }
            SeqReadFmt::Auto => {
                let input_fmt = infer_raw_input_auto(p);
                self.summarize_minimal(p, &input_fmt)
            }
        }
    }

    fn summarize_default(&self, path: &Path) -> FastqSummary {
        let mut summary = FastqSummary::new(path);
        summary.summarize(self.input_fmt);
        summary
    }

    fn summarize_complete(&self, path: &Path) -> FastqSummary {
        let mut summary = FastqSummary::new(path);
        let mapped_records = summary.summarize_map(self.input_fmt);
        let writer = ReadSummaryWriter::new(self.output);
        writer.write_per_read_records(path, &mapped_records.reads, &mapped_records.qscores);
        summary
    }

    fn print_output_info(&self) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Dir", self.output.display());
    }
}
