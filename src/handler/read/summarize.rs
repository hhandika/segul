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

        if self.mode == &SummaryMode::Minimal {
            let counts = self.count_read();
            let writer = ReadSummaryWriter::new(self.output);
            spin.set_message("Writing records\n");
            writer
                .write_read_count_only(&counts)
                .expect("Failed writing to file");
        } else {
            self.summarize_other_mode(&spin);
        }
        spin.finish_with_message("Finished processing fastq files\n");
        self.print_output_info();
    }

    fn summarize_other_mode(&mut self, spin: &ProgressBar) {
        let mut records = self.par_summarize();
        // Sort records by file name
        records.sort_by(|a, b| a.path.cmp(&b.path));
        spin.set_message("Writing records\n");
        let writer = ReadSummaryWriter::new(self.output);
        writer.write(&records).expect("Failed writing to file");
    }

    fn count_read(&self) -> BTreeMap<String, usize> {
        let (sender, receiver) = channel();

        self.inputs.par_iter().for_each_with(sender, |s, p| {
            let input_fmt = if self.input_fmt == &SeqReadFmt::Auto {
                infer_raw_input_auto(p)
            } else {
                self.input_fmt.clone()
            };
            let count = match input_fmt {
                SeqReadFmt::Fastq => {
                    let mut buff = files::open_file(p);
                    fastq::count_reads(&mut buff)
                }
                SeqReadFmt::Gzip => {
                    let mut decoder = files::decode_gzip(p);
                    fastq::count_reads(&mut decoder)
                }
                _ => unreachable!("Unsupported input format"),
            };

            s.send((p.display().to_string(), count))
                .expect("Failed parallel processing fastq files");
        });

        receiver.iter().collect()
    }

    fn par_summarize(&self) -> Vec<FastqSummary> {
        let (sender, receiver) = channel();

        self.inputs.par_iter().for_each_with(sender, |s, p| {
            let record = self.parse_fastq(p);

            s.send(record)
                .expect("Failed parallel processing fastq files");
        });

        receiver.iter().collect()
    }

    fn parse_fastq(&self, path: &Path) -> FastqSummary {
        let mut summary = FastqSummary::new(path);
        summary.summarize(self.input_fmt);
        summary
    }

    // fn map_record<R: BufRead>(&self, buff: &mut R, path: &Path) -> (ReadRecord, ReadQScore) {
    //     let mut summary = FastqSummary::new();
    //     let mut mapped_records = summary.compute_mapped(buff);
    //     let writer = ReadSummaryWriter::new(self.output);
    //     writer.write_per_read_records(path, &mapped_records.reads, &mapped_records.qscores);
    //     (mapped_records.reads, mapped_records.qscores)
    // }

    fn print_output_info(&self) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Dir", self.output.display());
    }
}
