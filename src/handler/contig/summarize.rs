use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;

use colored::Colorize;
use rayon::prelude::*;

use crate::writer::contigs::ContigSummaryWriter;
use crate::{
    helper::{types::ContigFmt, utils::set_spinner},
    stats::contigs::ContigSummary,
};

pub struct ContigSummaryHandler<'a> {
    files: &'a [PathBuf],
    input_fmt: &'a ContigFmt,
    output: &'a Path,
}

impl<'a> ContigSummaryHandler<'a> {
    pub fn new(files: &'a [PathBuf], input_fmt: &'a ContigFmt, output: &'a Path) -> Self {
        Self {
            files,
            input_fmt,
            output,
        }
    }

    pub fn summarize(&self) {
        let spin = set_spinner();
        spin.set_message("Calculating summary of contig files");
        let contig_summary = self.summarize_contigs();
        let writer = ContigSummaryWriter::new(&contig_summary, self.output);
        spin.set_message("Writing records\n");
        writer.write().expect("Failed writing to file");
        spin.finish_with_message("Finished processing contig files\n");
        self.print_input_info();
    }

    fn summarize_contigs(&self) -> Vec<ContigSummary> {
        if self.files.is_empty() {
            panic!("No contig files found")
        }

        if self.files.len() == 1 {
            return vec![self.process_contigs(&self.files[0])];
        }

        let (sender, receiver) = channel();

        self.files.par_iter().for_each_with(sender, |s, p| {
            let summary = self.process_contigs(p);
            s.send(summary).expect("Failed sending data");
        });

        receiver.iter().collect()
    }

    fn process_contigs(&self, input: &Path) -> ContigSummary {
        let mut summary = ContigSummary::new();
        summary.summarize(input, self.input_fmt);
        summary
    }

    fn print_input_info(&self) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Dir", self.output.display());
    }
}
