//! A handler for summarizing raw sequence data.

use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufRead, BufReader, Result, Write},
    path::{Path, PathBuf},
    sync::mpsc::channel,
};

use colored::Colorize;

use indicatif::ProgressBar;
use rayon::prelude::*;

use crate::{
    helper::{
        files,
        types::{infer_raw_input_auto, RawReadFmt, SummaryMode},
        utils::set_spinner,
    },
    parser::fastq::{self, FastqSummaryParser},
    stats::raw::{FastqRecords, QScoreRecords, ReadQScore, ReadRecord},
    writer::raw::RawSummaryWriter,
};

/// Include support for any compressed or uncompressed fastq files.
pub struct RawSummaryHandler<'a> {
    pub inputs: &'a mut [PathBuf],
    pub input_fmt: &'a RawReadFmt,
    pub mode: &'a SummaryMode,
    pub output: &'a Path,
}

impl<'a> RawSummaryHandler<'a> {
    pub fn new(
        inputs: &'a mut [PathBuf],
        input_fmt: &'a RawReadFmt,
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

    pub fn summarize(&mut self, low_mem: bool) {
        let spin = set_spinner();
        spin.set_message("Calculating summary of fastq files");

        if self.mode == &SummaryMode::Minimal {
            let counts = self.count_read();
            let writer = RawSummaryWriter::new(self.output);
            spin.set_message("Writing records\n");
            writer
                .write_read_count_only(&counts)
                .expect("Failed writing to file");
        } else {
            self.summarize_other_mode(&spin, low_mem);
        }
        spin.finish_with_message("Finished processing fastq files\n");
        self.print_output_info();
    }

    fn summarize_other_mode(&mut self, spin: &ProgressBar, low_mem: bool) {
        if low_mem {
            self.summarize_lowmem()
                .expect("Failed summarizing fastq files");
        } else {
            let mut records = self.par_summarize();
            // Sort records by file name
            records.sort_by(|a, b| a.0.path.cmp(&b.0.path));
            spin.set_message("Writing records\n");
            let writer = RawSummaryWriter::new(self.output);
            writer.write(&records).expect("Failed writing to file");
        }
    }

    /// Use a single tread and write records to file as they are processed
    /// to reduce memory usage.
    pub fn summarize_lowmem(&mut self) -> Result<()> {
        self.inputs.par_sort();
        let handler = RawSummaryWriter::new(self.output);
        let mut writer = handler.write_append();
        self.inputs.iter().for_each(|path| {
            let records = self.parse_fastq(path);
            handler
                .write_records(&mut writer, &[records])
                .expect("Failed writing to file");
        });
        writer.flush()?;

        Ok(())
    }

    fn count_read(&self) -> BTreeMap<String, usize> {
        let (sender, receiver) = channel();

        self.inputs.par_iter().for_each_with(sender, |s, p| {
            let input_fmt = if self.input_fmt == &RawReadFmt::Auto {
                infer_raw_input_auto(p)
            } else {
                self.input_fmt.clone()
            };
            let count = match input_fmt {
                RawReadFmt::Fastq => {
                    let mut buff = files::open_file(p);
                    fastq::count_reads(&mut buff)
                }
                RawReadFmt::Gzip => {
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

    fn par_summarize(&self) -> Vec<(FastqRecords, QScoreRecords)> {
        let (sender, receiver) = channel();

        self.inputs.par_iter().for_each_with(sender, |s, p| {
            let record = self.parse_fastq(p);

            s.send(record)
                .expect("Failed parallel processing fastq files");
        });

        receiver.iter().collect()
    }

    fn parse_fastq(&self, path: &Path) -> (FastqRecords, QScoreRecords) {
        let input_fmt = if self.input_fmt == &RawReadFmt::Auto {
            infer_raw_input_auto(path)
        } else {
            self.input_fmt.clone()
        };
        match input_fmt {
            RawReadFmt::Fastq => {
                let file = File::open(path).expect("Failed opening fastq file");
                let mut buff = BufReader::new(file);
                if self.mode == &SummaryMode::Complete {
                    self.map_record(&mut buff, path)
                } else {
                    self.parse_record(&mut buff, path)
                }
            }
            RawReadFmt::Gzip => {
                let mut decoder = files::decode_gzip(path);
                if self.mode == &SummaryMode::Complete {
                    self.map_record(&mut decoder, path)
                } else {
                    self.parse_record(&mut decoder, path)
                }
            }
            _ => unreachable!("Unsupported input format"),
        }
    }

    fn parse_record<R: BufRead>(&self, buff: &mut R, path: &Path) -> (FastqRecords, QScoreRecords) {
        let mut records = FastqSummaryParser::new();
        records.parse_record(buff);
        self.summarize_records(path, &records.reads, &records.qscores)
    }

    fn map_record<R: BufRead>(&self, buff: &mut R, path: &Path) -> (FastqRecords, QScoreRecords) {
        let mut records = FastqSummaryParser::new();
        let mut mapped_records = records.parse_map_records(buff);
        let writer = RawSummaryWriter::new(self.output);
        writer.write_per_read_records(path, &mapped_records.reads, &mapped_records.qscores);
        mapped_records.reads.clear();
        mapped_records.qscores.clear();
        self.summarize_records(path, &records.reads, &records.qscores)
    }

    fn summarize_records(
        &self,
        path: &Path,
        reads: &[ReadRecord],
        qscores: &[ReadQScore],
    ) -> (FastqRecords, QScoreRecords) {
        let mut seq_records = FastqRecords::new(path);
        let mut q_records = QScoreRecords::new();

        seq_records.summarize(reads);
        q_records.summarize(qscores);

        (seq_records, q_records)
    }

    fn print_output_info(&self) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Dir", self.output.display());
    }
}
