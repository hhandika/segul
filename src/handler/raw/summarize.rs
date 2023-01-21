//! A handler for summarizing raw sequence data.

use std::{
    fs::{self, File},
    io::{BufWriter, Result, Write},
    path::{Path, PathBuf},
    sync::mpsc::channel,
};

use colored::Colorize;
use noodles::fastq::Reader;
use rayon::prelude::*;

use crate::{
    helper::{
        fastq::{FastqRecords, QScoreRecords, ReadQScore, ReadRecord},
        types::{RawReadFmt, SummaryMode},
        utils::set_spinner,
    },
    parser::{fastq::QScoreParser, gzip::decode_gzip},
};

const DEFAULT_OUTPUT: &str = "summary.tsv";

/// Include support for any compressed or uncompressed fastq files.
pub struct RawSummaryHandler<'a> {
    pub inputs: &'a [PathBuf],
    pub input_fmt: &'a RawReadFmt,
    pub mode: &'a SummaryMode,
    pub output: &'a Path,
}

impl<'a> RawSummaryHandler<'a> {
    pub fn new(
        inputs: &'a [PathBuf],
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

    pub fn summarize(&self, low_mem: bool) {
        let spin = set_spinner();
        spin.set_message("Calculating summary of fastq files");

        if low_mem {
            self.summarize_lowmem();
        } else {
            let mut records = self.par_summarize();

            spin.set_message("Writing records\n");
            let mut writer = self.create_output_file(self.output, DEFAULT_OUTPUT);
            self.write_records(&mut writer, &records)
                .expect("Failed writing to file");
            // Sort records by file name
            records.sort_by(|a, b| a.0.path.cmp(&b.0.path));
        }
        spin.finish_with_message("Finished processing fastq files\n");
        self.print_output_info();
    }

    fn summarize_lowmem(&self) {
        self.inputs.iter().for_each(|path| {
            let records = self.parse_record(path);
            let output_dir = self.output.join("tsv");
            let fname = format!(
                "{}_{}",
                path.file_name()
                    .expect("Failed getting file name")
                    .to_str()
                    .expect("Failed converting file name to string"),
                DEFAULT_OUTPUT
            );
            let mut writer = self.create_output_file(&output_dir, &fname);
            self.write_records(&mut writer, &[records])
                .expect("Failed writing to file");
        })
    }

    fn par_summarize(&self) -> Vec<(FastqRecords, QScoreRecords)> {
        let (sender, receiver) = channel();

        self.inputs.par_iter().for_each_with(sender, |s, p| {
            let record = self.parse_record(p);

            s.send(record)
                .expect("Failed parallel processing fastq files");
        });

        receiver.iter().collect()
    }

    fn parse_record(&self, path: &Path) -> (FastqRecords, QScoreRecords) {
        let gzip_buff = decode_gzip(path);
        let mut reader = Reader::new(gzip_buff);
        let mut reads = Vec::new();
        let mut q_records = Vec::new();

        reader.records().for_each(|r| match r {
            Ok(record) => {
                let mut read_records = ReadRecord::new();
                read_records.summarize(record.sequence());
                reads.push(read_records);
                let mut read_qscores = ReadQScore::new();
                let qrecord = record.quality_scores();
                let qscores = self.parse_qscores(qrecord);
                read_qscores.summarize(&qscores);
                q_records.push(read_qscores);
            }
            Err(e) => {
                log::error!("Error parsing fastq record: {}", e);
            }
        });

        self.summarize_records(path, &reads, &q_records)
    }

    fn parse_qscores(&self, qscore: &[u8]) -> Vec<u8> {
        let mut qscores = Vec::with_capacity(qscore.len());
        let parser = QScoreParser::new(qscore);
        parser.into_iter().for_each(|q| {
            if let Some(q) = q {
                qscores.push(q);
            }
        });

        qscores
    }

    fn create_output_file(&self, output_dir: &Path, fname: &str) -> BufWriter<File> {
        fs::create_dir_all(output_dir).expect("Failed to create output directory");
        let fpath = output_dir.join(fname);
        let file = File::create(fpath).expect("Failed to create summary file");
        BufWriter::new(file)
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
        log::info!("{:18}: {}", "Summary file", "summary.tsv")
    }

    fn write_records<W: Write>(
        &self,
        writer: &mut W,
        records: &[(FastqRecords, QScoreRecords)],
    ) -> Result<()> {
        match self.mode {
            SummaryMode::Minimal => {
                writeln!(
                    writer,
                    "File\tNumReads\tNumBases\tMinReadLen\tMeanReadLen\tMaxReadLen"
                )?;
                for (seq, _) in records {
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        seq.path.display(),
                        seq.num_reads,
                        seq.num_bases,
                        seq.min_read_len,
                        seq.mean_read_len,
                        seq.max_read_len
                    )?;
                }
            }
            SummaryMode::Complete => {
                writeln!(writer,"File\tNumReads\tNumBases\tMinReadLen\tMeanReadLen\tMaxReadLen\tLowQ\tMean\tMin\tMax")?;
                for (seq, q) in records {
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        seq.path.display(),
                        seq.num_reads,
                        seq.num_bases,
                        seq.min_read_len,
                        seq.mean_read_len,
                        seq.max_read_len,
                        q.low_q,
                        q.mean,
                        q.min,
                        q.max
                    )?;
                }
            }
            SummaryMode::Default => {
                writeln!(
                    writer,
                    "File\tNumReads\tNumBases\t\
                    MinReadLen\tMeanReadLen\tMaxReadLen\t\
                    GCcount\tGCcontent\tATcount\tATContent\t\
                    Ncount\tNcontent\t\
                    LowQ\tMean\tMin\tMax\
                    "
                )?;
                for (seq, q) in records {
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        seq.path.display(),
                        seq.num_reads,
                        seq.num_bases,
                        seq.min_read_len,
                        seq.mean_read_len,
                        seq.max_read_len,
                        seq.gc_count,
                        seq.gc_content,
                        seq.at_count,
                        seq.at_content,
                        seq.n_count,
                        seq.n_content,
                        q.low_q,
                        q.mean,
                        q.min,
                        q.max,
                    )?;
                }
            }
        }
        Ok(())
    }
}
