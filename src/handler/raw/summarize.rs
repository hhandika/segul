//! A handler for summarizing raw sequence data.

use std::{
    path::{Path, PathBuf},
    sync::mpsc::channel,
};

use noodles::fastq;
use rayon::prelude::*;

use crate::{
    helper::{
        fastq::{FastqRecords, QScoreRecords, ReadQScore, ReadRecord},
        types::{RawReadFmt, SummaryMode},
        utils::set_spinner,
    },
    parser::{fastq::QScoreParser, gzip::decode_gzip},
};

/// Include support for any compressed or uncompressed fastq files.
pub struct RawSummaryHandler<'a> {
    pub inputs: &'a [PathBuf],
    pub input_fmt: &'a RawReadFmt,
    pub mode: &'a SummaryMode,
}

impl<'a> RawSummaryHandler<'a> {
    pub fn new(inputs: &'a [PathBuf], input_fmt: &'a RawReadFmt, mode: &'a SummaryMode) -> Self {
        Self {
            inputs,
            input_fmt,
            mode,
        }
    }

    pub fn summarize(&self) {
        let spin = set_spinner();
        spin.set_message("Calculating summary of fastq files");
        let records = self.par_summarize();
        spin.finish_with_message("Finished processing fastq files\n");
        self.print_records(&records);
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
        let mut reader = fastq::Reader::new(gzip_buff);
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

    fn print_records(&self, records: &[(FastqRecords, QScoreRecords)]) {
        match self.mode {
            SummaryMode::Minimal => {
                println!("File\tNumReads\tNumBases\tMinReadLen\tMeanReadLen\tMaxReadLen");
                for (seq, _) in records {
                    println!(
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        seq.path.display(),
                        seq.num_reads,
                        seq.num_bases,
                        seq.min_read_len,
                        seq.mean_read_len,
                        seq.max_read_len
                    );
                }
            }
            SummaryMode::Complete => {
                println!("File\tNumReads\tNumBases\tMinReadLen\tMeanReadLen\tMaxReadLen\tLowQ\tSum\tMean\tMin\tMax");
                for (seq, q) in records {
                    println!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        seq.path.display(),
                        seq.num_reads,
                        seq.num_bases,
                        seq.min_read_len,
                        seq.mean_read_len,
                        seq.max_read_len,
                        q.low_q,
                        q.sum,
                        q.mean,
                        q.min,
                        q.max
                    );
                }
            }
            SummaryMode::Default => {
                println!("File\tNumReads\tNumBases\tMinReadLen\tMeanReadLen\tMaxReadLen\tLowQ\tSum\tMean\tMin\tMax");
                for (seq, q) in records {
                    println!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        seq.path.display(),
                        seq.num_reads,
                        seq.num_bases,
                        seq.min_read_len,
                        seq.mean_read_len,
                        seq.max_read_len,
                        q.low_q,
                        q.sum,
                        q.mean,
                        q.min,
                        q.max
                    );
                }
            }
        }
    }
}
