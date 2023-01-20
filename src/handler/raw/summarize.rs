//! A handler for summarizing raw sequence data.

use std::{
    path::{Path, PathBuf},
    sync::mpsc::channel,
};

use noodles::fastq;
use rayon::prelude::*;

use crate::{
    helper::{
        qscores::QScoreParser,
        types::{RawReadFmt, SummaryMode},
        utils::set_spinner,
    },
    parser::gzip::decode_gzip,
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
        let q = b"IIKW";
        let qrecs = QScoreParser::new(q);
        qrecs.into_iter().for_each(|q| {
            println!("{:?}", q);
        });
        records.iter().for_each(|(p, c)| {
            println!("{}: {}", p.display(), c);
        });
    }

    fn par_summarize(&self) -> Vec<(PathBuf, usize)> {
        let (sender, receiver) = channel();

        self.inputs.par_iter().for_each_with(sender, |s, p| {
            let record = self.parse_record(p);

            s.send((p.to_path_buf(), record))
                .expect("Failed parallel processing fastq files");
        });

        receiver.iter().collect()
    }

    fn parse_record(&self, path: &Path) -> usize {
        let gzip_buff = decode_gzip(path);
        let mut reader = fastq::Reader::new(gzip_buff);
        let mut count = 0;
        reader.records().for_each(|_| count += 1);
        count
    }
}
