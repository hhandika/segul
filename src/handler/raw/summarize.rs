//! A handler for summarizing raw sequence data.

use std::{path::PathBuf, time::Instant};

use noodles::fastq;

use crate::{
    helper::types::{RawReadFmt, SummaryMode},
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
        let dur = Instant::now();
        self.summarize_fastq();
        let dur = dur.elapsed();
        println!("Elapsed time: {:?}", dur);
    }

    fn summarize_fastq(&self) {
        self.inputs.iter().for_each(|p| {
            let gzip_buff = decode_gzip(p);
            let mut reader = fastq::Reader::new(gzip_buff);
            let mut count: usize = 0;
            for result in reader.records() {
                let _record = result.unwrap();
                count += 1;
            }
            println!("{}: {}", p.display(), count);
        });
    }
}
