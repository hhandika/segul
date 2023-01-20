//! A handler for summarizing raw sequence data.

use std::{
    path::PathBuf,
    sync::{
        atomic::{AtomicUsize, Ordering},
        mpsc::channel,
        Arc, Mutex,
    },
    time::Instant,
};

use indicatif::ProgressBar;
use noodles::fastq;
use rayon::prelude::*;

use crate::{
    helper::{
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
        let dur = Instant::now();
        self.summarize_fastq();
        let dur = dur.elapsed();
        println!("Elapsed time: {:?}", dur);
        let dur = Instant::now();
        self.summarize_par_fastq();
        let dur = dur.elapsed();
        println!("Elapsed time: {:?}", dur);
        let dur = Instant::now();
        self.par_sender();
        let dur = dur.elapsed();
        println!("Elapsed time: {:?}", dur);
    }

    fn summarize_par_fastq(&self) {
        let spin = set_spinner();
        spin.set_message("Reading fastq records");
        let fcount = AtomicUsize::new(0);
        let records: Vec<(PathBuf, usize)> = Vec::new();
        // Insert record count here
        let mutex_vec = Arc::new(Mutex::new(records));
        self.inputs.par_iter().for_each(|p| {
            let gzip_buff = decode_gzip(p);
            let mut reader = fastq::Reader::new(gzip_buff);
            spin.set_message(format!(
                "Reading fastq {} of {} files",
                fcount.load(Ordering::Relaxed),
                self.inputs.len()
            ));
            let mut data = mutex_vec.lock().unwrap();
            let mut count = 0;
            for result in reader.records() {
                let _record = result.unwrap();
                count += 1;
            }
            data.push((p.to_path_buf(), count));
            fcount.fetch_add(1, Ordering::Relaxed);
        });
        spin.finish_with_message("Finished reading fastq records");
        println!("Count {}", fcount.load(Ordering::Relaxed));
        let data = mutex_vec.lock().unwrap();
        for (p, c) in data.iter() {
            println!("{}: {}", p.display(), c);
        }
    }

    fn par_sender(&self) {
        let (sender, receiver) = channel();
        let progress = ProgressBar::new(self.inputs.len() as u64);
        self.inputs.par_iter().for_each_with(sender, |s, p| {
            let gzip_buff = decode_gzip(p);
            let mut reader = fastq::Reader::new(gzip_buff);
            let mut count = 0;
            for result in reader.records() {
                let _record = result.unwrap();
                count += 1;
            }
            progress.inc(1);
            s.send((p.to_path_buf(), count)).unwrap();
        });
        for (p, c) in receiver {
            println!("{}: {}", p.display(), c);
        }
    }

    fn summarize_fastq(&self) {
        let spin = set_spinner();
        spin.set_message("Reading fastq records");
        let mut fcount = 0;
        let mut records: Vec<(PathBuf, usize)> = Vec::new();
        for p in self.inputs.iter() {
            let gzip_buff = decode_gzip(p);
            let mut reader = fastq::Reader::new(gzip_buff);
            spin.set_message(format!(
                "Reading fastq {} of {} files",
                fcount,
                self.inputs.len()
            ));
            let mut count = 0;
            for result in reader.records() {
                let _record = result.unwrap();
                count += 1;
            }
            records.push((p.to_path_buf(), count));
            fcount += 1;
        }
        spin.finish_with_message("Finished reading fastq records");
        println!("Count {}", fcount);
        for (p, c) in records.iter() {
            println!("{}: {}", p.display(), c);
        }
    }
}
