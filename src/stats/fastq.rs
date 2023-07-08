//! Fastq parser only for Illumina 1.8+ and Sanger quality scores

use std::{
    collections::BTreeMap,
    io::BufRead,
    path::{Path, PathBuf},
};

use noodles::fastq::Reader;

use crate::{
    helper::{
        files,
        types::{infer_raw_input_auto, SeqReadFmt},
    },
    parser::qscores::QScoreParser,
    stats::read::ReadRecord,
};

use super::{qscores::ReadQScore, read::ReadSummary};

macro_rules! summarize_reads {
    ($self: ident, $record: ident, $sequence: ident) => {
        let $sequence = $record.sequence();
        $self.reads.summarize($sequence);
        $self.read_summary.summarize(&$self.reads);
    };
}

macro_rules! summarize_qscores {
    ($self: ident, $record: ident, $qrecord: ident) => {
        let $qrecord = $record.quality_scores();
        let qscores = $self.parse_qscores($qrecord);
        $self.qscores.summarize(&qscores);
    };
}

macro_rules! compute_stats {
    ($self: ident, $fmt: ident, $compute: ident) => {
        match $fmt {
            SeqReadFmt::Fastq => {
                let mut buff = files::open_file(&$self.path);
                $self.$compute(&mut buff)
            }
            SeqReadFmt::Gzip => {
                let mut decoder = files::decode_gzip(&$self.path);
                $self.$compute(&mut decoder)
            }
            _ => unreachable!("Unsupported input format"),
        }
    };
}

trait SeqRead {
    fn parse_input_fmt(&self, path: &Path, file_fmt: &SeqReadFmt) -> SeqReadFmt {
        if file_fmt == &SeqReadFmt::Auto {
            infer_raw_input_auto(path)
        } else {
            *file_fmt
        }
    }
}

pub struct FastqSummaryMin {
    pub path: PathBuf,
    pub read_count: usize,
}

impl SeqRead for FastqSummaryMin {}

impl FastqSummaryMin {
    pub fn new(path: &Path) -> Self {
        Self {
            read_count: 0,
            path: path.to_path_buf(),
        }
    }

    pub fn summarize(&mut self, file_fmt: &SeqReadFmt) {
        let input_fmt = self.parse_input_fmt(&self.path, file_fmt);
        compute_stats!(self, input_fmt, count_reads);
    }

    fn count_reads<R: BufRead>(&mut self, buff: &mut R) {
        let mut reader = Reader::new(buff);
        self.read_count = reader.records().count()
    }
}

pub struct FastqSummary {
    /// Input file path
    pub path: PathBuf,
    pub reads: ReadRecord,
    pub read_summary: ReadSummary,
    pub qscores: ReadQScore,
}

impl SeqRead for FastqSummary {}

impl FastqSummary {
    pub fn new(path: &Path) -> Self {
        Self {
            path: path.to_path_buf(),
            reads: ReadRecord::new(),
            read_summary: ReadSummary::new(),
            qscores: ReadQScore::new(),
        }
    }

    pub fn summarize(&mut self, file_fmt: &SeqReadFmt) {
        let input_fmt = self.parse_input_fmt(&self.path, file_fmt);
        compute_stats!(self, input_fmt, compute_default);
    }

    pub fn summarize_map(&mut self, file_fmt: &SeqReadFmt) -> FastqMappedRead {
        let input_fmt = self.parse_input_fmt(&self.path, file_fmt);
        compute_stats!(self, input_fmt, compute_mapped)
    }

    fn compute_default<R: BufRead>(&mut self, buff: &mut R) {
        let mut reader = Reader::new(buff);
        reader.records().for_each(|r| match r {
            Ok(record) => {
                summarize_reads!(self, record, sequence);
                summarize_qscores!(self, record, qrecord);
            }
            Err(e) => {
                log::error!("Error parsing fastq record: {}", e);
            }
        });
    }

    fn compute_mapped<R: BufRead>(&mut self, buff: &mut R) -> FastqMappedRead {
        let mut reader = Reader::new(buff);
        let mut map_records = FastqMappedRead::new();

        reader.records().for_each(|r| match r {
            Ok(record) => {
                summarize_reads!(self, record, sequence);
                // Map reads to their index
                self.map_reads(&mut map_records, sequence);
                // Map quality scores to their index
                summarize_qscores!(self, record, qrecord);
                self.map_qscores(&mut map_records, qrecord);
            }
            Err(e) => {
                log::error!("Error parsing fastq record: {}", e);
            }
        });

        map_records
    }

    fn map_reads(&self, records: &mut FastqMappedRead, sequence: &[u8]) {
        let mut index = 1;
        sequence.iter().for_each(|s| {
            if let Some(read) = records.reads.get_mut(&index) {
                read.add(s);
            } else {
                let mut read = ReadRecord::new();
                read.add(s);
                records.reads.insert(index, read);
            }
            index += 1;
        });
    }

    fn map_qscores(&self, records: &mut FastqMappedRead, values: &[u8]) {
        let mut index = 1;
        let qscores = self.parse_qscores(values);
        qscores.iter().for_each(|s| {
            if let Some(qscore) = records.qscores.get_mut(&index) {
                qscore.update(s);
            } else {
                let mut qscore = ReadQScore::new();
                qscore.update(s);
                records.qscores.insert(index, qscore);
            }
            index += 1;
        });
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
}

/// Data structure for storing mapped read records
pub struct FastqMappedRead {
    pub reads: BTreeMap<i32, ReadRecord>,
    pub qscores: BTreeMap<i32, ReadQScore>,
}

impl Default for FastqMappedRead {
    fn default() -> Self {
        Self::new()
    }
}

impl FastqMappedRead {
    /// Create a new FastqMappedRead instance
    pub fn new() -> Self {
        Self {
            reads: BTreeMap::new(),
            qscores: BTreeMap::new(),
        }
    }
}

#[cfg(test)]
mod test {
    use std::path::Path;

    use super::*;

    #[test]
    fn test_read_count() {
        let path = Path::new("tests/files/raw/read_1.fastq");
        let file_fmt = SeqReadFmt::Fastq;
        let mut summary = FastqSummaryMin::new(path);
        summary.summarize(&file_fmt);
        assert_eq!(summary.path, PathBuf::from(path));
        assert_eq!(summary.read_count, 2);
    }
}
