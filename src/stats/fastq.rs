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
    stats::read::{QScoreStream, ReadQScore, ReadRecord},
};

use super::read::ReadSummary;

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

pub fn summarize_minimal<R: BufRead>(buff: &mut R) -> usize {
    let mut reader = Reader::new(buff);
    reader.records().count()
}

pub struct FastqSummary {
    /// Input file path
    pub path: PathBuf,
    pub reads: ReadRecord,
    pub read_summary: ReadSummary,
    pub qscores: ReadQScore,
}

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
        let input_fmt = self.parse_input_fmt(file_fmt);
        compute_stats!(self, input_fmt, count);
    }

    pub fn summarize_map(&mut self, file_fmt: &SeqReadFmt) -> FastqMappedRead {
        let input_fmt = self.parse_input_fmt(file_fmt);
        compute_stats!(self, input_fmt, compute_mapped)
    }

    fn parse_input_fmt(&self, file_fmt: &SeqReadFmt) -> SeqReadFmt {
        if file_fmt == &SeqReadFmt::Auto {
            infer_raw_input_auto(&self.path)
        } else {
            file_fmt.to_owned()
        }
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

    fn count<R: BufRead>(&mut self, buff: &mut R) {
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

    fn map_qscores(&self, records: &mut FastqMappedRead, qscores: &[u8]) {
        let mut index = 1;
        qscores.iter().for_each(|s| {
            if let Some(qscore) = records.qscores.get_mut(&index) {
                qscore.update(s);
            } else {
                let mut qscore = QScoreStream::new();
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

pub struct FastqMappedRead {
    pub reads: BTreeMap<i32, ReadRecord>,
    pub qscores: BTreeMap<i32, QScoreStream>,
}

impl FastqMappedRead {
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

    use crate::helper::files;

    use super::*;

    #[test]
    fn test_read_count() {
        let path = Path::new("tests/files/raw/read_1.fastq");
        let mut buff = files::open_file(path);
        let count = summarize_minimal(&mut buff);
        assert_eq!(count, 2);
    }
}
