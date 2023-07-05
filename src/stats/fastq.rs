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

macro_rules! update_records {
    ($self: ident, $record: ident, $sequence: ident, $qrecord: ident) => {
        let $sequence = $record.sequence();
        let $qrecord = $record.quality_scores();
        $self.reads.summarize($sequence);
        let qscores = $self.parse_qscores($qrecord);
        $self.qscores.summarize(&qscores);
    };
}

pub fn count_reads<R: BufRead>(buff: &mut R) -> usize {
    let mut reader = Reader::new(buff);
    reader.records().count()
}

pub struct FastqSummary {
    /// Input file path
    pub path: PathBuf,
    pub reads: ReadRecord,
    pub qscores: ReadQScore,
}

impl FastqSummary {
    pub fn new(path: &Path) -> Self {
        Self {
            path: path.to_path_buf(),
            reads: ReadRecord::new(),
            qscores: ReadQScore::new(),
        }
    }

    pub fn summarize(&mut self, file_fmt: &SeqReadFmt) {
        self.parse_file(file_fmt);
    }

    fn parse_file(&mut self, file_fmt: &SeqReadFmt) {
        match file_fmt {
            SeqReadFmt::Fastq => {
                let mut buff = files::open_file(&self.path);
                self.count(&mut buff)
            }
            SeqReadFmt::Gzip => {
                let mut decoder = files::decode_gzip(&self.path);
                self.count(&mut decoder)
            }
            SeqReadFmt::Auto => {
                let file_fmt = infer_raw_input_auto(&self.path);
                self.parse_file(&file_fmt)
            }
        }
    }

    fn count<R: BufRead>(&mut self, buff: &mut R) {
        let mut reader = Reader::new(buff);
        reader.records().for_each(|r| match r {
            Ok(record) => {
                update_records!(self, record, sequence, qrecord);
            }
            Err(e) => {
                log::error!("Error parsing fastq record: {}", e);
            }
        });
    }

    pub fn compute<R: BufRead>(&mut self, buff: &mut R) {
        let mut reader = Reader::new(buff);
        reader.records().for_each(|r| match r {
            Ok(record) => {
                update_records!(self, record, sequence, qrecord);
            }
            Err(e) => {
                log::error!("Error parsing fastq record: {}", e);
            }
        });
    }

    pub fn compute_mapped<R: BufRead>(&mut self, buff: &mut R) -> FastqMappedRead {
        let mut reader = Reader::new(buff);
        let mut map_records = FastqMappedRead::new();

        reader.records().for_each(|r| match r {
            Ok(record) => {
                update_records!(self, record, sequence, qrecord);
                // Map reads to their index
                self.map_reads(&mut map_records, sequence);
                // Map quality scores to their index
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
        let count = count_reads(&mut buff);
        assert_eq!(count, 2);
    }
}
