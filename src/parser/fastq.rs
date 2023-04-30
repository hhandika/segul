//! Fastq parser only for Illumina 1.8+ and Sanger quality scores

use std::{collections::BTreeMap, io::BufRead};

use noodles::fastq::Reader;

use crate::helper::stats::{QScoreStream, ReadQScore, ReadRecord};

macro_rules! insert_records {
    ($self: ident, $record: ident, $sequence: ident, $qrecord: ident) => {
        let mut read_records = ReadRecord::new();
        let $sequence = $record.sequence();
        let $qrecord = $record.quality_scores();
        read_records.summarize($sequence);
        $self.reads.push(read_records);
        let mut read_qscores = ReadQScore::new();
        let qscores = $self.parse_qscores($qrecord);
        read_qscores.summarize(&qscores);
        $self.qscores.push(read_qscores);
    };
}

pub fn count_reads<R: BufRead>(buff: &mut R) -> usize {
    let mut reader = Reader::new(buff);
    reader.records().count()
}

pub struct FastqSummaryParser {
    /// Input file path
    pub reads: Vec<ReadRecord>,
    pub qscores: Vec<ReadQScore>,
}

impl FastqSummaryParser {
    pub fn new() -> Self {
        Self {
            reads: Vec::new(),
            qscores: Vec::new(),
        }
    }

    pub fn parse_record<R: BufRead>(&mut self, buff: &mut R) {
        let mut reader = Reader::new(buff);
        reader.records().for_each(|r| match r {
            Ok(record) => {
                insert_records!(self, record, sequence, qrecord);
            }
            Err(e) => {
                log::error!("Error parsing fastq record: {}", e);
            }
        });
    }

    pub fn parse_map_records<R: BufRead>(&mut self, buff: &mut R) -> FastqMappedRead {
        let mut reader = Reader::new(buff);
        let mut map_records = FastqMappedRead::new();

        reader.records().for_each(|r| match r {
            Ok(record) => {
                insert_records!(self, record, sequence, qrecord);
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

/// Parse Illumina 1.8+ and Sanger quality scores
pub struct QScoreParser<'a> {
    /// Quality scores in ASCII format
    pub scores: &'a [u8],
    /// Index of the current quality score
    index: usize,
}

impl<'a> QScoreParser<'a> {
    /// Create a new QScoreParser
    pub fn new(scores: &'a [u8]) -> Self {
        Self { scores, index: 0 }
    }
}

impl<'a> Iterator for QScoreParser<'a> {
    type Item = Option<u8>;
    /// Read ASCII from vector bytes
    /// and convert to Illumina 1.8+ and Sanger quality scores
    fn next(&mut self) -> Option<Self::Item> {
        let q = self.scores.get(self.index);
        match q {
            Some(q) => {
                if q > &74 {
                    panic!("Invalid quality score: {}", q);
                }
                self.index += 1;
                Some(Some(q - 33))
            }
            None => None, // End iterator when index is out of bound
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    macro_rules! qscore_parser {
        ($scores:expr, $sum: ident) => {
            let records = QScoreParser::new($scores);
            let $sum: u8 = records
                .into_iter()
                .map(|x| match x {
                    Some(x) => x,
                    None => 0,
                })
                .sum();
        };
    }

    #[test]
    fn test_qscore_parser() {
        let scores = b"II";
        let scores_2 = b"00";
        qscore_parser!(scores, sum);
        qscore_parser!(scores_2, sum2);
        assert_eq!(80, sum);
        assert_eq!(30, sum2);
    }

    #[test]
    #[should_panic(expected = "Invalid quality score: 75")]
    fn test_qscore_parser_panic() {
        let scores = b"II!)K";
        let q = QScoreParser::new(scores);
        q.into_iter().for_each(|x| match x {
            Some(x) => println!("{}", x),
            None => println!("None"),
        });
    }

    #[test]
    fn test_iter_empty() {
        let scores = b"";
        let q = QScoreParser::new(scores);
        q.into_iter().for_each(|x| match x {
            Some(x) => println!("{}", x),
            None => println!("None"),
        });
    }
}
