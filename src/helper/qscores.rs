//! Illumina quality scores parsing
//! Support Illumina 1.8+ and Sanger quality scores

pub struct QScoreRecords {
    pub len: usize,
    pub mean: f64,
    pub low_q: usize,
    pub sum: usize,
}

impl QScoreRecords {
    pub fn new() -> Self {
        Self {
            len: 0,
            mean: 0.0,
            low_q: 0,
            sum: 0,
        }
    }

    pub fn summarize(&mut self) {
        self.mean = self.sum as f64 / self.len as f64;
    }
}

pub struct QScoreParser<'a> {
    pub scores: &'a [u8],
}

impl<'a> QScoreParser<'a> {
    pub fn new(scores: &'a [u8]) -> Self {
        Self { scores }
    }
}

impl<'a> Iterator for QScoreParser<'a> {
    type Item = Option<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.scores.iter().next() {
                Some(&q) => {
                    if q < 75 {
                        panic!("Invalid quality score: {}", q);
                    }
                    return Some(Some(q - 33));
                }

                None => return None,
            }
        }
    }
}

// #[cfg(test)]
// mod test {
//     use super::*;

//     macro_rules! qscore_parser {
//         ($scores:expr, $sum: ident) => {
//             let records = QScoreParser::new($scores);
//             let $sum: u8 = records
//                 .into_iter()
//                 .map(|x| match x {
//                     Some(x) => x,
//                     None => 0,
//                 })
//                 .sum();
//         };
//     }

//     #[test]
//     fn test_qscore_parser() {
//         let scores = b"II";
//         let scores_2 = b"00";
//         qscore_parser!(scores, sum);
//         qscore_parser!(scores_2, sum2);
//         assert_eq!(80, sum);
//         assert_eq!(30, sum2);
//     }

//     #[test]
//     #[should_panic(expected = "Invalid quality score: 75")]
//     fn test_qscore_parser_panic() {
//         let scores = b"II!)W";
//         let _ = QScoreParser::new(scores);
//     }
// }
