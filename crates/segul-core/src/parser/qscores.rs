//! Parse Illumina 1.8+ and Sanger quality scores
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

impl Iterator for QScoreParser<'_> {
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
        ($scores:expr_2021, $sum: ident) => {
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
