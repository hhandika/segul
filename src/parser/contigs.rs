use noodles::fasta::Reader;

pub struct ContigSummaryParser {
    pub contigs: Vec<ContigSummary>,
}

impl ContigSummaryParser {
    pub fn new() -> Self {
        Self {
            contigs: Vec::new(),
        }
    }

    pub fn parse<R: BufRead>(&mut self, buff: &mut R) {
        let mut reader = Reader::new(buff);
        reader.records().for_each(|r| match r {
            Ok(record) => {
                let contig = ContigSummary::new(record.id().to_string(), record.sequence().len());
                self.contigs.push(contig);
            }
            Err(e) => {
                log::error!("Error parsing fasta record: {}", e);
            }
        });
    }
}
