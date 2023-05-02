pub struct ContigSummaryParser {
    pub contigs: Vec<ContigSummary>,
}

impl ContigSummaryParser {
    pub fn new() -> Self {
        Self {
            contigs: Vec::new(),
        }
    }
}
