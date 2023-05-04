use std::io::BufRead;

use super::fasta::FastaReader;

pub struct ContigSummaryParser {
    pub g_count: usize,
    pub c_count: usize,
    pub a_count: usize,
    pub t_count: usize,
    pub n_count: usize,
    pub contig_len: Vec<usize>,
}

impl ContigSummaryParser {
    pub fn new() -> Self {
        Self {
            g_count: 0,
            c_count: 0,
            a_count: 0,
            t_count: 0,
            n_count: 0,
            contig_len: Vec::new(),
        }
    }

    pub fn parse<R: BufRead>(&mut self, buff: &mut R) {
        let reader = FastaReader::new(buff);

        reader.into_iter().for_each(|r| {
            self.contig_len.push(r.seq.len());
            self.g_count += r.seq.matches('G').count();
            self.c_count += r.seq.matches('C').count();
            self.a_count += r.seq.matches('A').count();
            self.t_count += r.seq.matches('T').count();
            self.n_count += r.seq.matches('N').count();
        })
    }
}
