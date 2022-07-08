// A module for sequence statistics.

use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;

use indexmap::IndexSet;
use rayon::prelude::*;

use crate::helper::finder::IDs;
use crate::helper::sequence::SeqParser;
use crate::helper::stats::{CharSummary, Chars, Completeness, SiteSummary, Sites, Taxa};
use crate::helper::types::{DataType, InputFmt};
use crate::helper::utils;
use crate::writer::summary::{CsvWriter, SummaryWriter};

pub struct SeqStats<'a> {
    input_fmt: &'a InputFmt,
    output: &'a Path,
    datatype: &'a DataType,
    ntax: usize,
    interval: usize,
}

impl<'a> SeqStats<'a> {
    pub fn new(
        input_fmt: &'a InputFmt,
        output: &'a Path,
        interval: usize,
        datatype: &'a DataType,
    ) -> Self {
        Self {
            input_fmt,
            output,
            ntax: 0,
            interval,
            datatype,
        }
    }

    pub fn summarize_all(&mut self, files: &[PathBuf], prefix: &Option<String>) {
        self.check_datatype();
        let spin = utils::set_spinner();
        spin.set_message("Indexing alignments...");
        let ids = self.get_id(files);
        self.ntax = ids.len();
        spin.set_message("Computing alignment stats...");
        let mut stats: Vec<(Sites, Chars, Taxa)> = self.par_get_stats(files);
        stats.sort_by(|a, b| alphanumeric_sort::compare_path(&a.0.path, &b.0.path));
        let (sites, dna, complete) = self.summarize_dna(&stats);
        spin.finish_with_message("Finished computing summary stats!\n");
        let sum = SummaryWriter::new(&sites, &dna, &complete, self.datatype);
        sum.print_summary().expect("Failed writing to stdout");
        let csv = CsvWriter::new(self.output, prefix, self.datatype, &stats);
        csv.write_taxon_summary(&ids)
            .expect("Failed writing a taxon stats file");
        csv.write_locus_summary()
            .expect("Failed writing a per locus csv file");
    }

    fn get_id(&mut self, files: &[PathBuf]) -> IndexSet<String> {
        IDs::new(files, self.input_fmt, self.datatype).id_unique()
    }

    fn par_get_stats(&self, files: &[PathBuf]) -> Vec<(Sites, Chars, Taxa)> {
        let (send, rec) = channel();
        files.par_iter().for_each_with(send, |s, file| {
            s.send(self.get_stats(file)).unwrap();
        });
        rec.iter().collect()
    }

    fn check_datatype(&mut self) {
        if let DataType::Ignore = self.datatype {
            self.datatype = &DataType::Dna
        }
    }

    fn get_stats(&self, path: &Path) -> (Sites, Chars, Taxa) {
        let aln = SeqParser::new(path, self.datatype);
        let (matrix, header) = aln.get_alignment(self.input_fmt);
        let mut dna = Chars::new();
        dna.count_chars(&matrix, &header);
        let mut sites = Sites::new();
        sites.get_stats(path, &matrix, self.datatype);
        let mut taxa = Taxa::new();
        taxa.summarize_taxa(&matrix);

        (sites, dna, taxa)
    }

    fn summarize_dna(
        &self,
        stats: &[(Sites, Chars, Taxa)],
    ) -> (SiteSummary, CharSummary, Completeness) {
        let (sites, dna): (Vec<Sites>, Vec<Chars>) =
            stats.par_iter().map(|p| (p.0.clone(), p.1.clone())).unzip();
        let mut sum_sites = SiteSummary::new();
        sum_sites.summarize(&sites);
        let mut sum_dna = CharSummary::new();
        sum_dna.summarize(&dna);
        let mut mat_comp = Completeness::new(&self.ntax, self.interval);
        mat_comp.matrix_completeness(&dna);
        (sum_sites, sum_dna, mat_comp)
    }
}
