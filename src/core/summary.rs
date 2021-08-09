// A module for sequence statistics.
// use std::collections::BTreeMap;

use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;

use rayon::prelude::*;

use crate::helper::finder::IDs;
use crate::helper::sequence::Sequence;
use crate::helper::stats::{CharSummary, Chars, Completeness, SiteSummary, Sites};
use crate::helper::types::{DataType, InputFmt};
use crate::helper::utils;
use crate::writer::summary::{self, CsvWriter, SummaryWriter};

pub struct SeqStats<'a> {
    input_format: &'a InputFmt,
    output: &'a Path,
    datatype: &'a DataType,
    ntax: usize,
    interval: usize,
}

impl<'a> SeqStats<'a> {
    pub fn new(
        input_format: &'a InputFmt,
        output: &'a Path,
        interval: usize,
        datatype: &'a DataType,
    ) -> Self {
        Self {
            input_format,
            output,
            ntax: 0,
            interval,
            datatype,
        }
    }

    pub fn get_seq_stats_file(&mut self, path: &Path) {
        self.check_datatype();
        let spin = utils::set_spinner();
        spin.set_message("Getting alignments...");
        let (site, dna) = self.get_stats(path);
        spin.finish_with_message("Finished getting alignments!\n");
        CsvWriter::new(self.output, self.datatype)
            .write_summary_file(&site, &dna)
            .expect("CANNOT WRITE PER LOCUS SUMMARY STATS");
        summary::print_stats(&site, &dna);
    }

    pub fn get_stats_dir(&mut self, files: &[PathBuf]) {
        self.check_datatype();
        let spin = utils::set_spinner();
        spin.set_message("Indexing alignments...");
        self.get_ntax(files);
        spin.set_message("Computing summary stats...");
        let mut stats: Vec<(Sites, Chars)> = self.par_get_stats(files);
        stats.sort_by(|a, b| alphanumeric_sort::compare_path(&a.0.path, &b.0.path));
        let (sites, dna, complete) = self.get_summary_dna(&stats);
        spin.finish_with_message("Finished computing summary stats!\n");
        let sum = SummaryWriter::new(&sites, &dna, &complete, self.datatype);
        sum.print_summary().expect("Failed writing to stdout");
        CsvWriter::new(self.output, self.datatype)
            .write_summary_dir(&stats)
            .expect("Failed writing a per locus csv file");
    }

    fn get_ntax(&mut self, files: &[PathBuf]) {
        self.ntax = IDs::new(files, self.input_format, self.datatype)
            .get_id_unique()
            .len();
    }

    fn par_get_stats(&self, files: &[PathBuf]) -> Vec<(Sites, Chars)> {
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

    fn get_stats(&self, path: &Path) -> (Sites, Chars) {
        let aln = Sequence::new(path, self.datatype);
        let (matrix, header) = aln.get_alignment(self.input_format);
        let mut dna = Chars::new();
        dna.count_chars(&matrix, &header);
        let mut sites = Sites::new();
        sites.get_stats(path, &matrix, self.datatype);

        (sites, dna)
    }

    fn get_summary_dna(
        &self,
        stats: &[(Sites, Chars)],
    ) -> (SiteSummary, CharSummary, Completeness) {
        let (sites, dna): (Vec<Sites>, Vec<Chars>) =
            stats.par_iter().map(|p| (p.0.clone(), p.1.clone())).unzip();
        let mut sum_sites = SiteSummary::new();
        sum_sites.get_summary(&sites);
        let mut sum_dna = CharSummary::new();
        sum_dna.get_summary(&dna);
        let mut mat_comp = Completeness::new(&self.ntax, self.interval);
        mat_comp.matrix_completeness(&dna);
        (sum_sites, sum_dna, mat_comp)
    }
}
