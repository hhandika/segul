// A module for sequence statistics.
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;

use colored::Colorize;
use indexmap::IndexSet;
use rayon::prelude::*;

use crate::helper::finder::IDs;
use crate::helper::sequence::SeqParser;
use crate::helper::types::{DataType, InputFmt, TaxonRecords};
use crate::helper::utils;
use crate::stats::sequence::{CharMatrix, CharSummary, Completeness, SiteSummary, Sites, Taxa};
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
        let mut stats: Vec<(Sites, CharMatrix, Taxa)> = self.par_get_stats(files);
        stats.sort_by(|a, b| alphanumeric_sort::compare_path(&a.0.path, &b.0.path));
        let (sites, dna, complete) = self.summarize_char_matrix(&stats);
        let taxon_records = self.summarize_taxa(&ids, &stats);
        spin.finish_with_message("Finished computing summary stats!\n");
        let sum = SummaryWriter::new(&sites, &dna, &complete, self.datatype);
        sum.print_summary().expect("Failed writing to stdout");
        let csv = CsvWriter::new(self.output, prefix, self.datatype);
        csv.write_taxon_summary(&taxon_records)
            .expect("Failed writing to csv");
        csv.write_locus_summary(&stats)
            .expect("Failed writing a per locus csv file");
    }

    pub fn summarize_locus(&mut self, files: &[PathBuf], prefix: &Option<String>) {
        self.check_datatype();
        let spin = utils::set_spinner();
        spin.set_message("Computing per locus summary...");
        files.iter().for_each(|file| {
            let (matrix, _) = SeqParser::new(file, self.datatype).get_alignment(self.input_fmt);
            let mut taxa = Taxa::new();
            taxa.summarize_taxa(&matrix, self.datatype);
            let csv = CsvWriter::new(self.output, prefix, self.datatype);
            csv.write_per_locus_summary(file, &taxa)
                .expect("Failed writing a taxon stats file");
        });
        spin.finish_with_message("Finished computing per locus summary!\n");
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Output dir", self.output.display());
    }

    fn summarize_taxa(
        &self,
        ids: &IndexSet<String>,
        stats: &[(Sites, CharMatrix, Taxa)],
    ) -> BTreeMap<String, TaxonRecords> {
        let mut taxon_summary: BTreeMap<String, TaxonRecords> = BTreeMap::new();
        stats.iter().for_each(|(_, _, taxa)| {
            ids.iter().for_each(|id| {
                if let Some(char_counts) = taxa.records.get(id) {
                    match taxon_summary.get_mut(id) {
                        Some(taxon) => {
                            char_counts.chars.iter().for_each(|(c, count)| {
                                *taxon.char_counts.entry(*c).or_insert(0) += count;
                            });
                            taxon.locus_counts += 1;
                            taxon.total_chars += char_counts.total_chars;
                            taxon.missing_data += char_counts.missing_data;
                            if DataType::Dna == *self.datatype {
                                taxon.gc_count += char_counts.gc_count;
                                taxon.at_count += char_counts.at_count;
                                taxon.nucleotides += char_counts.nucleotides;
                            }
                        }
                        None => {
                            let mut taxon = TaxonRecords::new();
                            taxon.char_counts = char_counts.chars.clone();
                            taxon.locus_counts = 1;
                            taxon.total_chars = char_counts.total_chars;
                            taxon.missing_data = char_counts.missing_data;
                            if DataType::Dna == *self.datatype {
                                taxon.gc_count = char_counts.gc_count;
                                taxon.at_count = char_counts.at_count;
                                taxon.nucleotides = char_counts.nucleotides;
                            }
                            taxon_summary.insert(id.to_string(), taxon);
                        }
                    }
                }
            });
        });

        taxon_summary
    }

    fn get_id(&mut self, files: &[PathBuf]) -> IndexSet<String> {
        IDs::new(files, self.input_fmt, self.datatype).id_unique()
    }

    fn par_get_stats(&self, files: &[PathBuf]) -> Vec<(Sites, CharMatrix, Taxa)> {
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

    fn get_stats(&self, path: &Path) -> (Sites, CharMatrix, Taxa) {
        let aln = SeqParser::new(path, self.datatype);
        let (matrix, header) = aln.get_alignment(self.input_fmt);
        let mut dna = CharMatrix::new();
        dna.count_chars(&matrix, &header, self.datatype);
        let mut sites = Sites::new(path);
        sites.get_stats(&matrix, self.datatype);
        let mut taxa = Taxa::new();
        taxa.summarize_taxa(&matrix, self.datatype);

        (sites, dna, taxa)
    }

    fn summarize_char_matrix(
        &self,
        stats: &[(Sites, CharMatrix, Taxa)],
    ) -> (SiteSummary, CharSummary, Completeness) {
        let (sites, dna): (Vec<Sites>, Vec<CharMatrix>) =
            stats.par_iter().map(|p| (p.0.clone(), p.1.clone())).unzip();
        let mut sum_sites = SiteSummary::new();
        sum_sites.summarize(&sites);
        let mut sum_dna = CharSummary::new();
        sum_dna.summarize(&dna, self.datatype);
        let mut mat_comp = Completeness::new(&self.ntax, self.interval);
        mat_comp.matrix_completeness(&dna);
        (sum_sites, sum_dna, mat_comp)
    }
}
