// A module for sequence statistics.
use std::collections::BTreeMap;

use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;

use indexmap::IndexSet;
use rayon::prelude::*;

use crate::helper::finder::IDs;
use crate::helper::sequence::Sequence;
use crate::helper::stats::{CharSummary, Chars, Completeness, SiteSummary, Sites};
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

    #[allow(unused_variables)]
    pub fn get_stats_dir(&mut self, files: &[PathBuf]) {
        self.check_datatype();
        let spin = utils::set_spinner();
        spin.set_message("Indexing alignments...");
        let ids = self.get_id(files);
        self.ntax = ids.len();
        spin.set_message("Computing summary stats...");
        let mut stats: Vec<(Sites, Chars)> = self.par_get_stats(files);
        stats.sort_by(|a, b| alphanumeric_sort::compare_path(&a.0.path, &b.0.path));
        let (sites, dna, complete) = self.get_summary_dna(&stats);
        let taxon_stats = self.compute_taxon_stats(files, &ids);
        spin.finish_with_message("Finished computing summary stats!\n");
        let sum = SummaryWriter::new(&sites, &dna, &complete, self.datatype);
        sum.print_summary().expect("Failed writing to stdout");
        CsvWriter::new(self.output, self.datatype)
            .write_summary_dir(&stats)
            .expect("Failed writing a per locus csv file");
    }

    fn count_sites(&self, seq: &str) -> usize {
        let mut seq = seq.to_string();
        seq.retain(|c| !"?-".contains(c));
        seq.len()
    }

    fn compute_taxon_stats(
        &self,
        locus_files: &[PathBuf],
        ids: &IndexSet<String>,
    ) -> BTreeMap<String, Vec<usize>> {
        let mut taxon_stats: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        locus_files.iter().for_each(|file| {
            ids.iter().for_each(|id| {
                let (seq, _) = Sequence::new(file, self.datatype).get(self.input_fmt);
                match seq.get(id) {
                    Some(seq) => {
                        let site_counts = self.count_sites(seq);
                        taxon_stats
                            .entry(id.to_string())
                            .or_insert(vec![])
                            .push(site_counts);
                    }
                    None => (),
                }
            })
        });
        taxon_stats
    }

    fn get_id(&mut self, files: &[PathBuf]) -> IndexSet<String> {
        IDs::new(files, self.input_fmt, self.datatype).get_id_unique()
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
        let (matrix, header) = aln.get_alignment(self.input_fmt);
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

// #[cfg(test)]
// mod test {
//     use super::*;

//     const INPUT_FMT: InputFmt = InputFmt::Fasta;
//     const DATATYPE: DataType = DataType::Dna;

//     #[test]
//     fn test_site_counts() {
//         let seq = "AGTCT-?";
//         let id = Id::new(Path::new("."), &INPUT_FMT, &DATATYPE);
//         let count = id.count_sites(seq);
//         assert_eq!(count, 5);
//     }
// }
