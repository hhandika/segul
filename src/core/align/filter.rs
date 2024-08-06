//! Filtering alignment based on selected criteria.
//!
//! The alignment filtering work on the alignment level.
//! Currently available options are:
//! 1. Minimum number of taxa. Filter alignment that has an equal or more number of taxa than the specified number.
//! 2. Minimum alignment length. Filter alignment that has an equal or more alignment length than the specified number.
//! 3. Minimum parsimony informative sites. Filter alignment that has an equal or more parsimony informative sites than the specified number.
//! 4. Percentage of parsimony informative sites. Filter alignment that has an equal or more percentage of parsimony informative sites than the specified number.
//! 4. Taxon all. Filter alignment that contains all specified taxa.
//!
use std::fs;
use std::io::Result;
use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;

use colored::Colorize;
use indexmap::{IndexMap, IndexSet};
use rayon::prelude::*;

use crate::core::align::concat::AlignmentConcatenation;
use crate::helper::concat::ConcatParams;
use crate::helper::sequence::SeqParser;
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, PartitionFmt, SeqMatrix};
use crate::helper::{files, utils};
use crate::parser::fasta;
use crate::parser::nexus::Nexus;
use crate::parser::phylip::Phylip;
use crate::stats::sequence;

pub enum FilteringParameters {
    /// Filtered by minimum number of taxa.
    MinTax(usize),
    /// Filtered by minimum alignment length.
    AlnLen(usize),
    /// Filtered by minimum parsimony informative sites.
    ParsInf(usize),
    /// Filtered by percentage of parsimony informative sites.
    PercInf(f64),
    /// Filtered by taxa proportion of missing data. The value is in percentage.
    /// Missing data is defined as a gap "-" or missing data "?"
    MissingData(f64),
    /// Filtered by taxa that contains all specified taxa.
    TaxonAll(Vec<String>),
}

pub struct AlignmentFiltering<'a> {
    files: &'a [PathBuf],
    input_fmt: &'a InputFmt,
    datatype: &'a DataType,
    output: &'a Path,
    params: &'a FilteringParameters,
    concat: Option<ConcatParams>,
}

impl<'a> AlignmentFiltering<'a> {
    /// Create a new alignment filtering instance.
    pub fn new(
        files: &'a [PathBuf],
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        output: &'a Path,
        params: &'a FilteringParameters,
    ) -> Self {
        Self {
            files,
            input_fmt,
            datatype,
            output,
            params,
            concat: None,
        }
    }

    /// Filter the alignment based on the selected criteria.
    pub fn filter(&mut self) {
        let mut ftr_aln: Vec<PathBuf> = if let FilteringParameters::PercInf(perc_inf) = self.params
        {
            self.par_ftr_perc_inf(perc_inf)
        } else {
            self.par_ftr_aln()
        };

        assert!(!ftr_aln.is_empty(), "No alignments left after filtering!");

        match &self.concat {
            Some(concat) => self.concat_results(
                &mut ftr_aln,
                &concat.part_fmt,
                &concat.output_fmt,
                &concat.prefix,
            ),
            None => {
                let spin = utils::set_spinner();
                fs::create_dir_all(self.output).expect("CANNOT CREATE A TARGET DIRECTORY");
                spin.set_message("Copying matching alignments...");
                self.par_copy_files(&ftr_aln);
                spin.finish_with_message("Finished copying files!\n");
                self.print_output(ftr_aln.len());
            }
        }
    }

    /// Set the concatenation parameters. The resulting alignment will be concatenated.
    pub fn set_concat(
        &mut self,
        output_fmt: &'a OutputFmt,
        part_fmt: &'a PartitionFmt,
        prefix: &'a Path,
    ) {
        self.concat = Some(ConcatParams::new(output_fmt, part_fmt, prefix));
    }

    fn par_ftr_perc_inf(&self, perc_inf: &f64) -> Vec<PathBuf> {
        let spin = utils::set_spinner();
        spin.set_message("Counting parsimony informative sites...");
        let (send, rx) = channel();
        self.files.par_iter().for_each_with(send, |s, file| {
            s.send({
                let pinf = self.get_pars_inf(file);
                (PathBuf::from(file), pinf)
            })
            .unwrap()
        });
        spin.set_message("Finding maximum parsimony informative sites...");
        let ftr_aln: Vec<(PathBuf, usize)> = rx.iter().collect();
        let max_pinf = ftr_aln
            .iter()
            .map(|(_, pinf)| pinf)
            .max()
            .expect("Pinf contain none values");
        spin.finish_with_message("Finished counting pars. inf. sites!\n");
        let min_pinf = self.count_min_pinf(max_pinf, perc_inf);
        log::info!("{:18}: {}", "Max pinf. sites", max_pinf);
        log::info!("{:18}: {}\n", "Min pinf. sites", min_pinf);
        ftr_aln
            .iter()
            .filter(|(_, pinf)| *pinf >= min_pinf)
            .map(|(aln, _)| PathBuf::from(aln))
            .collect()
    }

    fn par_ftr_aln(&self) -> Vec<PathBuf> {
        let spin = utils::set_spinner();
        spin.set_message("Filtering alignments...");
        let (send, rx) = channel();
        self.files
            .par_iter()
            .for_each_with(send, |s, file| match self.params {
                FilteringParameters::MinTax(min_taxa) => {
                    let header = self.get_header(file);
                    if header.ntax >= *min_taxa {
                        s.send(file.to_path_buf()).expect("FAILED GETTING FILES");
                    }
                }
                FilteringParameters::AlnLen(nchar) => {
                    let header = self.get_header(file);
                    if header.nchar >= *nchar {
                        s.send(file.to_path_buf()).expect("FAILED GETTING FILES");
                    }
                }
                FilteringParameters::ParsInf(pars_inf) => {
                    let pars = self.get_pars_inf(file);
                    if pars >= *pars_inf {
                        s.send(file.to_path_buf()).expect("FAILED GETTING FILES");
                    }
                }
                FilteringParameters::TaxonAll(taxon_id) => {
                    let ids = self.parse_id(file);
                    if taxon_id.iter().all(|id| ids.contains(id)) {
                        s.send(file.to_path_buf()).expect("FAILED GETTING FILES");
                    }
                }
                FilteringParameters::MissingData(perc) => {
                    let missing_data = self.calculate_prop_missing_data(file);
                    if missing_data <= *perc {
                        s.send(file.to_path_buf()).expect("FAILED GETTING FILES");
                    }
                }
                _ => (),
            });

        let ftr_aln = rx.iter().collect();
        spin.finish_with_message("Finished filtering alignments!\n");
        ftr_aln
    }

    fn par_copy_files(&self, match_path: &[PathBuf]) {
        match_path.par_iter().for_each(|path| {
            self.copy_files(path).expect("Failed copying files");
        });
    }

    fn count_min_pinf(&self, max_inf: &usize, perc_inf: &f64) -> usize {
        (*max_inf as f64 * perc_inf).floor() as usize
    }

    fn concat_results(
        &self,
        ftr_files: &mut [PathBuf],
        part_fmt: &PartitionFmt,
        output_fmt: &OutputFmt,
        prefix: &Path,
    ) {
        let output_dir = files::create_output_fname(self.output, prefix, output_fmt);
        let mut concat =
            AlignmentConcatenation::new(self.input_fmt, &output_dir, output_fmt, part_fmt, prefix);
        concat.concat(ftr_files, self.datatype);
    }

    fn copy_files(&self, origin: &Path) -> Result<()> {
        let fname = origin.file_name().unwrap();
        let destination = self.output.join(fname);

        fs::copy(origin, destination)?;

        Ok(())
    }

    fn print_output(&self, fcounts: usize) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&fcounts));
        log::info!("{:18}: {}", "Dir", self.output.display());
    }

    fn calculate_prop_missing_data(&self, file: &Path) -> f64 {
        let (matrix, header) = self.get_alignment(file);
        let total_chars = header.nchar * header.ntax;
        let missing_data = self.count_missing_data(&matrix);
        missing_data as f64 / total_chars as f64
    }

    fn count_missing_data(&self, matrix: &SeqMatrix) -> usize {
        matrix
            .iter()
            .map(|(_, seq)| seq.bytes().filter(|&c| c == b'-' || c == b'?').count())
            .sum()
    }

    fn get_pars_inf(&self, file: &Path) -> usize {
        let (matrix, _) = self.get_alignment(file);
        sequence::get_pars_inf(&matrix, self.datatype)
    }

    fn parse_id(&self, file: &Path) -> IndexSet<String> {
        match self.input_fmt {
            InputFmt::Fasta => fasta::parse_only_id(file),
            InputFmt::Nexus => Nexus::new(file, self.datatype).parse_only_id(),
            InputFmt::Phylip => Phylip::new(file, self.datatype).parse_only_id(),
            _ => unreachable!("Auto format is not supported. Please, specify input format"),
        }
    }

    fn get_header(&self, file: &Path) -> Header {
        let (_, header) = self.get_alignment(file);
        header
    }

    fn get_alignment(&self, file: &Path) -> (IndexMap<String, String>, Header) {
        let aln = SeqParser::new(file, self.datatype);
        aln.get_alignment(self.input_fmt)
    }
}

#[cfg(test)]
mod test {
    use assert_approx_eq::assert_approx_eq;

    use super::*;
    use crate::helper::finder::SeqFileFinder;

    const PATH: &str = "tests/files/pinf/";
    const INPUT_FMT: InputFmt = InputFmt::Fasta;

    #[test]
    fn test_min_pinf() {
        let path = Path::new(PATH);
        let files = SeqFileFinder::new(path).find(&INPUT_FMT);
        let ftr = AlignmentFiltering::new(
            &files,
            &INPUT_FMT,
            &DataType::Dna,
            Path::new("test"),
            &FilteringParameters::PercInf(0.9),
        );

        let pinf = 4;
        let percent = 0.9;
        let percent_2 = 0.5;
        let ftr_aln = ftr.par_ftr_perc_inf(&percent);
        let ftr_aln_2 = ftr.par_ftr_perc_inf(&percent_2);
        assert_eq!(3, ftr.count_min_pinf(&pinf, &percent));
        assert_eq!(1, ftr_aln.len());
        assert_eq!(4, ftr_aln_2.len());
    }

    #[test]
    fn test_all_id() {
        let ids = vec!["1", "2", "3", "4"];
        let id_2 = vec!["1", "2", "3"];
        let id_3 = vec!["1", "2", "3", "4"];
        assert!(!ids.iter().all(|id| id_2.contains(id)));
        assert!(ids.iter().all(|id| id_3.contains(id)));
    }

    #[test]
    fn test_missing_data() {
        let path = Path::new("tests/files/gappy/");
        let input_fmt = InputFmt::Nexus;
        let files = SeqFileFinder::new(path).find(&input_fmt);
        let ftr = AlignmentFiltering::new(
            &files,
            &input_fmt,
            &DataType::Dna,
            &Path::new("test"),
            &FilteringParameters::MissingData(0.5),
        );
        let file = path.join("gene_1.nex");
        let missing_data = ftr.calculate_prop_missing_data(&file);
        assert_approx_eq!(0.27, missing_data, 2f64);
    }
}
