use std::fs;
use std::io::Result;
use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;

use ansi_term::Colour::Yellow;
use indexmap::IndexMap;
use rayon::prelude::*;

use crate::core::concat::ConcatHandler;
use crate::helper::sequence::Sequence;
use crate::helper::stats;
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, PartitionFmt};
use crate::helper::utils;

pub enum Params {
    MinTax(usize),
    AlnLen(usize),
    ParsInf(usize),
    PercInf(f64),
}

pub struct SeqFilter<'a> {
    files: &'a [PathBuf],
    input_fmt: &'a InputFmt,
    datatype: &'a DataType,
    output: &'a Path,
    params: &'a Params,
    concat: Option<(&'a OutputFmt, &'a PartitionFmt)>,
}

impl<'a> SeqFilter<'a> {
    pub fn new(
        files: &'a [PathBuf],
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        output: &'a Path,
        params: &'a Params,
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

    pub fn filter_aln(&mut self) {
        let mut ftr_aln: Vec<PathBuf> = if let Params::PercInf(perc_inf) = self.params {
            self.par_ftr_perc_inf(perc_inf)
        } else {
            self.par_ftr_aln()
        };

        assert!(!ftr_aln.is_empty(), "No sequences left after filtering");

        match self.concat {
            Some((output_fmt, part_fmt)) => self.concat_results(&mut ftr_aln, output_fmt, part_fmt),
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

    pub fn set_concat(
        &mut self,
        output: &'a Path,
        output_fmt: &'a OutputFmt,
        part_fmt: &'a PartitionFmt,
    ) {
        self.output = output;
        self.concat = Some((output_fmt, part_fmt))
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
                Params::MinTax(min_taxa) => {
                    let header = self.get_header(file);
                    if header.ntax >= *min_taxa {
                        s.send(file.to_path_buf()).expect("FAILED GETTING FILES");
                    }
                }
                Params::AlnLen(nchar) => {
                    let header = self.get_header(file);
                    if header.nchar >= *nchar {
                        s.send(file.to_path_buf()).expect("FAILED GETTING FILES");
                    }
                }
                Params::ParsInf(pars_inf) => {
                    let pars = self.get_pars_inf(file);
                    if pars >= *pars_inf {
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
        output_fmt: &OutputFmt,
        part_fmt: &PartitionFmt,
    ) {
        let mut concat = ConcatHandler::new(self.input_fmt, self.output, output_fmt, part_fmt);
        concat.concat_alignment(ftr_files, self.datatype);
    }

    fn copy_files(&self, origin: &Path) -> Result<()> {
        let fname = origin.file_name().unwrap();
        let destination = self.output.join(fname);

        fs::copy(origin, destination)?;

        Ok(())
    }

    fn print_output(&self, fcounts: usize) {
        log::info!("{}", Yellow.paint("Output"));
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&fcounts));
        log::info!("{:18}: {}", "Dir", self.output.display());
    }

    fn get_pars_inf(&self, file: &Path) -> usize {
        let (matrix, _) = self.get_alignment(file);
        stats::get_pars_inf(&matrix, self.datatype)
    }

    fn get_header(&self, file: &Path) -> Header {
        let (_, header) = self.get_alignment(file);
        header
    }

    fn get_alignment(&self, file: &Path) -> (IndexMap<String, String>, Header) {
        let aln = Sequence::new(file, self.datatype);
        aln.get_alignment(self.input_fmt)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::helper::finder::Files;

    const PATH: &str = "test_files/pinf/";
    const INPUT_FMT: InputFmt = InputFmt::Fasta;

    #[test]
    fn test_min_pinf() {
        let path = Path::new(PATH);
        let files = Files::new(path, &INPUT_FMT).get_files();
        let ftr = SeqFilter::new(
            &files,
            &INPUT_FMT,
            &DataType::Dna,
            Path::new("test"),
            &Params::PercInf(0.9),
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
}
