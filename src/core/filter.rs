use std::fs;
use std::io::{self, BufWriter, Result, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;

use rayon::prelude::*;

use crate::core::stats;
use crate::helper::alignment::Alignment;
use crate::helper::common::{Header, SeqFormat};
use crate::helper::utils;

pub enum Params {
    MinTax(usize),
    AlnLen(usize),
    ParsInf(usize),
}

// TODO:
// 1. Add support to concat the result
pub struct SeqFilter<'a> {
    files: &'a [PathBuf],
    input_format: &'a SeqFormat,
    output_dir: &'a Path,
    params: &'a Params,
}

impl<'a> SeqFilter<'a> {
    pub fn new(
        files: &'a [PathBuf],
        input_format: &'a SeqFormat,
        output_dir: &'a Path,
        params: &'a Params,
    ) -> Self {
        Self {
            files,
            input_format,
            output_dir,
            params,
        }
    }

    pub fn get_min_taxa(&mut self) {
        fs::create_dir_all(self.output_dir).expect("CANNOT CREATE A TARGET DIRECTORY");
        let match_aln = self.par_match_aln();
        self.par_copy_files(&match_aln);
        self.display_output(match_aln.len())
            .expect("CANNOT DISPLAY TO STDOUT");
    }

    fn par_match_aln(&self) -> Vec<PathBuf> {
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
            });

        rx.iter().collect()
    }

    fn par_copy_files(&self, match_path: &[PathBuf]) {
        match_path.par_iter().for_each(|path| {
            self.copy_files(path).expect("CANNOT COPY FILES");
        });
    }

    fn copy_files(&self, origin: &Path) -> Result<()> {
        let fname = origin.file_name().unwrap();
        let destination = self.output_dir.join(fname);

        fs::copy(origin, destination)?;

        Ok(())
    }

    fn display_output(&self, fcounts: usize) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "\n\x1b[0;33mOutput\x1b[0m")?;
        writeln!(writer, "File count\t: {}", utils::fmt_num(&fcounts))?;
        writeln!(writer, "Dir\t\t: {}", self.output_dir.display())?;

        Ok(())
    }

    fn get_pars_inf(&self, file: &Path) -> usize {
        let aln = self.get_alignment(file);
        stats::get_pars_inf(&aln.alignment)
    }

    fn get_header(&self, file: &Path) -> Header {
        let aln = self.get_alignment(file);
        aln.header
    }

    fn get_alignment(&self, file: &Path) -> Alignment {
        let mut aln = Alignment::new();
        aln.get_aln_any(file, self.input_format);
        aln
    }
}

// #[cfg(test)]
// mod test {
//     use super::*;

//     #[test]
//     fn min_taxa_test() {
//         let mut files = [PathBuf::from(".")];
//         let mut pick = SeqFilter::new(&mut files, &SeqFormat::Nexus, Path::new("."), 0.65);
//         pick.ntax = 10;
//         assert_eq!(6, pick.count_min_tax());
//     }
// }
