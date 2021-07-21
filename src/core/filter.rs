use std::fs;
use std::io::{self, BufWriter, Result, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;

use indexmap::IndexMap;
use rayon::prelude::*;

use crate::core::msa::MSAlignment;
use crate::helper::sequence::Sequence;
use crate::helper::stats;
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, PartitionFmt};
use crate::helper::utils;

pub enum Params {
    MinTax(usize),
    AlnLen(usize),
    ParsInf(usize),
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
        let mut ftr_aln = self.par_ftr_aln();
        match self.concat {
            Some((output_fmt, part_fmt)) => self.concat_results(&mut ftr_aln, output_fmt, part_fmt),
            None => {
                fs::create_dir_all(self.output).expect("CANNOT CREATE A TARGET DIRECTORY");
                self.par_copy_files(&ftr_aln);
                self.print_output(ftr_aln.len())
                    .expect("CANNOT DISPLAY TO STDOUT");
            }
        }
    }

    pub fn set_concat(&mut self, output_fmt: &'a OutputFmt, part_fmt: &'a PartitionFmt) {
        self.concat = Some((output_fmt, part_fmt))
    }

    fn par_ftr_aln(&self) -> Vec<PathBuf> {
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

    fn concat_results(
        &self,
        ftr_files: &mut [PathBuf],
        output_fmt: &OutputFmt,
        part_fmt: &PartitionFmt,
    ) {
        let output = self.output.to_string_lossy();
        let concat = MSAlignment::new(self.input_fmt, &output, output_fmt, part_fmt);
        concat.concat_alignment(ftr_files, &self.datatype);
    }

    fn copy_files(&self, origin: &Path) -> Result<()> {
        let fname = origin.file_name().unwrap();
        let destination = self.output.join(fname);

        fs::copy(origin, destination)?;

        Ok(())
    }

    fn print_output(&self, fcounts: usize) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "\x1b[0;33mOutput\x1b[0m")?;
        writeln!(writer, "File count\t: {}", utils::fmt_num(&fcounts))?;
        writeln!(writer, "Dir\t\t: {}", self.output.display())?;

        Ok(())
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
        let aln = Sequence::new(file, &self.datatype);
        aln.get_alignment(self.input_fmt)
    }
}
