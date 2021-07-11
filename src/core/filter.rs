// use core::sync::atomic::{AtomicUsize, Ordering};
use std::fs;
use std::io::{self, BufWriter, Result, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;

use rayon::prelude::*;

use crate::helper::alignment::Alignment;
use crate::helper::common::{Header, SeqFormat};
use crate::helper::finder::IDs;
use crate::helper::utils;

#[allow(dead_code)]
pub enum Params {
    MinTax(usize),
    ParInf(usize),
    Nchar(usize),
}

// TODO:
// 1. Add support to concat the result
// 2. Allow more parameters, such as min aln length
pub struct SeqFilter<'a> {
    files: &'a mut [PathBuf],
    input_format: &'a SeqFormat,
    output_dir: &'a Path,
    percent: f64,
    min_taxa: usize,
    ntax: usize,
}

impl<'a> SeqFilter<'a> {
    pub fn new(
        files: &'a mut [PathBuf],
        input_format: &'a SeqFormat,
        output_dir: &'a Path,
        percent: f64,
    ) -> Self {
        Self {
            files,
            input_format,
            output_dir,
            percent,
            min_taxa: 0,
            ntax: 0,
        }
    }

    pub fn get_min_taxa(&mut self) {
        self.ntax = IDs::new(self.files, self.input_format).get_id_all().len();
        self.min_taxa = self.count_min_tax();
        self.display_input().expect("CANNOT DISPLAY TO STDOUT");
        fs::create_dir_all(self.output_dir).expect("CANNOT CREATE A TARGET DIRECTORY");
        let match_aln = self.par_match_aln();
        self.par_copy_files(&match_aln);
        self.display_output(match_aln.len())
            .expect("CANNOT DISPLAY TO STDOUT");
    }

    fn par_match_aln(&self) -> Vec<PathBuf> {
        let (send, rx) = channel();
        self.files.par_iter().for_each_with(send, |s, file| {
            let header = self.get_header(file);
            if header.ntax >= self.min_taxa {
                s.send(file.to_path_buf()).expect("FAILED GETTING FILES");
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

    fn display_input(&self) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(
            writer,
            "File count\t: {}",
            utils::fmt_num(&self.files.len())
        )?;
        writeln!(writer, "Taxon count\t: {}", self.ntax)?;
        writeln!(writer, "Percent\t\t: {}%", self.percent * 100.0)?;
        writeln!(writer, "Min tax\t\t: {}", self.min_taxa)?;
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

    fn get_header(&self, file: &Path) -> Header {
        let mut aln = Alignment::new();
        aln.get_aln_any(file, self.input_format);
        aln.header
    }

    fn count_min_tax(&self) -> usize {
        (self.ntax as f64 * self.percent).floor() as usize
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn min_taxa_test() {
        let mut files = [PathBuf::from(".")];
        let mut pick = SeqFilter::new(&mut files, &SeqFormat::Nexus, Path::new("."), 0.65);
        pick.ntax = 10;
        assert_eq!(6, pick.count_min_tax());
    }
}
