use std::fs;
use std::io::{self, BufWriter, Result, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use rayon::prelude::*;

use crate::common::{Header, SeqFormat};
use crate::fasta::Fasta;
use crate::finder::IDs;
use crate::nexus::Nexus;
use crate::phylip::Phylip;
use crate::utils;

pub struct Picker<'a> {
    files: &'a mut [PathBuf],
    file_counts: usize,
    input_format: &'a SeqFormat,
    output_dir: &'a Path,
    percent: f64,
    min_taxa: usize,
    ntax: usize,
}

impl<'a> Picker<'a> {
    pub fn new(
        files: &'a mut [PathBuf],
        input_format: &'a SeqFormat,
        output_dir: &'a Path,
        percent: f64,
    ) -> Self {
        Self {
            files,
            file_counts: 0,
            input_format,
            output_dir,
            percent,
            min_taxa: 0,
            ntax: 0,
        }
    }

    pub fn get_min_taxa(&mut self) {
        self.ntax = IDs::new(self.files, self.input_format).get_id_all().len();
        self.file_counts = self.files.len();
        let fcounts = Arc::new(Mutex::new(0));
        self.min_taxa = self.count_min_tax();
        self.display_input().expect("CANNOT DISPLAY TO STDOUT");
        fs::create_dir_all(self.output_dir).expect("CANNOT CREATE A TARGET DIRECTORY");
        self.files.par_iter().for_each(|file| {
            let header = self.get_header(file);
            if header.ntax >= self.min_taxa {
                self.copy_files(file).expect("CANNOT COPY FILES");
                let mut fcounts = fcounts.lock().unwrap();
                *fcounts += 1;
            }
        });

        self.display_output(*fcounts.lock().unwrap())
            .expect("CANNOT DISPLAY TO STDOUT");
    }

    fn display_input(&self) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(
            writer,
            "File count\t: {}",
            utils::fmt_thousand_sep(&self.file_counts)
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
        writeln!(
            writer,
            "File count\t: {}",
            utils::fmt_thousand_sep(&fcounts)
        )?;
        writeln!(writer, "Dir\t\t: {}", self.output_dir.display())?;

        Ok(())
    }

    fn get_header(&self, file: &Path) -> Header {
        match self.input_format {
            SeqFormat::Fasta | SeqFormat::FastaInt => self.get_fas_header(file),
            SeqFormat::Nexus | SeqFormat::NexusInt => self.get_nex_header(file),
            SeqFormat::Phylip => self.get_phy_header(file, false),
            SeqFormat::PhylipInt => self.get_phy_header(file, true),
        }
    }

    fn copy_files(&self, origin: &Path) -> Result<()> {
        let fname = origin.file_name().unwrap();
        let destination = self.output_dir.join(fname);

        fs::copy(origin, destination)?;

        Ok(())
    }

    fn get_nex_header(&self, file: &Path) -> Header {
        let mut nex = Nexus::new(file);
        nex.read().unwrap();
        nex.header
    }

    fn get_phy_header(&self, file: &Path, interleave: bool) -> Header {
        let mut phy = Phylip::new(file, interleave);
        phy.read().unwrap();
        phy.header
    }

    fn get_fas_header(&self, file: &Path) -> Header {
        let mut fas = Fasta::new(file);
        fas.read();
        fas.header
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
        let mut pick = Picker::new(&mut files, &SeqFormat::Nexus, Path::new("."), 0.65);
        pick.ntax = 10;
        assert_eq!(6, pick.count_min_tax());
    }
}
