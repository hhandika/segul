use std::fs;
use std::io::Result;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use rayon::prelude::*;

use crate::common::{Header, SeqFormat};
use crate::fasta::Fasta;
use crate::finder::IDs;
use crate::nexus::Nexus;
use crate::phylip::Phylip;

pub struct Picker<'a> {
    files: &'a mut [PathBuf],
    input_format: &'a SeqFormat,
    output_dir: &'a Path,
    percent: f64,
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
            input_format,
            output_dir,
            percent,
        }
    }

    pub fn get_min_taxa(&self) {
        let ntax = IDs::new(self.files, self.input_format).get_id_all().len();
        let file_counts = self.files.len();
        let count = Arc::new(Mutex::new(0));
        let min_taxa = self.count_min_tax(ntax);
        fs::create_dir_all(self.output_dir).expect("CANNOT CREATE A TARGET DIRECTORY");
        self.files.par_iter().for_each(|file| {
            let header = self.get_header(file);
            if header.ntax >= min_taxa {
                self.copy_files(file).expect("CANNOT COPY FILES");
                let mut count = count.lock().unwrap();
                *count += 1;
            }
        });

        println!("File origin: {}", file_counts);
        println!("File final: {}", *count.lock().unwrap());
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

    fn count_min_tax(&self, ntax: usize) -> usize {
        (ntax as f64 * self.percent).floor() as usize
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn min_taxa_test() {
        let ntax = 10;
        let mut files = [PathBuf::from(".")];
        let pick = Picker::new(&mut files, &SeqFormat::Nexus, Path::new("."), 0.65);
        assert_eq!(6, pick.count_min_tax(ntax));
    }
}
