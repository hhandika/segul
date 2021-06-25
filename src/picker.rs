use std::fs;
use std::io::Result;
use std::path::Path;
use std::sync::{Arc, Mutex};

use rayon::prelude::*;

use crate::common::{Header, SeqFormat};
use crate::fasta::Fasta;
use crate::finder::{Files, IDs};
use crate::nexus::Nexus;
use crate::phylip::Phylip;

// Get Alignment
pub fn get_min_taxa(dir: &str, input_format: &SeqFormat, percent: f32, output_dir: &Path) {
    let files = Files::new(dir, input_format).get_files();
    let ntax = IDs::new(&files, input_format).get_id_all().len();
    let file_counts = files.len();
    let count = Arc::new(Mutex::new(0));
    fs::create_dir_all(output_dir).expect("CANNOT CREATE A TARGET DIRECTORY");
    files.par_iter().for_each(|file| {
        let min_taxa = count_min_tax(ntax, percent);
        let header = get_header(input_format, file);
        if header.ntax >= min_taxa {
            copy_files(file, output_dir).expect("CANNOT COPY FILES");
            let mut count = count.lock().unwrap();
            *count += 1;
        }
    });

    println!("File origin: {}", file_counts);
    println!("File final: {}", *count.lock().unwrap());
}

fn get_header(input_format: &SeqFormat, file: &Path) -> Header {
    match input_format {
        SeqFormat::Fasta | SeqFormat::FastaInt => get_fas_header(file),
        SeqFormat::Nexus | SeqFormat::NexusInt => get_nex_header(file),
        SeqFormat::Phylip => get_phy_header(file, false),
        SeqFormat::PhylipInt => get_phy_header(file, true),
    }
}

fn copy_files(origin: &Path, output_dir: &Path) -> Result<()> {
    let fname = origin.file_name().unwrap();
    let destination = output_dir.join(fname);

    fs::copy(origin, destination)?;

    Ok(())
}

fn get_nex_header(file: &Path) -> Header {
    let mut nex = Nexus::new(file);
    nex.read().unwrap();
    nex.header
}

fn get_phy_header(file: &Path, interleave: bool) -> Header {
    let mut phy = Phylip::new(file, interleave);
    phy.read().unwrap();
    phy.header
}

fn get_fas_header(file: &Path) -> Header {
    let mut fas = Fasta::new(file);
    fas.read();
    fas.header
}

fn count_min_tax(ntax: usize, percent: f32) -> usize {
    (ntax as f32 * percent).floor() as usize
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn min_taxa_test() {
        let ntax = 10;
        let percent = 0.65;
        assert_eq!(6, count_min_tax(ntax, percent));
    }
}
