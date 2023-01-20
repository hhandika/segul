use std::{fs::File, io::BufReader, path::Path};

use flate2::read::MultiGzDecoder;

pub fn decode_gzip(path: &Path) -> BufReader<MultiGzDecoder<File>> {
    let file = File::open(path).expect("Failed to open file");
    let decoder = MultiGzDecoder::new(file);
    BufReader::new(decoder)
}
