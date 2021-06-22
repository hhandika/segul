mod alignment;
mod cli;
mod common;
mod converter;
mod fasta;
mod finder;
mod msa;
mod nexus;
mod phylip;
mod stats;
mod utils;
mod writer;

use std::io::{BufWriter, Write};
use std::time::Instant;

use clap::crate_version;

fn main() {
    let version = crate_version!();
    let time = Instant::now();
    cli::parse_cli(&version);
    let duration = time.elapsed();
    let io = std::io::stdout();
    let mut writer = BufWriter::new(io);
    writeln!(writer).unwrap();
    writeln!(writer, "Execution time {:?}", duration).unwrap();
}
