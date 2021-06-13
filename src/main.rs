mod alignment;
mod cli;
mod common;
mod fasta;
mod nexus;
mod phylip;
mod writer;

use std::time::Instant;

use clap::crate_version;

fn main() {
    let version = crate_version!();
    let time = Instant::now();
    cli::parse_cli(&version);
    let duration = time.elapsed();
    println!("Execution time {:?}", duration);
}
