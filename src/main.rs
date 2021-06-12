mod cli;
mod common;
mod converter;
mod fasta;
mod nexus;
mod phylip;

use std::time::Instant;

use clap::crate_version;

fn main() {
    let version = crate_version!();
    let time = Instant::now();
    cli::parse_cli(&version);
    let duration = time.elapsed();
    println!("Execution time {:?}", duration);
}
