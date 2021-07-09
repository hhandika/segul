mod alignment;
mod cli;
mod common;
mod converter;
mod filter;
mod finder;
mod msa;
mod parser;
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
    if duration.as_secs() < 60 {
        writeln!(writer, "Execution time: {:?}", duration).unwrap();
    } else {
        utils::print_formatted_duration(&mut writer, duration.as_secs()).unwrap();
    }
}
