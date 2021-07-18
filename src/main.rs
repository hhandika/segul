use std::io::{BufWriter, Write};
use std::time::Instant;

use clap::crate_version;

#[macro_use]
extern crate lazy_static;

mod cli;
mod core;
mod helper;
mod parser;
mod writer;

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
        helper::utils::print_formatted_duration(&mut writer, duration.as_secs()).unwrap();
    }
}
