use std::panic;
use std::time::Instant;

use clap::crate_version;

mod cli;
mod handler;
mod helper;
mod parser;
mod writer;

fn main() {
    // We ignore backtrace for now. It does
    // not seem useful for most cases.
    let time = Instant::now();
    panic::set_hook(Box::new(move |panic_info| {
        log::error!("{}", panic_info);
    }));

    let version = crate_version!();
    cli::parse_cli(version);
    let duration = time.elapsed();
    log::info!("{:18}: {}", "Log file", cli::LOG_FILE);
    println!();
    if duration.as_secs() < 60 {
        log::info!("{:18}: {:?}", "Execution time", duration);
    } else {
        let time = helper::utils::parse_duration(duration.as_secs());
        log::info!("{:18}: {}", "Execution time (HH:MM:SS)", time);
    }
}
