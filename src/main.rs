use std::time::Instant;

use segul::cli;
use segul::helper;

use segul::helper::logger;

#[cfg(not(tarpaulin_include))]
fn main() {
    let time = Instant::now();
    cli::parse_cli();

    let duration = time.elapsed();
    log::info!("{:18}: {}", "Log file", logger::LOG_FILE);
    println!();
    if duration.as_secs() < 60 {
        log::info!("{:18}: {:?}", "Execution time", duration);
    } else {
        let time = helper::utils::parse_duration(duration.as_secs());
        log::info!("{:18}: {}", "Execution time (HH:MM:SS)", time);
    }
}
