use std::io::Result;
use std::panic;
use std::time::Instant;

use clap::crate_version;
use log::LevelFilter;
use log4rs::append::console::ConsoleAppender;
use log4rs::append::file::FileAppender;
use log4rs::config::{Appender, Config, Root};
use log4rs::encode::pattern::PatternEncoder;

mod cli;
mod core;
mod helper;
mod parser;
mod writer;

fn main() {
    panic::set_hook(Box::new(move |panic_info| {
        log::error!("{}", panic_info);
    }));

    setup_logger().expect("Failed setting up a log file.");
    let version = crate_version!();
    let time = Instant::now();
    cli::parse_cli(&version);
    let duration = time.elapsed();
    println!();
    if duration.as_secs() < 60 {
        log::info!("Execution time: {:?}", duration);
    } else {
        helper::utils::print_formatted_duration(duration.as_secs());
    }
}

fn setup_logger() -> Result<()> {
    let log_dir = std::env::current_dir()?;
    let target = log_dir.join("segul.log");
    let tofile = FileAppender::builder()
        .encoder(Box::new(PatternEncoder::new(
            "{d(%Y-%m-%d %H:%M:%S %Z)} - {l} - {m}\n",
        )))
        .build(target)?;

    let stdout = ConsoleAppender::builder()
        .encoder(Box::new(PatternEncoder::new("{m}\n")))
        .build();

    let config = Config::builder()
        .appender(Appender::builder().build("stdout", Box::new(stdout)))
        .appender(Appender::builder().build("logfile", Box::new(tofile)))
        .build(
            Root::builder()
                .appender("stdout")
                .appender("logfile")
                .build(LevelFilter::Info),
        )
        .unwrap();

    log4rs::init_config(config).unwrap();

    Ok(())
}
