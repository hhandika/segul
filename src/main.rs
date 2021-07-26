use std::io::Result;
use std::io::{BufWriter, Write};
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

    setup_logging().expect("Failed setting up a log file.");
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

fn setup_logging() -> Result<()> {
    let log_dir = std::env::current_dir()?;
    let target = log_dir.join("myte.log");
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
