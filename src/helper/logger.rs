//! Functions to setup the logger
use std::io::Result;
use std::path::Path;

use log::LevelFilter;
use log4rs::append::console::ConsoleAppender;
use log4rs::append::file::FileAppender;
use log4rs::config::{Appender, Config, Root};
use log4rs::encode::pattern::PatternEncoder;

/// Log file name
pub const LOG_FILE: &str = "segul.log";

/// Setup the logger
///
/// In the main function, call this function to setup the logger.
/// It initializes the logger with a console and file appender.
///
/// # Example
///
/// ```
/// use std::path::Path;
///
/// use segul::helper::logger;
///
/// fn main() {
///    let log_path = Path::new("segul.log");
///    logger::setup_logger(log_path).expect("Failed setting up logger");
/// }
/// ```
pub fn setup_logger(file_path: &Path) -> Result<()> {
    create_dir(file_path)?;
    let target = file_path.with_extension("log");
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
        .expect("Failed building log configuration");

    log4rs::init_config(config).expect("Cannot initiate log configuration");

    Ok(())
}

/// Setup to only log to file
/// # Example
/// ```
/// use std::path::Path;
/// use segul::helper::logger;
///
/// fn main() {
///   let output_dir = Path::new("output");
///   logger::setup_file_logger(output_dir);
/// }
pub fn setup_file_logger(output_dir: &Path) -> Result<()> {
    create_dir(output_dir)?;
    let target = output_dir.join(LOG_FILE);
    let tofile = FileAppender::builder()
        .encoder(Box::new(PatternEncoder::new(
            "{d(%Y-%m-%d %H:%M:%S %Z)} - {l} - {m}\n",
        )))
        .build(target)
        .expect("Failed building log file");

    let config = Config::builder()
        .appender(Appender::builder().build("logfile", Box::new(tofile)))
        .build(Root::builder().appender("logfile").build(LevelFilter::Info))
        .expect("Failed building log configuration");

    log4rs::init_config(config).expect("Cannot initiate log configuration");

    Ok(())
}

fn create_dir(file_path: &Path) -> Result<()> {
    let dir = file_path.parent().expect("Failed getting parent directory");
    if !dir.exists() {
        std::fs::create_dir_all(dir)?;
    }
    Ok(())
}
