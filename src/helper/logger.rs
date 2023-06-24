//! Functions to setup the logger
use std::io::Result;
use std::path::{Path, PathBuf};

use log::LevelFilter;
use log4rs::append::console::ConsoleAppender;
use log4rs::append::file::FileAppender;
use log4rs::config::{Appender, Config, Root};
use log4rs::encode::pattern::PatternEncoder;

use crate::helper::types::{ContigFmt, DataType, InputFmt, SeqReadFmt};
use crate::helper::utils;

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
///    logger::init_logger(log_path).expect("Failed setting up logger");
/// }
/// ```
pub fn init_logger(file_path: &Path) -> Result<()> {
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
///   logger::init_file_logger(output_dir);
/// }
pub fn init_file_logger(output_dir: &Path) -> Result<()> {
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

trait InputLogger {
    fn log_input_info(&self) {}
}

impl InputLogger for ReadLogger<'_> {
    fn log_input_info(&self) {
        if let Some(input) = self.input {
            log::info!("{:18}: {}", "Input dir", &input.display());
        } else {
            log::info!("{:18}: {}", "Input path", "STDIN");
        }
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&self.fcounts));
    }
}

pub struct ReadLogger<'a> {
    pub input: &'a Option<PathBuf>,
    pub input_fmt: &'a SeqReadFmt,
    pub task_desc: &'a str,
    pub fcounts: usize,
}

impl<'a> ReadLogger<'a> {
    pub fn new(
        input: &'a Option<PathBuf>,
        input_fmt: &'a SeqReadFmt,
        task_desc: &'a str,
        fcounts: usize,
    ) -> Self {
        Self {
            input,
            input_fmt,
            task_desc,
            fcounts,
        }
    }

    pub fn log(&self) {
        self.log_input_info();
        log::info!("{:18}: {}\n", "Input format:", self.input_fmt);
        log::info!("{:18}: {}\n", "Task", self.task_desc);
    }
}

impl InputLogger for ContigLogger<'_> {
    fn log_input_info(&self) {
        if let Some(input) = self.input {
            log::info!("{:18}: {}", "Input dir", &input.display());
        } else {
            log::info!("{:18}: {}", "Input path", "STDIN");
        }
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&self.fcounts));
    }
}

pub struct ContigLogger<'a> {
    pub input: &'a Option<PathBuf>,
    pub input_fmt: &'a ContigFmt,
    pub task_desc: &'a str,
    pub fcounts: usize,
}

impl<'a> ContigLogger<'a> {
    pub fn new(
        input: &'a Option<PathBuf>,
        input_fmt: &'a ContigFmt,
        task_desc: &'a str,
        fcounts: usize,
    ) -> Self {
        Self {
            input,
            input_fmt,
            task_desc,
            fcounts,
        }
    }

    pub fn log(&self) {
        self.log_input_info();
        log::info!("{:18}: {}\n", "Input format:", self.input_fmt);
        log::info!("{:18}: {}\n", "Task", self.task_desc);
    }
}

impl InputLogger for AlignSeqLogger<'_> {
    fn log_input_info(&self) {
        if let Some(input) = self.input {
            log::info!("{:18}: {}", "Input dir", &input.display());
        } else {
            log::info!("{:18}: {}", "Input path", "STDIN");
        }
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&self.fcounts));
    }
}

pub struct AlignSeqLogger<'a> {
    pub input: &'a Option<PathBuf>,
    pub input_fmt: &'a InputFmt,
    pub datatype: &'a DataType,
    pub task_desc: &'a str,
    pub fcounts: usize,
}

impl<'a> AlignSeqLogger<'a> {
    pub fn new(
        input: &'a Option<PathBuf>,
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        task_desc: &'a str,
        fcounts: usize,
    ) -> Self {
        Self {
            input,
            input_fmt,
            datatype,
            task_desc,
            fcounts,
        }
    }

    pub fn log(&self) {
        self.log_input_info();
        self.log_seq_input_fmt();
        self.log_seq_datatype();
        log::info!("{:18}: {}\n", "Task", self.task_desc);
    }

    fn log_seq_datatype(&self) {
        match self.datatype {
            DataType::Aa => log::info!("{:18}: {}", "Data type", "Amino Acid"),
            DataType::Dna => log::info!("{:18}: {}", "Data type", "DNA"),
            DataType::Ignore => log::info!("{:18}: {}", "Data type", "Ignore"),
        }
    }

    fn log_seq_input_fmt(&self) {
        match self.input_fmt {
            InputFmt::Auto => log::info!("{:18}: {}", "Input format", "Auto"),
            InputFmt::Fasta => log::info!("{:18}: {}", "Input format", "FASTA"),
            InputFmt::Nexus => log::info!("{:18}: {}", "Input format", "NEXUS"),
            InputFmt::Phylip => log::info!("{:18}: {}", "Input format", "PHYLIP"),
        }
    }
}

pub fn log_input_partition(input: &Path, input_counts: usize) {
    log::info!("{:18}: {}", "Input dir", input.display());
    log::info!("{:18}: {}", "File counts", utils::fmt_num(&input_counts));
    log::info!("{:18}: {}\n", "Task", "Converting partitions");
}
