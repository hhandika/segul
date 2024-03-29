//! Functions to setup the logger
use std::fs::create_dir_all;
use std::io::Result;
use std::path::Path;

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
/// It initialize logger to write to file and terminal output.
///
/// # Example
///
/// ```
/// use std::path::Path;
///
/// use segul::helper::logger;
///
/// let log_path = Path::new("segul.log");
/// logger::init_logger(log_path).expect("Failed setting up logger");
/// ```
pub fn init_logger(file_path: &Path) -> Result<()> {
    if let Some(dir) = file_path.parent() {
        create_dir_all(dir)?;
    }
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
/// let output_dir = Path::new("output");
/// logger::init_file_logger(output_dir);
pub fn init_file_logger(file_path: &Path) -> Result<()> {
    if let Some(dir) = file_path.parent() {
        create_dir_all(dir)?;
    }
    let target = file_path.with_extension("log");
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

/// Regex function to strip the terminal color text from file logger
pub fn strip_color(text: &str) -> String {
    let re = regex::Regex::new(r"\x1B\[([0-9]{1,2}(;[0-9]{1,2})?)?[m|K]").expect("Failed regex");
    re.replace_all(text, "").to_string()
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

/// Log input information for sequence reads
///
///
/// # Example
/// ```
/// use std::path::{Path, PathBuf};
/// use segul::helper::logger::{self, ReadLogger};
/// use segul::helper::types::SeqReadFmt;
///
/// let log_path = Path::new("segul.log");
/// logger::init_logger(log_path).expect("Failed setting up logger");
/// let input = Some(PathBuf::from("input"));
/// let input_fmt = SeqReadFmt::Fastq;
/// let task = "Counting reads";
/// let fcounts = 2;
/// let logger = ReadLogger::new(input.as_deref(), &input_fmt, fcounts);
/// logger.log(task);
/// ```
pub struct ReadLogger<'a> {
    pub input: Option<&'a Path>,
    pub input_fmt: &'a SeqReadFmt,
    pub fcounts: usize,
}

impl<'a> ReadLogger<'a> {
    pub fn new(input: Option<&'a Path>, input_fmt: &'a SeqReadFmt, fcounts: usize) -> Self {
        Self {
            input,
            input_fmt,
            fcounts,
        }
    }

    pub fn log(&self, task: &str) {
        self.log_input_info();
        log::info!("{:18}: {}\n", "Input format:", self.input_fmt);
        log::info!("{:18}: {}\n", "Task", task);
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

/// Log input information for contigs
///
/// # Example
/// ```
/// use std::path::{Path, PathBuf};
/// use segul::helper::logger::{self, ContigLogger};
/// use segul::helper::types::ContigFmt;
///
/// let log_path = Path::new("segul.log");
/// logger::init_logger(log_path).expect("Failed setting up logger");
/// let input = Some(Path::new("input"));
/// let input_fmt = ContigFmt::Fasta;
/// let task = "Counting contigs";
/// let fcounts = 2;
/// let logger = ContigLogger::new(input, &input_fmt,  fcounts);
/// logger.log(task);
/// ```
pub struct ContigLogger<'a> {
    pub input: Option<&'a Path>,
    pub input_fmt: &'a ContigFmt,
    pub fcounts: usize,
}

impl<'a> ContigLogger<'a> {
    pub fn new(input: Option<&'a Path>, input_fmt: &'a ContigFmt, fcounts: usize) -> Self {
        Self {
            input,
            input_fmt,
            fcounts,
        }
    }

    pub fn log(&self, task: &str) {
        self.log_input_info();
        log::info!("{:18}: {}\n", "Input format:", self.input_fmt);
        log::info!("{:18}: {}\n", "Task", task);
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

/// Log input information for aligned sequences
///
/// # Example
/// ```
/// use std::path::{Path, PathBuf};
/// use segul::helper::logger::{self, AlignSeqLogger};
/// use segul::helper::types::{InputFmt, DataType};
///
/// let log_path = Path::new("segul.log");
/// logger::init_logger(log_path).expect("Failed setting up logger");
/// let input = Some(Path::new("input"));
/// let input_fmt = InputFmt::Fasta;
/// let datatype = DataType::Dna;
/// let task = "Counting aligned sequences";
/// let fcounts = 2;
/// let logger = AlignSeqLogger::new(input, &input_fmt, &datatype, fcounts);
/// logger.log(task);
/// ```
pub struct AlignSeqLogger<'a> {
    pub input: Option<&'a Path>,
    pub input_fmt: &'a InputFmt,
    pub datatype: &'a DataType,
    pub fcounts: usize,
}

impl<'a> AlignSeqLogger<'a> {
    pub fn new(
        input: Option<&'a Path>,
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
        fcounts: usize,
    ) -> Self {
        Self {
            input,
            input_fmt,
            datatype,
            fcounts,
        }
    }

    pub fn log(&self, task: &str) {
        self.log_input_info();
        self.log_seq_input_fmt();
        self.log_seq_datatype();
        log::info!("{:18}: {}\n", "Task", task);
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

/// Log input information for partition conversion
///
/// # Example
/// ```
/// use std::path::Path;
/// use segul::helper::logger;
/// use segul::helper::utils;
///
/// let log_path = Path::new("segul.log");
/// logger::init_logger(log_path).expect("Failed setting up logger");
/// let input = Path::new("input");
/// let input_counts = 2;
/// logger::log_input_partition(Some(input), input_counts);
/// ```
pub fn log_input_partition(input: Option<&Path>, input_counts: usize) {
    match input {
        Some(input) => log::info!("{:18}: {}", "Input dir", input.display()),
        None => log::info!("{:18}: {}", "Input path", "STDIN"),
    }
    log::info!("{:18}: {}", "File counts", utils::fmt_num(&input_counts));
    log::info!("{:18}: {}\n", "Task", "Converting partitions");
}
