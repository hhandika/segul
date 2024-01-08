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
/// let output_dir = Path::new("output");
/// logger::init_file_logger(output_dir);
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
/// let task_desc = "Counting reads";
/// let fcounts = 2;
/// let logger = ReadLogger::new(&input, &input_fmt, &task_desc, fcounts);
/// logger.log();
/// ```
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
/// let input = Some(PathBuf::from("input"));
/// let input_fmt = ContigFmt::Fasta;
/// let task_desc = "Counting contigs";
/// let fcounts = 2;
/// let logger = ContigLogger::new(&input, &input_fmt, &task_desc, fcounts);
/// logger.log();
/// ```
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
/// let input = Some(PathBuf::from("input"));
/// let input_fmt = InputFmt::Fasta;
/// let datatype = DataType::Dna;
/// let task_desc = "Counting aligned sequences";
/// let fcounts = 2;
/// let logger = AlignSeqLogger::new(&input, &input_fmt, &datatype, &task_desc, fcounts);
/// logger.log();
/// ```

pub struct AlignSeqLogger<'a> {
    pub input: &'a Option<PathBuf>,
    pub input_fmt: &'a InputFmt,
    pub datatype: &'a DataType,
    pub fcounts: usize,
}

impl<'a> AlignSeqLogger<'a> {
    pub fn new(
        input: &'a Option<PathBuf>,
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
/// logger::log_input_partition(input, input_counts);
/// ```
pub fn log_input_partition(input: &Path, input_counts: usize) {
    log::info!("{:18}: {}", "Input dir", input.display());
    log::info!("{:18}: {}", "File counts", utils::fmt_num(&input_counts));
    log::info!("{:18}: {}\n", "Task", "Converting partitions");
}
