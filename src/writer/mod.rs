//! Writer trait and its implementations.
//!
//! The writer trait is used to write the output of the program to a file.
//! The trait is implemented for the following types:
//! 1. `SeqWriter`: write sequence data to a file.
//! 2. `PartitionWriter`: write partition data to a file.
//! 3. `SummaryWriter`: write summary data to a file.
pub mod partition;
pub mod read;
pub mod sequences;
pub mod summary;

use std::fs::{self, File, OpenOptions};
use std::io::BufWriter;
use std::path::Path;

use anyhow::{Context, Result};

trait FileWriter {
    fn create_output_file(&self, path: &Path) -> Result<BufWriter<File>> {
        let dir_name = path.parent().expect("Failed creating parent directory");
        fs::create_dir_all(dir_name).with_context(|| {
            format!("Failed creating an output directory for {}", path.display())
        })?;
        let file = OpenOptions::new().write(true).create_new(true).open(path);
        match file {
            Ok(writer) => Ok(BufWriter::new(writer)),
            Err(error) => panic!("Failed writing to {}: {}", path.display(), error),
        }
    }

    fn append_output_file(&self, path: &Path) -> Result<BufWriter<File>> {
        let file = OpenOptions::new().append(true).create(true).open(path);
        match file {
            Ok(writer) => Ok(BufWriter::new(writer)),
            Err(error) => panic!("Failed appending to {}: {}", path.display(), error),
        }
    }
}
