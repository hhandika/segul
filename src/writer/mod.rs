//! Writer trait and its implementations.
//!
//! The writer trait is used to write the output of the program to a file.
//! The trait is implemented for the following types:
//! 1. `SeqWriter`: write sequence data to a file.
//! 2. `PartitionWriter`: write partition data to a file.
//! 3. `SummaryWriter`: write summary data to a file.
//! 4. `ContigSummaryWriter`: write contig summary data to a file.
//! 5. `ArchiveWriter`: write archive data to a file.
//! 6. `ReadSummaryWriter`: write read summary data to a file.
pub mod archive;
pub mod contigs;
pub mod partition;
pub mod read;
pub mod sequences;
pub mod summary;
pub mod text;

use std::fs::{self, File, OpenOptions};
use std::io::BufWriter;
use std::path::Path;

use anyhow::{Context, Result};
use zip::ZipWriter;

trait FileWriter {
    fn create_output_file(&self, path: &Path) -> Result<BufWriter<File>> {
        create_parent_directory(path)?;
        let file = OpenOptions::new().write(true).create_new(true).open(path);
        match file {
            Ok(writer) => Ok(BufWriter::new(writer)),
            Err(error) => panic!("Failed writing to {}: {}", path.display(), error),
        }
    }

    /// File trait to write zip files.
    fn create_zip_writer(&self, output: &Path) -> Result<ZipWriter<File>> {
        create_parent_directory(output)?;
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output)?;
        let zip: ZipWriter<File> = ZipWriter::new(file);

        Ok(zip)
    }

    fn append_output_file(&self, path: &Path) -> Result<BufWriter<File>> {
        create_parent_directory(path)?;
        let file = OpenOptions::new().append(true).create(true).open(path);
        match file {
            Ok(writer) => Ok(BufWriter::new(writer)),
            Err(error) => panic!("Failed appending to {}: {}", path.display(), error),
        }
    }
}

fn create_parent_directory(path: &Path) -> Result<()> {
    let dir_name = path.parent().expect("Failed creating parent directory");
    fs::create_dir_all(dir_name)
        .with_context(|| format!("Failed creating an output directory for {}", path.display()))?;
    Ok(())
}
