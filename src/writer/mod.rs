pub mod partition;
pub mod sequences;
pub mod summary;

use std::fs::{self, File, OpenOptions};
use std::io::BufWriter;
use std::path::Path;

use anyhow::{Context, Result};

trait FileWriter {
    fn create_output_file(&self, path: &Path) -> Result<BufWriter<File>> {
        let dir_name = path.parent().expect("Failed creating parent directory");
        fs::create_dir_all(&dir_name).with_context(|| {
            format!("Failed creating an output directory for {}", path.display())
        })?;
        let file = OpenOptions::new().write(true).create_new(true).open(path);
        match file {
            Ok(writer) => Ok(BufWriter::new(writer)),
            Err(error) => panic!("Failed writing to {}: {}", path.display(), error),
        }
    }

    fn append_output_file(&self, path: &Path) -> Result<BufWriter<File>> {
        let file = OpenOptions::new().append(true).open(path);
        match file {
            Ok(writer) => Ok(BufWriter::new(writer)),
            Err(error) => panic!("Failed appending to {}: {}", path.display(), error),
        }
    }
}
