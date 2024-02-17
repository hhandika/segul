//! An archive writer module.
use std::{
    io::BufReader,
    path::{Path, PathBuf},
};

/// Archive multiple files into a single zip file.
pub struct Archive<'a> {
    /// The output path for the archive.
    pub output_path: &'a Path,
    /// The input files to be archived.
    pub input_files: &'a [PathBuf],
}

impl<'a> Archive<'a> {
    /// Create a new Archive instance.
    pub fn new(output_path: &'a Path, input_files: &'a [PathBuf]) -> Self {
        Self {
            output_path,
            input_files,
        }
    }

    /// Archive the input files into a single zip file.
    pub fn zip(&self) -> Result<(), Box<dyn std::error::Error>> {
        let mut archive = zip::ZipWriter::new(std::fs::File::create(self.output_path)?);
        let options =
            zip::write::FileOptions::default().compression_method(zip::CompressionMethod::Deflated);

        for file in self.input_files {
            archive.start_file(
                file.file_name()
                    .expect("Failed getting file name")
                    .to_str()
                    .expect("Failed converting file name to string"),
                options,
            )?;
            let input_file = std::fs::File::open(file)?;
            let mut buff = BufReader::new(input_file);
            std::io::copy(&mut buff, &mut archive)?;
        }

        archive.finish()?;

        Ok(())
    }
}
