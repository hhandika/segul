use std::{
    io::BufReader,
    path::{Path, PathBuf},
};

pub struct Archive<'a> {
    pub output_path: &'a Path,
    pub input_directory: &'a str,
    pub input_files: &'a [PathBuf],
}

impl<'a> Archive<'a> {
    pub fn new(
        output_path: &'a Path,
        input_directory: &'a str,
        input_files: &'a [PathBuf],
    ) -> Self {
        Self {
            output_path,
            input_directory,
            input_files,
        }
    }

    pub fn zip(&self) -> Result<(), Box<dyn std::error::Error>> {
        let mut archive = zip::ZipWriter::new(std::fs::File::create(&self.output_path)?);
        let options =
            zip::write::FileOptions::default().compression_method(zip::CompressionMethod::Deflated);
        archive.add_directory(self.input_directory, options)?;

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
