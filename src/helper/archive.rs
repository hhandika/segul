use std::{io::BufReader, path::Path};

pub struct Archive {
    pub output_path: String,
    pub input_directory: String,
    pub input_files: Vec<String>,
}

impl Archive {
    pub fn new(output_path: String, input_directory: String, input_files: Vec<String>) -> Self {
        Self {
            output_path,
            input_directory,
            input_files,
        }
    }

    pub fn zip(&self) {
        let mut archive = zip::ZipWriter::new(std::fs::File::create(&self.output_path)?);
        let options =
            zip::write::FileOptions::default().compression_method(zip::CompressionMethod::Deflated);
        archive
            .add_directory(self.input_directory, options)
            .expect("Failed adding directory");

        for input in &self.input_files {
            let file = Path::new(input);
            archive.start_file(
                file.file_name()
                    .expect("Failed getting file name")
                    .to_str()
                    .expect("Failed converting file name to string"),
                options,
            )?;
            let mut input_file = std::fs::File::open(file).expect("Failed opening file");
            let mut buff = BufReader::new(input_file);
            std::io::copy(&mut buff, &mut archive).expect("Failed copying file");
        }

        archive.finish().expect("Failed finishing archive");
    }
}
