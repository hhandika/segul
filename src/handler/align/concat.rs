//! Multi-Sequence Alignment Module.
//! Parse user input and call the appropriate functions
//! to concatenate sequence alignments.

use std::path::{Path, PathBuf};

use colored::Colorize;

use crate::helper::concat::Concat;
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, PartitionFmt};
use crate::helper::utils;
use crate::writer::{partition::PartWriter, sequences::SeqWriter};

pub struct ConcatHandler<'a> {
    input_fmt: &'a InputFmt,
    output: &'a Path,
    output_fmt: &'a OutputFmt,
    part_fmt: &'a PartitionFmt,
}

impl<'a> ConcatHandler<'a> {
    pub fn new(
        input_fmt: &'a InputFmt,
        output: &'a Path,
        output_fmt: &'a OutputFmt,
        part_fmt: &'a PartitionFmt,
    ) -> Self {
        Self {
            input_fmt,
            output,
            output_fmt,
            part_fmt,
        }
    }

    pub fn concat_alignment(&mut self, files: &mut [PathBuf], datatype: &DataType) {
        let mut concat = Concat::new(files, self.input_fmt, datatype);
        let spin = utils::set_spinner();
        concat.concat_alignment(&spin);
        let mut seq_writer = SeqWriter::new(self.output, &concat.alignment, &concat.header);
        let part_fname = self.construct_part_fpath();
        let part_writer = PartWriter::new(&part_fname, &concat.partition, self.part_fmt, datatype);
        spin.set_message("Writing output files...");
        seq_writer
            .write_sequence(self.output_fmt)
            .expect("Failed writing the output file");
        part_writer.write_partition();
        spin.finish_with_message("Finished concatenating alignments!\n");
        self.print_output_info(concat.partition.len(), &concat.header, &part_fname);
    }

    fn construct_part_fpath(&mut self) -> PathBuf {
        match self.part_fmt {
            PartitionFmt::Charset | PartitionFmt::CharsetCodon => PathBuf::from(self.output),
            PartitionFmt::Nexus | PartitionFmt::NexusCodon => self.get_part_fname("nex"),
            PartitionFmt::Raxml | PartitionFmt::RaxmlCodon => self.get_part_fname("txt"),
        }
    }

    fn get_part_fname(&self, ext: &str) -> PathBuf {
        let fname = format!(
            "{}_partition.{}",
            self.output
                .file_stem()
                .expect("Failed getting file name for partition")
                .to_string_lossy(),
            ext
        );

        self.output
            .parent()
            .expect("Failed getting output parent directory")
            .join(Path::new(&fname))
    }

    fn print_output_info(&self, count: usize, header: &Header, part_file: &Path) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Taxa", utils::fmt_num(&header.ntax));
        log::info!("{:18}: {}", "Loci", utils::fmt_num(&count));
        log::info!(
            "{:18}: {}",
            "Alignment length",
            utils::fmt_num(&header.nchar)
        );
        log::info!("{:18}: {}", "Alignment file", self.output.display());
        log::info!("{:18}: {}", "Partition file", &part_file.display(),);
    }
}
