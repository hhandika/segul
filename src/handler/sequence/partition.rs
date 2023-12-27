//! Convert partitioned sequence file to another format.
use std::path::Path;

use colored::Colorize;

use crate::handler::PartitionPrint;
use crate::helper::types::{DataType, PartitionFmt};
use crate::helper::utils;
use crate::parser::partition::PartitionParser;
use crate::writer::partition::PartWriter;

impl PartitionPrint for PartConverter<'_> {}

/// Convert partitioned sequence file to another format.
///
/// It accepts input in in-file NEXUS, NEXUS, and RaXML format.
pub struct PartConverter<'a> {
    /// Input path in in-file NEXUS, NEXUS, or RaXML format.
    input: &'a Path,
    input_partition_fmt: &'a PartitionFmt,
    output: &'a Path,

    output_partition_fmt: &'a PartitionFmt,
}

impl<'a> PartConverter<'a> {
    /// Create a new PartConverter instance.
    pub fn new(
        input: &'a Path,
        input_partition_fmt: &'a PartitionFmt,
        output: &'a Path,
        output_partition_fmt: &'a PartitionFmt,
    ) -> Self {
        Self {
            input,
            input_partition_fmt,
            output,
            output_partition_fmt,
        }
    }

    /// Convert partitioned sequence file to another format.
    /// It accepts input in in-file NEXUS, NEXUS, and RaXML format.
    /// The output format can be in in-file NEXUS (Charset), NEXUS, or RaXML format.
    /// The output file will be named as `output_partition_fmt`.
    /// For example, if the output format is `Nexus`, the output file will be `partition.nex`.
    /// If the output format is `Raxml`, the output file will be `partition.txt`.
    /// If the output format is `Charset`, the output file will be `partition.nex`.
    ///
    /// # Example
    /// ```
    /// use std::path::{Path, PathBuf};
    /// use segul::handler::sequence::partition::PartConverter;
    /// use segul::helper::partition::construct_partition_path;
    /// use segul::helper::types::{DataType, PartitionFmt};
    /// use tempdir::TempDir;
    ///
    /// let input = Path::new("tests/files/partition/partition.nex");
    /// let output_dir = TempDir::new("temp").unwrap();
    /// let output_file = output_dir.path().join("partition");
    /// let final_output = construct_partition_path(&output_file, &PartitionFmt::Raxml);
    /// let mut converter = PartConverter::new(
    ///    &input,
    ///    &PartitionFmt::Charset,
    ///    &final_output,
    ///    &PartitionFmt::Raxml,
    /// );
    /// converter.convert(&DataType::Dna, false);
    /// ```
    pub fn convert(&self, datatype: &DataType, is_uncheck: bool) {
        let partitions =
            PartitionParser::new(self.input, self.input_partition_fmt, is_uncheck).parse();
        self.print_partition_info(self.input, &partitions.len());
        let spin = utils::set_spinner();
        spin.set_message("Converting partitions...");
        let writer = PartWriter::new(
            self.output,
            &partitions,
            self.output_partition_fmt,
            datatype,
        );
        writer.write_partition();
        spin.finish_with_message("Finished converting partitions!\n");
        self.print_output_info()
    }

    fn print_output_info(&self) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Output path", self.output.display());
    }
}
