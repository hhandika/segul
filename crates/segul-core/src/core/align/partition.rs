//! Convert partitioned sequence file to another format.
use std::path::Path;

use anyhow::Result;
use colored::Colorize;

use crate::core::PartitionPrint;
use crate::helper::types::{DataType, PartitionFmt};
use crate::helper::utils;
use crate::parser::partition::PartitionParser;
use crate::writer::partition::PartWriter;

impl PartitionPrint for PartitionConverter<'_> {}

/// Convert partitioned sequence file to another format.
///
/// It accepts input in in-file NEXUS, NEXUS, and RaXML format.
pub struct PartitionConverter<'a> {
    /// Input path in in-file NEXUS, NEXUS, or RaXML format.
    input: &'a Path,
    input_partition_fmt: &'a PartitionFmt,
    output: &'a Path,
    output_partition_fmt: &'a PartitionFmt,
}

impl<'a> PartitionConverter<'a> {
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

    /// Convert partitioned sequence file to another format, panicking on failure.
    ///
    /// # Example
    /// ```
    /// use std::path::Path;
    /// use segul_core::core::align::partition::PartitionConverter;
    /// use segul_core::helper::partition::construct_partition_path;
    /// use segul_core::helper::types::{DataType, PartitionFmt};
    /// use tempfile::TempDir;
    ///
    /// let input = Path::new("tests/files/partition/partition.nex");
    /// let output_dir = TempDir::with_prefix("temp").unwrap();
    /// let output_file = output_dir.path().join("partition");
    /// let final_output = construct_partition_path(&output_file, &PartitionFmt::Raxml);
    /// let handle = PartitionConverter::new(
    ///    input,
    ///    &PartitionFmt::Charset,
    ///    &final_output,
    ///    &PartitionFmt::Raxml,
    /// );
    /// handle.convert(&DataType::Dna, false);
    /// ```
    pub fn convert(&self, datatype: &DataType, is_uncheck: bool) {
        self.try_convert(datatype, is_uncheck)
            .unwrap_or_else(|err| panic!("{err}"));
    }

    /// Convert partitioned sequence file to another format, returning an error on failure.
    pub fn try_convert(&self, datatype: &DataType, is_uncheck: bool) -> Result<()> {
        let partitions =
            PartitionParser::new(self.input, self.input_partition_fmt, is_uncheck).try_parse()?;
        self.print_partition_info(self.input, &partitions.len());
        let spin = utils::set_spinner();
        spin.set_message("Converting partitions...");
        let writer = PartWriter::new(
            self.output,
            &partitions,
            self.output_partition_fmt,
            datatype,
        );
        writer.try_write_partition()?;
        spin.finish_with_message("Finished converting partitions!\n");
        self.print_output_info();
        Ok(())
    }

    fn print_output_info(&self) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Output path", self.output.display());
    }
}
