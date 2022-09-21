use std::path::Path;

use colored::Colorize;

use crate::handler::PartitionPrint;
use crate::helper::types::{DataType, PartitionFmt};
use crate::helper::utils;
use crate::parser::partition::PartitionParser;
use crate::writer::partition::PartWriter;

impl PartitionPrint for PartConverter<'_> {}

pub struct PartConverter<'a> {
    input: &'a Path,
    partition_fmt: &'a PartitionFmt,
    output: &'a Path,
    out_part_fmt: &'a PartitionFmt,
}

impl<'a> PartConverter<'a> {
    pub fn new(
        input: &'a Path,
        partition_fmt: &'a PartitionFmt,
        output: &'a Path,
        out_part_fmt: &'a PartitionFmt,
    ) -> Self {
        Self {
            input,
            partition_fmt,
            output,
            out_part_fmt,
        }
    }

    pub fn convert(&self, datatype: &DataType, is_uncheck: bool) {
        let partitions = PartitionParser::new(self.input, self.partition_fmt, is_uncheck).parse();
        self.print_partition_info(self.input, &partitions.len());
        let spin = utils::set_spinner();
        spin.set_message("Converting partitions...");
        let writer = PartWriter::new(self.output, &partitions, self.out_part_fmt, datatype);
        writer.write_partition();
        spin.finish_with_message("Finished converting partitions!\n");
        self.print_output_info()
    }

    fn print_output_info(&self) {
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Output path", self.output.display());
    }
}
