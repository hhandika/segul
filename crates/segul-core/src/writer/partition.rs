//! Partition writer module
use std::io::prelude::*;
use std::path::Path;

use anyhow::Result;

use crate::helper::types::{DataType, Partition, PartitionFmt};
use crate::writer::FileWriter;

impl FileWriter for PartWriter<'_> {}

macro_rules! write_partition {
    ($self:ident, $write_type:ident, $partition:ident, $is_codon: expr_2021) => {{
        let mut writer = $self
            .$write_type($self.fpath)
            .expect("Failed creating/appending a partition file");
        $self
            .$partition(&mut writer, $is_codon)
            .expect("Failed writing a partition file");
    }};
}

pub struct PartWriter<'a> {
    fpath: &'a Path,
    partition: &'a [Partition],
    part_fmt: &'a PartitionFmt,
    datatype: &'a DataType,
}

impl<'a> PartWriter<'a> {
    pub fn new(
        fpath: &'a Path,
        partition: &'a [Partition],
        part_fmt: &'a PartitionFmt,
        datatype: &'a DataType,
    ) -> Self {
        Self {
            fpath,
            partition,
            part_fmt,
            datatype,
        }
    }

    pub fn write_partition(&self) {
        match self.part_fmt {
            PartitionFmt::Charset => {
                write_partition!(self, append_output_file, write_part_charset, false);
            }
            PartitionFmt::CharsetCodon => {
                write_partition!(self, append_output_file, write_part_charset, true);
            }
            PartitionFmt::Nexus => {
                write_partition!(self, create_output_file, write_part_nexus, false);
            }
            PartitionFmt::NexusCodon => {
                write_partition!(self, create_output_file, write_part_nexus, true);
            }
            PartitionFmt::Raxml => {
                write_partition!(self, create_output_file, write_part_raxml, false);
            }
            PartitionFmt::RaxmlCodon => {
                write_partition!(self, create_output_file, write_part_raxml, true)
            }
        }
    }

    fn write_part_nexus<W: Write>(&self, writer: &mut W, is_codon: bool) -> Result<()> {
        writeln!(writer, "#nexus")?;
        self.write_part_charset(writer, is_codon)?;
        Ok(())
    }

    fn write_part_raxml<W: Write>(&self, writer: &mut W, is_codon: bool) -> Result<()> {
        let dtype = if DataType::Dna == *self.datatype {
            "DNA, "
        } else {
            ""
        };
        self.partition.iter().for_each(|part| {
            if is_codon {
                self.write_raxml_codon(writer, part, dtype).unwrap();
            } else {
                writeln!(
                    writer,
                    "{}{} = {}-{}",
                    dtype, part.gene, part.start, part.end
                )
                .expect("Failed writing a partition file");
            }
        });
        writer.flush()?;
        Ok(())
    }

    fn write_part_charset<W: Write>(&self, writer: &mut W, is_codon: bool) -> Result<()> {
        writeln!(writer, "begin sets;")?;
        self.partition.iter().for_each(|part| {
            if is_codon {
                self.write_nex_codon(writer, part).unwrap();
            } else {
                writeln!(
                    writer,
                    "charset {} = {}-{};",
                    self.get_gene_name(&part.gene),
                    part.start,
                    part.end
                )
                .unwrap();
            }
        });
        writeln!(writer, "end;")?;
        writer.flush()?;
        Ok(())
    }

    fn get_gene_name(&self, name: &str) -> String {
        if name.contains('-') {
            format!("'{}'", name)
        } else {
            name.to_string()
        }
    }

    fn write_raxml_codon<W: Write>(
        &self,
        writer: &mut W,
        part: &Partition,
        dtype: &str,
    ) -> Result<()> {
        writeln!(
            writer,
            "{}{}_Subset1 = {}-{}\\3",
            dtype, part.gene, part.start, part.end
        )?;
        writeln!(
            writer,
            "{}{}_Subset2 = {}-{}\\3",
            dtype,
            part.gene,
            part.start + 1,
            part.end
        )?;
        writeln!(
            writer,
            "{}{}_Subset3 = {}-{}\\3",
            dtype,
            part.gene,
            part.start + 2,
            part.end
        )?;

        Ok(())
    }

    fn write_nex_codon<W: Write>(&self, writer: &mut W, part: &Partition) -> Result<()> {
        writeln!(
            writer,
            "charset {}_Subset1 = {}-{}\\3;",
            part.gene, part.start, part.end
        )?;
        writeln!(
            writer,
            "charset {}_Subset2 = {}-{}\\3;",
            part.gene,
            part.start + 1,
            part.end
        )?;
        writeln!(
            writer,
            "charset {}_Subset3 = {}-{}\\3;",
            part.gene,
            part.start + 2,
            part.end
        )?;

        Ok(())
    }
}
