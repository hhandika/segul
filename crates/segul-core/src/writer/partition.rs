//! Partition writer module
use std::io::prelude::*;
use std::path::Path;

use anyhow::Result;

use crate::helper::types::{DataType, Partition, PartitionFmt, SiteRange};
use crate::writer::FileWriter;

impl FileWriter for PartWriter<'_> {}

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

    /// Write partitions to disk, panicking on I/O failure.
    pub fn write_partition(&self) {
        self.try_write_partition()
            .expect("Failed writing a partition file");
    }

    /// Write partitions to disk, returning any I/O error.
    pub fn try_write_partition(&self) -> Result<()> {
        match self.part_fmt {
            PartitionFmt::Charset => {
                let mut writer = self.append_output_file(self.fpath)?;
                self.write_part_charset(&mut writer, false)
            }
            PartitionFmt::CharsetCodon => {
                let mut writer = self.append_output_file(self.fpath)?;
                self.write_part_charset(&mut writer, true)
            }
            PartitionFmt::Nexus => {
                let mut writer = self.create_output_file(self.fpath)?;
                self.write_part_nexus(&mut writer, false)
            }
            PartitionFmt::NexusCodon => {
                let mut writer = self.create_output_file(self.fpath)?;
                self.write_part_nexus(&mut writer, true)
            }
            PartitionFmt::Raxml => {
                let mut writer = self.create_output_file(self.fpath)?;
                self.write_part_raxml(&mut writer, false)
            }
            PartitionFmt::RaxmlCodon => {
                let mut writer = self.create_output_file(self.fpath)?;
                self.write_part_raxml(&mut writer, true)
            }
        }
    }

    fn write_part_nexus<W: Write>(&self, writer: &mut W, is_codon: bool) -> Result<()> {
        writeln!(writer, "#nexus")?;
        self.write_part_charset(writer, is_codon)?;
        Ok(())
    }

    fn write_part_raxml<W: Write>(&self, writer: &mut W, is_codon: bool) -> Result<()> {
        for part in self.partition {
            let prefix = self.raxml_prefix(part);
            let ranges = part.effective_ranges();
            if is_codon && !part.has_stride() {
                for offset in 0..3usize {
                    let shifted = shift_codon_ranges(&ranges, offset);
                    if !shifted.is_empty() {
                        writeln!(
                            writer,
                            "{}{}_Subset{} = {}",
                            prefix,
                            part.gene,
                            offset + 1,
                            format_ranges(&shifted, ", ")
                        )?;
                    }
                }
            } else {
                writeln!(
                    writer,
                    "{}{} = {}",
                    prefix,
                    part.gene,
                    format_ranges(&ranges, ", ")
                )?;
            }
        }
        writer.flush()?;
        Ok(())
    }

    fn write_part_charset<W: Write>(&self, writer: &mut W, is_codon: bool) -> Result<()> {
        writeln!(writer, "begin sets;")?;
        let mut model_entries: Vec<(String, Option<&str>)> = Vec::new();
        for part in self.partition {
            let ranges = part.effective_ranges();
            let model = part.model.as_deref();
            if is_codon && !part.has_stride() {
                for offset in 0..3usize {
                    let shifted = shift_codon_ranges(&ranges, offset);
                    if !shifted.is_empty() {
                        let subset_name = format!("{}_Subset{}", part.gene, offset + 1);
                        let quoted = self.get_gene_name(&subset_name);
                        writeln!(
                            writer,
                            "charset {} = {};",
                            quoted,
                            format_ranges(&shifted, " ")
                        )?;
                        model_entries.push((quoted, model));
                    }
                }
            } else {
                let quoted = self.get_gene_name(&part.gene);
                writeln!(
                    writer,
                    "charset {} = {};",
                    quoted,
                    format_ranges(&ranges, " ")
                )?;
                model_entries.push((quoted, model));
            }
        }
        if model_entries.iter().any(|(_, m)| m.is_some()) {
            self.write_charpartition(writer, &model_entries)?;
        }
        writeln!(writer, "end;")?;
        writer.flush()?;
        Ok(())
    }

    fn write_charpartition<W: Write>(
        &self,
        writer: &mut W,
        entries: &[(String, Option<&str>)],
    ) -> Result<()> {
        let fallback = match self.datatype {
            DataType::Dna => "DNA",
            DataType::Aa => "AA",
            DataType::Ignore => "AUTO",
        };
        writeln!(writer, "charpartition models =")?;
        for (idx, (gene, model)) in entries.iter().enumerate() {
            let label = format_model_label(model.unwrap_or(fallback));
            let sep = if idx + 1 == entries.len() { ";" } else { "," };
            writeln!(writer, "    {}: {}{}", label, gene, sep)?;
        }
        Ok(())
    }

    fn raxml_prefix(&self, part: &Partition) -> String {
        if let Some(model) = &part.model {
            format!("{model}, ")
        } else if DataType::Dna == *self.datatype {
            "DNA, ".to_string()
        } else {
            String::new()
        }
    }

    fn get_gene_name(&self, name: &str) -> String {
        if name.contains(['-', '+', '.', '\'', ' ']) {
            format!("'{}'", name)
        } else {
            name.to_string()
        }
    }
}

fn shift_codon_ranges(ranges: &[SiteRange], offset: usize) -> Vec<SiteRange> {
    ranges
        .iter()
        .filter_map(|r| {
            let start = r.start + offset;
            (start <= r.end).then_some(SiteRange::with_step(start, r.end, 3))
        })
        .collect()
}

fn format_ranges(ranges: &[SiteRange], sep: &str) -> String {
    ranges
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(sep)
}

fn format_model_label(label: &str) -> String {
    if label.contains([' ', ';', ':', ',']) {
        format!("'{label}'")
    } else {
        label.to_string()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::parser::partition::PartitionParser;
    use tempfile::TempDir;

    #[test]
    fn test_write_merged_raxml_and_nexus_roundtrip() {
        let dir = TempDir::with_prefix("segul_part_writer").unwrap();
        let nex_path = dir.path().join("out.nex");
        let rax_path = dir.path().join("out.txt");

        let parts = vec![
            Partition::from_ranges(
                "Subset1",
                vec![SiteRange::new(1, 5), SiteRange::new(11, 15)],
                Some("GTR+I+G".to_string()),
            ),
            Partition::from_ranges(
                "Subset2",
                vec![SiteRange::with_step(6, 10, 3)],
                Some("HKY+G".to_string()),
            ),
        ];

        PartWriter::new(&nex_path, &parts, &PartitionFmt::Nexus, &DataType::Dna)
            .try_write_partition()
            .unwrap();
        let nex_text = std::fs::read_to_string(&nex_path).unwrap();
        assert!(nex_text.contains("charset Subset1 = 1-5 11-15;"));
        assert!(nex_text.contains("charset Subset2 = 6-10\\3;"));
        assert!(nex_text.contains("charpartition models ="));
        assert!(nex_text.contains("GTR+I+G: Subset1,"));
        assert!(nex_text.contains("HKY+G: Subset2;"));

        let parsed_nex = PartitionParser::new(&nex_path, &PartitionFmt::Nexus, true)
            .try_parse()
            .unwrap();
        assert_eq!(parsed_nex, parts);

        PartWriter::new(&rax_path, &parsed_nex, &PartitionFmt::Raxml, &DataType::Dna)
            .try_write_partition()
            .unwrap();
        let rax_text = std::fs::read_to_string(&rax_path).unwrap();
        assert_eq!(
            rax_text.lines().collect::<Vec<_>>(),
            vec!["GTR+I+G, Subset1 = 1-5, 11-15", "HKY+G, Subset2 = 6-10\\3",]
        );
    }

    #[test]
    fn test_write_codon_expansion_on_merged_unstrided() {
        let dir = TempDir::with_prefix("segul_part_codon").unwrap();
        let rax_path = dir.path().join("codon.txt");
        let parts = vec![Partition::from_ranges(
            "geneA",
            vec![SiteRange::new(1, 6), SiteRange::new(10, 12)],
            None,
        )];
        PartWriter::new(&rax_path, &parts, &PartitionFmt::RaxmlCodon, &DataType::Dna)
            .try_write_partition()
            .unwrap();
        let rax_text = std::fs::read_to_string(&rax_path).unwrap();
        assert_eq!(
            rax_text.lines().collect::<Vec<_>>(),
            vec![
                "DNA, geneA_Subset1 = 1-6\\3, 10-12\\3",
                "DNA, geneA_Subset2 = 2-6\\3, 11-12\\3",
                "DNA, geneA_Subset3 = 3-6\\3, 12-12\\3",
            ]
        );
    }
}
