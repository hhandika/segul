use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

use crate::helper::types::{Partition, PartitionFmt};

#[allow(dead_code)]
pub struct PartitionParser<'a> {
    path: &'a Path,
    partition_fmt: &'a PartitionFmt,
}

#[allow(dead_code)]
impl<'a> PartitionParser<'a> {
    pub fn new(path: &'a Path, partition_fmt: &'a PartitionFmt) -> Self {
        Self {
            path,
            partition_fmt,
        }
    }

    pub fn parse(&self) -> Vec<Partition> {
        match self.partition_fmt {
            PartitionFmt::Nexus => self.parse_nexus(),
            PartitionFmt::Raxml => self.parse_raxml(),
            _ => panic!("Unsupported partition format."),
        }
    }

    fn parse_raxml(&self) -> Vec<Partition> {
        let file = File::open(self.path).expect("Unable to open file");
        let reader = BufReader::new(file);
        let mut partitions = Vec::new();
        reader.lines().filter_map(|ok| ok.ok()).for_each(|line| {
            if !line.contains(',') {
                panic!("Invalid partition file format.");
            }
            let parts = line.trim().split('=').collect::<Vec<&str>>();
            partitions.push(self.parse_partition(&parts[0].trim(), &parts[1].trim()));
        });

        partitions
    }

    fn parse_nexus(&self) -> Vec<Partition> {
        let file = File::open(self.path).expect("Unable to open file");
        let reader = BufReader::new(file);
        let mut partitions = Vec::new();
        reader.lines().filter_map(|ok| ok.ok()).for_each(|line| {
            let nex_line = line.trim();
            if nex_line.to_lowercase().starts_with("charset") {
                let parts = line.split('=').collect::<Vec<&str>>();
                partitions.push(
                    self.parse_partition(&parts[0].trim(), &parts[1].replace(";", "").trim()),
                );
            }
        });
        partitions
    }

    fn parse_partition(&self, part_gene: &str, part_pos: &str) -> Partition {
        let mut partition = Partition::new();
        let gene_line = part_gene.split_whitespace().collect::<Vec<&str>>();
        partition.gene = gene_line[1].to_string();
        let pos = part_pos.split('-').collect::<Vec<&str>>();
        partition.start = pos[0]
            .trim()
            .parse::<usize>()
            .expect("Failed parsing gene start location");
        partition.end = pos[1]
            .trim()
            .parse::<usize>()
            .expect("Failed parsing gene end location");
        partition
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_parse_partition_raxml() {
        let path = Path::new("test_files/partition/partition.txt");
        let parser = PartitionParser::new(path, &PartitionFmt::Raxml);
        let partitions = parser.parse_raxml();
        assert_eq!(partitions.len(), 3);
        assert_eq!(partitions[0].gene, "Subset1");
        assert_eq!(partitions[0].start, 1);
        assert_eq!(partitions[0].end, 5);
        assert_eq!(partitions[1].gene, "Subset2");
        assert_eq!(partitions[1].start, 6);
        assert_eq!(partitions[1].end, 10);
        assert_eq!(partitions[2].gene, "Subset3");
        assert_eq!(partitions[2].start, 11);
        assert_eq!(partitions[2].end, 15);
    }

    #[test]
    fn test_parse_partition_nexus() {
        let path = Path::new("test_files/partition/partition.nex");
        let parser = PartitionParser::new(path, &PartitionFmt::Nexus);
        let partitions = parser.parse_nexus();
        assert_eq!(partitions.len(), 3);
        assert_eq!(partitions[0].gene, "Subset1");
        assert_eq!(partitions[0].start, 1);
        assert_eq!(partitions[0].end, 5);
        assert_eq!(partitions[1].gene, "Subset2");
        assert_eq!(partitions[1].start, 6);
        assert_eq!(partitions[1].end, 10);
        assert_eq!(partitions[2].gene, "Subset3");
        assert_eq!(partitions[2].start, 11);
        assert_eq!(partitions[2].end, 15);
    }
}
