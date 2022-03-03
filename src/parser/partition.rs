use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

use crate::helper::types::{Partition, PartitionFmt};

pub struct PartitionParser<'a> {
    path: &'a Path,
    partition_fmt: &'a PartitionFmt,
}

impl<'a> PartitionParser<'a> {
    pub fn new(path: &'a Path, partition_fmt: &'a PartitionFmt) -> Self {
        Self {
            path,
            partition_fmt,
        }
    }

    pub fn parse(&self) -> Vec<Partition> {
        let file = File::open(self.path).expect("Unable to open file");
        let mut reader = BufReader::new(file);
        match self.partition_fmt {
            PartitionFmt::Nexus => self.parse_nexus(&mut reader),
            PartitionFmt::Raxml => self.parse_raxml(&mut reader),
            _ => panic!("Unsupported partition format."),
        }
    }

    /// This function will parse two kinds of raxml partition file:
    /// ```Text
    /// DNA, locus_1 = 1-100
    /// DNA, locus_2 = 101-200
    /// ```
    /// or
    /// ```Text
    /// locus_1 = 1-100
    /// locus_2 = 101-200
    fn parse_raxml<R: BufRead>(&self, reader: &mut R) -> Vec<Partition> {
        let mut partitions = Vec::new();
        reader.lines().filter_map(|ok| ok.ok()).for_each(|line| {
            let parts = line.trim().split('=').collect::<Vec<&str>>();
            let mut gene_name = parts[0].to_string();
            let pos: &str = parts[1].trim();

            // Raxml have so many different datatype written in the partition.
            // Here we will just brute force to remove it.
            if gene_name.contains(',') {
                let raxml_type = parts[0]
                    .find(',')
                    .expect("Error in parsing raxml partition file");
                gene_name.replace_range(..=raxml_type, "");
            }
            partitions.push(self.parse_partition(&gene_name, pos));
        });

        partitions
    }

    fn parse_nexus<R: BufRead>(&self, reader: &mut R) -> Vec<Partition> {
        let mut partitions = Vec::new();
        reader.lines().filter_map(|ok| ok.ok()).for_each(|line| {
            let nex_line = line.trim();
            if nex_line.to_lowercase().starts_with("charset") {
                let parts = line.split('=').collect::<Vec<&str>>();
                // We will split the gene part of the partition.
                // and exlude the first element, which is the "charset".
                let gene_parts = parts[0].split_whitespace().collect::<Vec<&str>>();
                let gene_name = gene_parts[1];
                let pos = parts[1].trim().replace(";", "");
                partitions.push(self.parse_partition(gene_name, pos.trim()));
            }
        });
        partitions
    }

    fn parse_partition(&self, gene_name: &str, pos: &str) -> Partition {
        let mut partition = Partition::new();
        partition.gene = gene_name.trim().to_string();
        assert!(
            !partition.gene.contains(' '),
            "Failed parsing partition file. Gene name cannot contain spaces"
        );
        let part_pos = pos.split('-').collect::<Vec<&str>>();
        partition.start = part_pos[0]
            .trim()
            .parse::<usize>()
            .expect("Failed parsing gene start location");
        partition.end = part_pos[1]
            .trim()
            .parse::<usize>()
            .expect("Failed parsing gene end location");
        partition
    }
}

#[cfg(test)]
mod test {
    use super::*;

    macro_rules! test_partition_parser {
        ($input:expr, $format:ident) => {
            let path = Path::new($input);
            let parser = PartitionParser::new(path, &PartitionFmt::$format);
            let partitions = parser.parse();
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
        };
    }

    #[test]
    fn test_parse_partition_raxml() {
        test_partition_parser!("test_files/partition/partition.txt", Raxml);
    }

    #[test]
    fn test_parse_partition_nexus() {
        test_partition_parser!("test_files/partition/partition.nex", Nexus);
    }

    #[test]
    fn test_parse_partition_raxml_with_whitespaces() {
        let path = Path::new("test_files/partition/partition_whitespaces.txt");
        test_partition_parser!(path, Raxml);
    }

    #[test]
    fn test_partition_nexus_with_whitespaces() {
        let path = Path::new("test_files/partition/partition_whitespaces.nex");
        test_partition_parser!(path, Nexus);
    }

    #[test]
    fn test_partition_raxml_no_datatype() {
        let path = Path::new("test_files/partition/partition_no_datatype.txt");
        test_partition_parser!(path, Raxml);
    }

    #[test]
    #[should_panic]
    fn test_parse_partition_raxml_with_invalid_format() {
        let path = Path::new("test_files/partition/partition_invalid.txt");
        test_partition_parser!(path, Raxml);
    }
}
