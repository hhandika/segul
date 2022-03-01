use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

#[allow(unused_imports)]
use crate::helper::types::Partition;

#[allow(dead_code)]
pub struct PartitionParser<'a> {
    path: &'a Path,
}

#[allow(dead_code)]
impl<'a> PartitionParser<'a> {
    pub fn new(path: &'a Path) -> Self {
        Self { path }
    }

    fn parse_raxml(&self) -> Vec<Partition> {
        let file = File::open(self.path).expect("Unable to open file");
        let reader = BufReader::new(file);
        let mut partitions = Vec::new();
        reader.lines().filter_map(|ok| ok.ok()).for_each(|line| {
            if !line.contains(',') {
                panic!("Invalid partition file format.");
            }
            let parts = line.split_whitespace().collect::<Vec<&str>>();
            assert_eq!(parts.len(), 4);
            partitions.push(self.parse_partition(&parts));
        });
        partitions
    }

    fn parse_nexus(&self) -> Vec<Partition> {
        let file = File::open(self.path).expect("Unable to open file");
        let reader = BufReader::new(file);
        let mut partitions = Vec::new();
        reader.lines().filter_map(|ok| ok.ok()).for_each(|line| {
            let nex_line = line.trim().replace(";", "");
            if nex_line.to_lowercase().starts_with("charset") {
                let parts = nex_line.split_whitespace().collect::<Vec<&str>>();
                assert_eq!(parts.len(), 4);
                partitions.push(self.parse_partition(&parts));
            }
        });
        partitions
    }

    fn parse_partition(&self, parts: &[&str]) -> Partition {
        let mut partition = Partition::new();
        partition.gene = parts[1].to_string();
        let genes = parts[3].split('-').collect::<Vec<&str>>();
        partition.start = genes[0]
            .parse::<usize>()
            .expect("Failed parsing gene start location");
        partition.end = genes[1]
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
        let parser = PartitionParser::new(path);
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
        let parser = PartitionParser::new(path);
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
