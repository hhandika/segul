use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

use lazy_static::lazy_static;
use regex::Regex;

use crate::helper::types::{Partition, PartitionFmt};

macro_rules! parse_partition {
    ($self: ident, $pos: ident, $gene_name: ident, $partitions: ident, $current_start_pos: ident, $current_end_pos: ident) => {
        if $pos.contains(r#"\3"#) {
            let partition = $self.parse_partition($gene_name, $pos.trim(), true);
            let start_pos = partition.start;
            let end_pos = partition.end;
            assert!(
                start_pos != end_pos,
                "Invalid partition format. \
                Start and end position are the same."
            );
            if $current_end_pos != partition.end {
                if !$self.is_uncheck {
                    $self.check_partition_format(start_pos, $current_end_pos);
                }
                $partitions.push(partition);
            }
            $current_start_pos = start_pos;
            $current_end_pos = end_pos;
        } else {
            let partition = $self.parse_partition($gene_name, $pos.trim(), false);
            let start_pos = partition.start;
            let end_pos = partition.end;
            if !$self.is_uncheck {
                $self.check_partition_format(start_pos, $current_end_pos);
            }
            $partitions.push(partition);
            $current_start_pos = start_pos;
            $current_end_pos = end_pos;
        }
    };
}

macro_rules! assert_partition {
    ($self: ident, $partitions: ident) => {
        assert!(
            !$partitions.is_empty(),
            "Failed parsing partition. \
        No partition found."
        );
        if !$self.is_uncheck {
            assert_eq!(
                $partitions[0].start, 1,
                "Invalid partition input. \
                First partition start position is {} not 1.",
                $partitions[0].start
            );
        }
    };
}

pub struct PartitionParser<'a> {
    path: &'a Path,
    partition_fmt: &'a PartitionFmt,
    is_uncheck: bool,
}

impl<'a> PartitionParser<'a> {
    pub fn new(path: &'a Path, partition_fmt: &'a PartitionFmt, is_uncheck: bool) -> Self {
        Self {
            path,
            partition_fmt,
            is_uncheck,
        }
    }

    pub fn parse(&self) -> Vec<Partition> {
        let file = File::open(self.path).expect("Unable to open file");
        let mut reader = BufReader::new(file);
        match self.partition_fmt {
            PartitionFmt::Nexus
            | PartitionFmt::NexusCodon
            | PartitionFmt::Charset
            | PartitionFmt::CharsetCodon => self.parse_nexus(&mut reader),
            PartitionFmt::Raxml | PartitionFmt::RaxmlCodon => self.parse_raxml(&mut reader),
        }
    }

    fn check_partition_format(&self, start_pos: usize, curr_end_pos: usize) {
        assert_eq!(
            start_pos,
            curr_end_pos + 1,
            "Invalid partition format. \
        Start position ({}) is not the next position \
        after the previous end position ({}).",
            start_pos,
            curr_end_pos
        );
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
        let mut current_start_pos = 1;
        let mut current_end_pos = 0;
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
            parse_partition!(
                self,
                pos,
                gene_name,
                partitions,
                current_start_pos,
                current_end_pos
            );
        });
        assert_partition!(self, partitions);
        partitions
    }

    fn parse_nexus<R: BufRead>(&self, reader: &mut R) -> Vec<Partition> {
        let mut partitions = Vec::new();
        let mut current_start_pos = 1;
        let mut current_end_pos = 0;
        reader.lines().filter_map(|ok| ok.ok()).for_each(|line| {
            let nex_line = line.trim();
            if nex_line.to_lowercase().starts_with("charset") {
                let parts = line.split('=').collect::<Vec<&str>>();
                // We will split the gene part of the partition.
                // and exlude the first element, which is the "charset".
                let gene_parts = parts[0].split_whitespace().collect::<Vec<&str>>();
                let gene_name = gene_parts[1];
                let pos = parts[1].trim().replace(';', "");
                parse_partition!(
                    self,
                    pos,
                    gene_name,
                    partitions,
                    current_start_pos,
                    current_end_pos
                );
            }
        });
        assert_partition!(self, partitions);
        partitions
    }

    fn parse_partition<S: AsRef<str>>(&self, gene_name: S, pos: &str, is_codon: bool) -> Partition {
        let mut partition = Partition::new();
        partition.gene = gene_name.as_ref().trim().to_string();
        partition.gene.retain(|c| !r#"()/\,"';:?!"#.contains(c));
        if partition.gene.contains('.') {
            partition.gene = partition.gene.replace('.', "_");
        }
        assert!(
            !partition.gene.contains(' '),
            "Failed parsing partition file. Gene name cannot contain spaces"
        );
        let (start, end) = if is_codon {
            let codon = pos.replace(r#"\3"#, "");
            let subset = capture_subsets(&partition.gene);
            partition.gene = partition.gene.replace(&subset, "");
            self.parse_pos(&codon)
        } else {
            self.parse_pos(pos)
        };
        partition.start = start;
        partition.end = end;
        partition
    }

    fn parse_pos(&self, pos: &str) -> (usize, usize) {
        let parts = pos.split('-').collect::<Vec<&str>>();
        (
            parts[0]
                .trim()
                .parse::<usize>()
                .expect("Failed parsing gene start location"),
            parts[1]
                .trim()
                .parse::<usize>()
                .expect("Failed parsing gene end location"),
        )
    }
}

fn capture_subsets(text: &str) -> String {
    lazy_static! { // Match this formats: subset1, subset_1, or 1stpos
        static ref RE: Regex = Regex::new(r"(?i)(_subset_?\d|_\d\D{2}pos)").expect("Failed capturing partition subset");
    }
    match RE.captures(text) {
        Some(subset) => subset[0].to_string(),
        None => String::new(),
    }
}

#[cfg(test)]
mod test {
    use super::*;

    macro_rules! test_partition_parser {
        ($input:expr, $format:ident) => {
            let path = Path::new($input);
            let parser = PartitionParser::new(path, &PartitionFmt::$format, true);
            let partitions = parser.parse();
            assert_eq!(partitions.len(), 3);
            assert_eq!(partitions[0].gene, "locus1");
            assert_eq!(partitions[0].start, 1);
            assert_eq!(partitions[0].end, 5);
            assert_eq!(partitions[1].gene, "locus2");
            assert_eq!(partitions[1].start, 6);
            assert_eq!(partitions[1].end, 10);
            assert_eq!(partitions[2].gene, "locus3");
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

    #[test]
    fn test_parse_partition_raxml_codon() {
        let path = Path::new("test_files/partition/partition_codon.txt");
        test_partition_parser!(path, Raxml);
    }

    #[test]
    #[should_panic]
    fn test_parse_partition_raxml_inv() {
        let path = Path::new("test_files/partition/partition_inv_pos.txt");
        test_partition_parser!(path, Raxml);
    }

    #[test]
    #[should_panic]
    fn test_parse_partition_raxml_inv_start() {
        let path = Path::new("test_files/partition/partition_inv_start.txt");
        test_partition_parser!(path, Raxml);
    }

    #[test]
    fn test_subset_regex_match() {
        let subset_1 = "locus1_subset_1";
        let subset2 = "locus1_subset2";
        let subset_3 = "locus1_subset3";
        let stpos = "locus1_1stpos";
        let ndpos = "locus1_2ndpos";
        let rdpos = "locus1_3rdpos";
        assert_eq!(capture_subsets(subset_1), "_subset_1");
        assert_eq!(capture_subsets(subset2), "_subset2");
        assert_eq!(capture_subsets(subset_3), "_subset3");
        assert_eq!(capture_subsets(stpos), "_1stpos");
        assert_eq!(capture_subsets(ndpos), "_2ndpos");
        assert_eq!(capture_subsets(rdpos), "_3rdpos");
    }
}
