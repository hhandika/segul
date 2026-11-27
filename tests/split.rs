mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::SeqFileFinder;
use segul::helper::sequence::SeqParser;
use segul::helper::types::{DataType, InputFmt};

pub const DIR: &str = "temp";

macro_rules! initiate_split_cmd {
    ($cmd: ident, $tmp_dir: ident, $part: expr_2021) => {
        let $tmp_dir = utils::create_tmp_dir().unwrap();
        let dir = env::current_dir().unwrap().join("tests/files/partition");
        let input = dir.join("concat_part.fas");
        let partition = dir.join($part);
        let path = std::path::PathBuf::from($tmp_dir.path());
        let mut $cmd = utils::segul(&path);
        $cmd.arg("align")
            .arg("split")
            .arg("-i")
            .arg(input)
            .arg("-I")
            .arg(partition);
    };
}

#[test]
fn test_splitting_aln() {
    initiate_split_cmd!(cmd, tmp_dir, "partition.txt");
    cmd.assert().success();
    test_results!(3, tmp_dir, "Align-Split", Nexus);
}

#[test]
fn test_splitting_merged_raxml_and_codon() {
    initiate_split_cmd!(cmd, tmp_dir, "partition_merged.txt");
    cmd.arg("-F").arg("fasta").assert().success();
    let out_dir = tmp_dir.path().join("Align-Split");
    let (sub1, _) = SeqParser::new(&out_dir.join("Subset1.fas"), &DataType::Dna)
        .get_alignment(&InputFmt::Fasta);
    let (sub2, _) = SeqParser::new(&out_dir.join("Subset2.fas"), &DataType::Dna)
        .get_alignment(&InputFmt::Fasta);
    // Subset1 = 1-5\3 (1,4), 2-5\3 (2,5), 11-15 -> 1,2,4,5 + 11..=15
    assert_eq!(sub1.get("ABCG").unwrap(), b"aaaaccccc");
    // Subset2 = 3-5\3 (3), 6-10 -> 3 + 6..=10
    assert_eq!(sub2.get("ABCG").unwrap(), b"aggggg");
}

#[test]
fn test_splitting_iqtree_charpartition() {
    initiate_split_cmd!(cmd, tmp_dir, "partition_iqtree.nex");
    cmd.arg("-F").arg("fasta").assert().success();
    let out_dir = tmp_dir.path().join("Align-Split");
    let (merged, _) = SeqParser::new(&out_dir.join("locus1_locus3.fas"), &DataType::Dna)
        .get_alignment(&InputFmt::Fasta);
    let (single, _) =
        SeqParser::new(&out_dir.join("locus2.fas"), &DataType::Dna).get_alignment(&InputFmt::Fasta);
    assert_eq!(merged.get("ABCG").unwrap(), b"aaaaaccccc");
    assert_eq!(single.get("ABCG").unwrap(), b"ggggg");
}

#[test]
fn test_splitting_codon_file() {
    initiate_split_cmd!(cmd, tmp_dir, "partition_codon.txt");
    cmd.arg("-F").arg("fasta").assert().success();
    test_results!(9, tmp_dir, "Align-Split", Fasta);
}
