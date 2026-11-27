mod utils;

use std::env;
use std::path::PathBuf;

use predicates::Predicate;

use segul::parser::txt;

macro_rules! initiate_part_cmd {
    ($cmd: ident, $tmp_dir: ident, $input_name: expr_2021, $part: expr_2021) => {
        let $tmp_dir = utils::create_tmp_dir().unwrap();
        let dir = env::current_dir().unwrap().join("tests/files/partition");
        let input = dir.join($input_name);
        let path = PathBuf::from($tmp_dir.path());
        let tmp_input = PathBuf::from(path.join($input_name));
        std::fs::copy(&input, &tmp_input).unwrap();
        let mut $cmd = utils::segul(&path);
        $cmd.arg("partition")
            .arg("convert")
            .arg("-i")
            .arg(tmp_input)
            .arg("-P")
            .arg($part);
    };
}

#[test]
fn test_partition() {
    initiate_part_cmd!(cmd, tmp_dir, "partition.nex", "raxml");
    cmd.assert().success();
    let pred = predicates::path::is_file();
    let res_path = tmp_dir.path().join("partition_partition.txt");
    let part = txt::parse_text_file(&res_path);

    assert!(pred.eval(&res_path));
    assert_eq!(3, part.len());
}

#[test]
fn test_partition_codon() {
    initiate_part_cmd!(cmd, tmp_dir, "partition.nex", "raxml");
    cmd.arg("--codon").assert().success();
    let res_path = tmp_dir.path().join("partition_partition.txt");
    let part = txt::parse_text_file(&res_path);
    assert_eq!(9, part.len());
}

#[test]
fn test_partition_merged_raxml_to_nexus_and_output_dir() {
    initiate_part_cmd!(cmd, tmp_dir, "partition_merged.txt", "nexus");
    let out_dir = tmp_dir.path().join("out");
    cmd.arg("-o").arg(&out_dir).assert().success();
    let res_path = out_dir.join("partition_merged_partition.nex");
    let content = std::fs::read_to_string(&res_path).unwrap();
    assert!(content.contains("charset Subset1 = 1-5\\3 2-5\\3 11-15;"));
    assert!(content.contains("charset Subset2 = 3-5\\3 6-10;"));
    assert!(content.contains("GTR+I+G: Subset1,"));
    assert!(content.contains("HKY+G: Subset2;"));
}

#[test]
fn test_partition_merged_nexus_to_raxml() {
    initiate_part_cmd!(cmd, tmp_dir, "partition_merged.nex", "raxml");
    cmd.assert().success();
    let res_path = tmp_dir.path().join("partition_merged_partition.txt");
    let lines = txt::parse_text_file(&res_path);
    assert_eq!(
        lines,
        vec![
            "DNA, Subset1 = 1-5\\3, 2-5\\3, 11-15",
            "DNA, Subset2 = 3-5\\3, 6-10",
        ]
    );
}

#[test]
fn test_partition_iqtree_charpartition_to_raxml() {
    initiate_part_cmd!(cmd, tmp_dir, "partition_iqtree.nex", "raxml");
    cmd.assert().success();
    let res_path = tmp_dir.path().join("partition_iqtree_partition.txt");
    let lines = txt::parse_text_file(&res_path);
    assert_eq!(
        lines,
        vec![
            "GTR+F+G4, locus1_locus3 = 1-5, 11-15",
            "HKY+F+I, locus2 = 6-10",
        ]
    );
}
