mod utils;

use std::env;
use std::path::PathBuf;

use predicates::Predicate;

use segul::parser::txt;

macro_rules! initiate_part_cmd {
    ($cmd: ident, $tmp_dir: ident, $part: expr) => {
        let $tmp_dir = utils::create_tmp_dir().unwrap();
        let dir = env::current_dir().unwrap().join("tests/files/partition");
        let input = dir.join("partition.nex");
        let path = PathBuf::from($tmp_dir.path());
        let tmp_input = PathBuf::from(path.join("partition.nex"));
        std::fs::copy(&input, &tmp_input).unwrap();
        let mut $cmd = utils::segul(&path);
        $cmd.arg("partition")
            .arg("-i")
            .arg(tmp_input)
            .arg("-P")
            .arg($part);
    };
}

#[test]
fn test_partition() {
    initiate_part_cmd!(cmd, tmp_dir, "raxml");
    cmd.assert().success();
    let pred = predicates::path::is_file();
    let res_path = tmp_dir.path().join("partition_partition.txt");
    let part = txt::parse_text_file(&res_path);

    assert!(pred.eval(&res_path));
    assert_eq!(3, part.len());
}

#[test]
fn test_partition_codon() {
    initiate_part_cmd!(cmd, tmp_dir, "raxml");
    cmd.arg("--codon").assert().success();
    let res_path = tmp_dir.path().join("partition_partition.txt");
    let part = txt::parse_text_file(&res_path);
    assert_eq!(9, part.len());
}
