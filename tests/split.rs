mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::Files;
use segul::helper::types::InputFmt;
// use std::path::Path;

// use assert_cmd::Command;
// use tempdir::TempDir;

pub const DIR: &str = "temp";

macro_rules! initiate_split_cmd {
    ($cmd: ident, $tmp_dir: ident, $part: expr) => {
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
