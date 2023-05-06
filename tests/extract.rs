mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::SeqFileFinder;
use segul::helper::types::InputFmt;

#[test]
fn test_extract() {
    initiate_cmd!(cmd, "sequence", "extract", "tests/files/concat/", tmp_dir);
    cmd.arg("--id").arg("ABCD").assert().success();
    test_results!(4, tmp_dir, "Sequence-Extract", Nexus);
}

#[test]
#[should_panic]
fn test_conflicting_extract_cmd() {
    initiate_cmd!(cmd, "sequence", "extract", "tests/files/concat/", tmp_dir);
    cmd.arg("--id")
        .arg("ABCD")
        .arg("--re=^AB")
        .assert()
        .success();
}

#[test]
#[should_panic]
fn test_no_extract_cmd() {
    initiate_cmd!(cmd, "sequence", "extract", "tests/files/concat/", tmp_dir);
    cmd.arg("--id")
        .arg("ABCD")
        .arg("--file")
        .arg("tests/files/concat/concat.nex")
        .assert()
        .success();
}

#[test]
fn test_extract_re() {
    initiate_cmd!(cmd, "sequence", "extract", "tests/files/concat/", tmp_dir);
    cmd.arg("--re=^AB").assert().success();
    test_results!(4, tmp_dir, "Sequence-Extract", Nexus);
}
