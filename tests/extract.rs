mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::Files;
use segul::helper::types::InputFmt;

#[test]
fn test_extract() {
    initiate_cmd!(cmd, "extract", "tests/files/concat/", tmp_dir);
    cmd.arg("--id").arg("ABCD").assert().success();
    test_results!(4, tmp_dir, "SEGUL-Extract", Nexus);
}

#[test]
fn test_extract_re() {
    initiate_cmd!(cmd, "extract", "tests/files/concat/", tmp_dir);
    cmd.arg("--re=^AB").assert().success();
    test_results!(4, tmp_dir, "SEGUL-Extract", Nexus);
}
