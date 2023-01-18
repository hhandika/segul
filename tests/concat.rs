mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::Files;
use segul::helper::types::InputFmt;

#[test]
fn test_concat() {
    initiate_cmd!(cmd, "align", "concat", "tests/files/concat/", tmp_dir);
    cmd.arg("--part-format").arg("raxml").assert().success();
    test_results!(1, tmp_dir, "SEGUL-Concat", Nexus);
}

#[test]
fn test_concat_nexus_part() {
    initiate_cmd!(cmd, "align", "concat", "tests/files/concat/", tmp_dir);
    cmd.arg("--part-format").arg("nexus").assert().success();

    tmp_dir.close().unwrap();
}
