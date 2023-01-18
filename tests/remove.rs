mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::Files;
use segul::helper::types::InputFmt;

#[test]
fn test_remove_id() {
    initiate_cmd!(cmd, "sequence", "remove", "tests/files/concat/", tmp_dir);
    cmd.arg("--id").arg("ABCD").assert().success();
    test_results!(3, tmp_dir, "SEGUL-Remove", Nexus);
}

#[test]
fn test_remove_re() {
    initiate_cmd!(cmd, "sequence", "remove", "tests/files/concat/", tmp_dir);
    cmd.arg("--re=E$").assert().success();
    test_results!(4, tmp_dir, "SEGUL-Remove", Nexus);
}
