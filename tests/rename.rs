mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::Files;
use segul::helper::types::InputFmt;

#[test]
fn test_extract() {
    initiate_cmd!(cmd, "rename", "tests/files/concat/", tmp_dir);
    cmd.arg("--remove=D").assert().success();
    test_results!(4, tmp_dir, "SEGUL-Rename", Nexus);
}
