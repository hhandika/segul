mod utils;

use std::env;

use predicates::Predicate;

#[test]
fn test_remove_id() {
    initiate_cmd!(cmd, "remove", tmp_dir);
    cmd.arg("--id").arg("ABCD").assert().success();
    tmp_dir.close().unwrap();
}
