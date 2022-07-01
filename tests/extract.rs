mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::Files;
use segul::helper::types::InputFmt;

#[test]
fn test_extract() {
    initiate_cmd!(cmd, "extract", tmp_dir);
    cmd.arg("--id").arg("ABCD").assert().success();

    let pred = predicates::path::is_dir();
    let res_path = tmp_dir.path().join("SEGUL-Extract");

    let files = Files::new(&res_path, &InputFmt::Nexus).find();

    assert!(pred.eval(&res_path));
    assert_eq!(4, files.len());

    tmp_dir.close().unwrap();
}
