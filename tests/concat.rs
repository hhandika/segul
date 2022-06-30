mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::Files;
use segul::helper::types::InputFmt;

#[test]
fn test_concat() {
    initiate_cmd!(cmd, "concat", tmp_dir);
    cmd.arg("--part").arg("raxml").assert().success();

    let pred = predicates::path::is_dir();
    let res_path = tmp_dir.path().join("SEGUL-Concat");

    let files = Files::new(&res_path, &InputFmt::Nexus).find();

    assert!(pred.eval(&res_path));
    assert_eq!(1, files.len());

    tmp_dir.close().unwrap();
}

#[test]
fn test_concat_nexus_part() {
    initiate_cmd!(cmd, "concat", tmp_dir);
    cmd.arg("--part").arg("nexus").assert().success();

    tmp_dir.close().unwrap();
}
