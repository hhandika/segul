mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::Files;
use segul::helper::types::InputFmt;

#[test]
fn test_remove_id() {
    initiate_cmd!(cmd, "remove", "tests/files/concat/", tmp_dir);
    cmd.arg("--id").arg("ABCD").assert().success();

    let pred = predicates::path::is_dir();
    let res_path = tmp_dir.path().join("SEGUL-Remove");

    let files = Files::new(&res_path, &InputFmt::Nexus).find();

    assert!(pred.eval(&res_path));
    assert_eq!(3, files.len());

    tmp_dir.close().unwrap();
}

#[test]
fn test_remove_re() {
    initiate_cmd!(cmd, "remove", "tests/files/concat/", tmp_dir);
    cmd.arg("--re=E$").assert().success();

    let pred = predicates::path::is_dir();
    let res_path = tmp_dir.path().join("SEGUL-Remove");

    let files = Files::new(&res_path, &InputFmt::Nexus).find();

    assert!(pred.eval(&res_path));
    assert_eq!(4, files.len());

    tmp_dir.close().unwrap();
}
