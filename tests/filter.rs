mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::Files;
use segul::helper::types::InputFmt;

#[test]
fn test_extract() {
    initiate_cmd!(cmd, "filter", tmp_dir);
    let output = tmp_dir.path().join("concat_50p");
    cmd.arg("--percent")
        .arg("0.5")
        .arg("-o")
        .arg(&output)
        .assert()
        .success();

    let pred = predicates::path::is_dir();
    let files = Files::new(&output, &InputFmt::Nexus).find();

    assert!(pred.eval(&output));
    assert_eq!(4, files.len());

    tmp_dir.close().unwrap();
}
