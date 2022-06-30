mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::Files;
use segul::helper::types::InputFmt;

#[test]
fn test_concat() {
    initiate_cmd!(cmd, "convert", tmp_dir);
    cmd.arg("-F").arg("phylip").assert().success();

    let pred = predicates::path::is_dir();
    let res_path = tmp_dir.path().join("SEGUL-Convert");

    let files = Files::new(&res_path, &InputFmt::Phylip).find();

    assert!(pred.eval(&res_path));
    assert_eq!(4, files.len());

    tmp_dir.close().unwrap();
}
