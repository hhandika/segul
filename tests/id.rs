mod utils;

use std::env;

use predicates::Predicate;

use segul::parser::txt;

#[test]
fn test_id_success() {
    initiate_cmd!(cmd, "id", "tests/files/concat/", tmp_dir);
    cmd.assert().success();
    let pred = predicates::path::is_file();
    let res_path = tmp_dir.path().join("id.txt");
    let ids = txt::parse_text_file(&res_path);

    assert!(pred.eval(&res_path));
    assert_eq!(3, ids.len());

    tmp_dir.close().unwrap();
}
