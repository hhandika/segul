mod utils;

#[test]
fn test_view_help() {
    utils::segul(utils::create_tmp_dir().unwrap().path())
        .arg("view")
        .arg("align")
        .arg("--help")
        .assert()
        .success()
        .stdout(predicates::str::contains("--theme"));
}

#[test]
fn test_view_missing_file() {
    utils::segul(utils::create_tmp_dir().unwrap().path())
        .arg("view")
        .arg("align")
        .arg("-i")
        .arg("missing.fas")
        .assert()
        .failure()
        .stderr(predicates::str::contains("missing.fas"));
}
