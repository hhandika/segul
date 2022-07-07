mod utils;

#[test]
fn test_version() {
    use clap::crate_version;
    let version = crate_version!();
    utils::segul(utils::create_tmp_dir().unwrap().path())
        .arg("-V")
        .assert()
        .success()
        .stdout(predicates::str::contains(version));
}
