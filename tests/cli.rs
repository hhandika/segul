mod utils;
use clap::crate_version;

#[test]
fn test_version() {
    let version = crate_version!();
    utils::segul(utils::create_tmp_dir().unwrap().path())
        .arg("-V")
        .assert()
        .success()
        .stdout(predicates::str::contains(version));
}
