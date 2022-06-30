use assert_cmd::Command;

pub fn segul() -> Command {
    let mut cmd = Command::cargo_bin("segul").unwrap();
    cmd.current_dir("tests/files");
    cmd
}

#[test]
fn test_version() {
    use clap::crate_version;
    let version = crate_version!();
    segul()
        .arg("-V")
        .assert()
        .success()
        .stdout(predicates::str::contains(version));
}
