use assert_cmd::Command;

fn segul() -> Command {
    let mut cmd = Command::cargo_bin("segul").unwrap();
    let dir = std::env::current_dir().unwrap().join("tests").join("files");
    cmd.current_dir(dir);
    cmd
}

#[test]
fn test_bin() {
    segul().arg("-V").assert().success();
}
