use assert_cmd::Command;

fn segul() -> Command {
    let mut cmd = Command::cargo_bin("segul").unwrap();
    cmd.current_dir("tests/files");
    cmd
}

#[test]
fn test_bin() {
    segul().arg("-V").assert().success();
}

// #[test]
// fn test_convert() {
//     segul().arg("convert").arg("-h").assert().success();
// }
