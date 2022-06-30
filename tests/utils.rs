use std::io::Result;
use std::path::Path;

use assert_cmd::Command;
use tempdir::TempDir;

pub const DIR: &str = "temp";

#[macro_export]
macro_rules! initiate_cmd {
    ($cmd: ident, $sub_cmd: expr, $tmp_dir: ident) => {
        let $tmp_dir = utils::create_tmp_dir().unwrap();
        let dir = env::current_dir().unwrap().join("tests/files/concat/");
        let path = std::path::PathBuf::from($tmp_dir.path());
        let mut $cmd = utils::segul(&path);
        $cmd.arg($sub_cmd).arg("-d").arg(dir).arg("-f").arg("nexus");
    };
}

pub fn segul(dir: &Path) -> Command {
    let mut cmd = Command::cargo_bin("segul").unwrap();
    cmd.current_dir(dir);
    cmd
}

pub fn create_tmp_dir() -> Result<TempDir> {
    let tmp_dir = TempDir::new(DIR)?;
    Ok(tmp_dir)
}
