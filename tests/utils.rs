use std::io::Result;
use std::path::Path;

use assert_cmd::{Command, cargo};
use tempdir::TempDir;

pub const DIR: &str = "temp";

#[macro_export]
macro_rules! initiate_cmd {
    ($cmd: ident, $sub_cmd1: expr_2021, $sub_cmd2: expr_2021, $dir: expr_2021, $tmp_dir: ident) => {
        let $tmp_dir = utils::create_tmp_dir().unwrap();
        let dir = env::current_dir().unwrap().join($dir);
        let path = std::path::PathBuf::from($tmp_dir.path());
        let mut $cmd = utils::segul(&path);
        $cmd.arg($sub_cmd1)
            .arg($sub_cmd2)
            .arg("-d")
            .arg(dir)
            .arg("-f")
            .arg("nexus");
    };
}

#[macro_export]
macro_rules! test_results {
    ($res: expr_2021, $tmp_dir: ident, $path: expr_2021, $fmt: ident) => {
        let pred = predicates::path::is_dir();
        let res_path = $tmp_dir.path().join($path);

        let files = SeqFileFinder::new(&res_path).find(&InputFmt::$fmt);

        assert!(pred.eval(&res_path));
        assert_eq!($res, files.len());

        $tmp_dir.close().unwrap();
    };
}

pub fn segul(dir: &Path) -> Command {
    let mut cmd = cargo::cargo_bin_cmd!("segul");
    cmd.current_dir(dir);
    cmd
}
pub fn create_tmp_dir() -> Result<TempDir> {
    let tmp_dir = TempDir::new(DIR)?;
    Ok(tmp_dir)
}
