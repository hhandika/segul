mod utils;

// use std::env;

// use predicates::Predicate;
// use segul::helper::finder::Files;
// use segul::helper::types::RawReadFmt;

macro_rules! initiate_raw_cmd {
    ($cmd: ident, $tmp_dir: ident) => {
        let $tmp_dir = utils::create_tmp_dir().unwrap();
        // let dir = env::current_dir().unwrap().join("tests/files/raw");
        let path = std::path::PathBuf::from($tmp_dir.path());
        let mut $cmd = utils::segul(&path);
        $cmd.arg("raw").arg("summary");
        // .arg("-d")
        // .arg(dir)
        // .arg("-f")
        // .arg("fastq");
    };
}

#[test]
fn test_raw_cmd() {
    initiate_raw_cmd!(cmd, tmp_dir);
    cmd.arg("--help").assert().success();
    // let pred = predicates::path::is_dir();
    // let res_path = tmp_dir.path().join("Raw-Summary");
    // let files = Files::new(&res_path).find_raw_read(&RawReadFmt::Fastq);

    // assert!(pred.eval(&res_path));
    // assert_eq!(2, files.len());
}
