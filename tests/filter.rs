mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::Files;
use segul::helper::types::InputFmt;

macro_rules! test_filter {
    ($test: ident, $arg: expr, $val: expr, $res: expr) => {
        #[test]
        fn $test() {
            initiate_cmd!(cmd, "filter", "tests/files/filter", tmp_dir);
            let path = "concat_50p";
            let output = tmp_dir.path().join(path);
            cmd.arg($arg)
                .arg($val)
                .arg("-o")
                .arg(&output)
                .assert()
                .success();
            test_results!($res, tmp_dir, path, Nexus);
        }
    };
}

test_filter! {test_filter_percent, "--percent", "0.5", 4}
test_filter! {test_percent_pinf, "--pinf", "3", 1}
test_filter! {test_percent_percent_inf, "--percent-inf", ".5", 4}
test_filter! {test_percent_len, "--len", "25", 4}

#[test]
#[should_panic]
fn test_filter_panic() {
    initiate_cmd!(cmd, "filter", "tests/files/filter", tmp_dir);
    let output = tmp_dir.path().join("concat_50p");
    cmd.arg("--len")
        .arg("27")
        .arg("-o")
        .arg(&output)
        .assert()
        .success();
    tmp_dir.close().unwrap();
}
