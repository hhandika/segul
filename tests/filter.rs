mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::SeqFileFinder;
use segul::helper::types::InputFmt;

macro_rules! test_filter {
    ($test: ident, $arg: expr, $val: expr, $output: expr, $res: expr) => {
        #[test]
        fn $test() {
            initiate_cmd!(cmd, "align", "filter", "tests/files/long-aln/", tmp_dir);
            cmd.arg($arg).arg($val).assert().success();
            test_results!($res, tmp_dir, $output, Nexus);
        }
    };
}

test_filter! {test_filter_percent, "--percent", "0.5", "Align-Filter_50p",4}
test_filter! {test_min_pinf, "--min-pinf", "3","Align-Filter_3pinf", 1}
test_filter! {test_max_pinf, "--max-pinf", "3","Align-Filter_3pinf", 3}
test_filter! {test_percent_percent_inf, "--percent-inf", ".5","Align-Filter_50percent_pinf", 4}
test_filter! {test_min_len, "--min-len", "25", "Align-Filter_25bp", 4}
test_filter! {test_max_len, "--max-len", "30", "Align-Filter_30bp", 4}

#[test]
fn test_filter_missing_data() {
    initiate_cmd!(cmd, "align", "filter", "tests/files/gappy/", tmp_dir);
    cmd.arg("--missing-data").arg("0.25").assert().success();
    let output = tmp_dir.path().join("Align-Filter_25percent_missing");
    test_results!(3, tmp_dir, output, Nexus);
}

#[test]
#[should_panic]
fn test_filter_panic() {
    initiate_cmd!(cmd, "align", "filter", "tests/files/long-aln/", tmp_dir);
    let output = tmp_dir.path().join("concat_50p");
    cmd.arg("--min-len")
        .arg("27")
        .arg("-o")
        .arg(&output)
        .assert()
        .success();
    tmp_dir.close().unwrap();
}
