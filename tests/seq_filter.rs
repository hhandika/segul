mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::SeqFileFinder;
use segul::helper::types::InputFmt;

#[test]
fn test_sequence_filter_gaps() {
    initiate_cmd!(cmd, "sequence", "filter", "tests/files/gappy/", tmp_dir);
    cmd.arg("--max-gap").arg(".5").assert().success();
    test_results!(4, tmp_dir, "Sequence-Filter", Nexus);
}

#[test]
fn test_sequence_filter_length() {
    initiate_cmd!(cmd, "sequence", "filter", "tests/files/gappy/", tmp_dir);
    cmd.arg("--min-length").arg("10").assert().success();
    test_results!(1, tmp_dir, "Sequence-Filter", Nexus);
}
