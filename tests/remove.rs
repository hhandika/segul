mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::SeqFileFinder;
use segul::helper::types::InputFmt;

#[test]
fn test_remove_id() {
    initiate_cmd!(cmd, "sequence", "remove", "tests/files/alignments/", tmp_dir);
    cmd.arg("--id=ABCD").assert().success();
    test_results!(3, tmp_dir, "Sequence-Remove", Nexus);
}

#[test]
fn test_remove_re() {
    initiate_cmd!(cmd, "sequence", "remove", "tests/files/alignments/", tmp_dir);
    cmd.arg("--re=E$").assert().success();
    test_results!(4, tmp_dir, "Sequence-Remove", Nexus);
}
