mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::SeqFileFinder;
use segul::helper::types::InputFmt;

#[test]
fn test_rename() {
    initiate_cmd!(cmd, "sequence", "rename", "tests/files/alignments/", tmp_dir);
    cmd.arg("--remove=D").assert().success();
    test_results!(4, tmp_dir, "Sequence-Rename", Nexus);
}
