mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::Files;
use segul::helper::types::InputFmt;

#[test]
fn test_dna_translation() {
    initiate_cmd!(cmd, "sequence", "translate", "tests/files/concat/", tmp_dir);
    cmd.assert().success();
    test_results!(4, tmp_dir, "Sequence-Translate", Nexus);
}
