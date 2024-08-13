mod utils;

use std::env;

use predicates::Predicate;

use segul::helper::finder::SeqFileFinder;
use segul::helper::types::InputFmt;

#[test]
fn test_convert() {
    initiate_cmd!(cmd, "align", "convert", "tests/files/alignments/", tmp_dir);
    cmd.arg("-F").arg("phylip").assert().success();
    test_results!(4, tmp_dir, "Align-Convert", Phylip);
}
