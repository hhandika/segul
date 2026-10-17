mod utils;

use std::env;

#[test]
fn test_sequence_transpose() {
    initiate_cmd!(
        cmd,
        "sequence",
        "transpose",
        "crates/segul-core/tests/files/alignments/",
        tmp_dir
    );
    cmd.assert().success();
    let res_path = tmp_dir.path().join("Sequence-Transpose");
    for sample in ["ABCD", "ABCE", "ABCF"] {
        let output = res_path.join(sample).join(format!("{}.fas", sample));
        assert!(output.is_file(), "Missing {}", output.display());
    }
    tmp_dir.close().unwrap();
}
