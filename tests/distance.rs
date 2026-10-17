mod utils;

use std::env;

const LOCI: [&str; 4] = ["gene_1", "gene_2", "gene_3", "gene_4"];

#[test]
fn test_align_distance() {
    initiate_cmd!(
        cmd,
        "align",
        "distance",
        "crates/segul-core/tests/files/alignments/",
        tmp_dir
    );
    cmd.assert().success();
    let res_path = tmp_dir.path().join("Align-Distance");
    for locus in LOCI {
        let output = res_path.join("per_locus").join(format!("{}.csv", locus));
        assert!(output.is_file(), "Missing {}", output.display());
    }
    assert!(res_path.join("distance_summary.csv").is_file());
    tmp_dir.close().unwrap();
}

#[test]
fn test_align_distance_k2p() {
    initiate_cmd!(
        cmd,
        "align",
        "distance",
        "crates/segul-core/tests/files/alignments/",
        tmp_dir
    );
    cmd.arg("--model").arg("k2p").arg("--prefix").arg("k2p");
    cmd.assert().success();
    let res_path = tmp_dir.path().join("Align-Distance");
    assert!(res_path.join("k2p_summary.csv").is_file());
    tmp_dir.close().unwrap();
}

#[test]
fn test_align_distance_invalid_model() {
    initiate_cmd!(
        cmd,
        "align",
        "distance",
        "crates/segul-core/tests/files/alignments/",
        tmp_dir
    );
    cmd.arg("--model").arg("poisson");
    cmd.assert().failure();
    tmp_dir.close().unwrap();
}
