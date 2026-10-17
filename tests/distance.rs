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

#[test]
fn test_align_distance_lower_triangle() {
    initiate_cmd!(
        cmd,
        "align",
        "distance",
        "crates/segul-core/tests/files/alignments/",
        tmp_dir
    );
    cmd.arg("--triangle").arg("lower");
    cmd.assert().success();
    let matrix = tmp_dir
        .path()
        .join("Align-Distance")
        .join("per_locus")
        .join("gene_2.csv");
    let content = std::fs::read_to_string(matrix).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines[1], "ABCD,0.000000,");
    assert_eq!(lines[2], "ABCF,0.125000,0.000000");
    tmp_dir.close().unwrap();
}

#[test]
fn test_align_distance_blank_diagonal() {
    initiate_cmd!(
        cmd,
        "align",
        "distance",
        "crates/segul-core/tests/files/alignments/",
        tmp_dir
    );
    cmd.arg("--triangle")
        .arg("lower")
        .arg("--diagonal")
        .arg("blank");
    cmd.assert().success();
    let matrix = tmp_dir
        .path()
        .join("Align-Distance")
        .join("per_locus")
        .join("gene_2.csv");
    let content = std::fs::read_to_string(matrix).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines[1], "ABCD,,");
    assert_eq!(lines[2], "ABCF,0.125000,");
    tmp_dir.close().unwrap();
}

#[test]
fn test_align_distance_group_long() {
    initiate_cmd!(
        cmd,
        "align",
        "distance",
        "crates/segul-core/tests/files/alignments/",
        tmp_dir
    );
    cmd.arg("--group-regex").arg("^ABC").arg("--long");
    cmd.assert().success();
    let res_path = tmp_dir.path().join("Align-Distance");
    for locus in LOCI {
        let output = res_path
            .join("per_locus_group")
            .join(format!("{}.csv", locus));
        assert!(output.is_file(), "Missing {}", output.display());
    }
    for output in [
        "distance_group_mean.csv",
        "distance_group_pairs.csv",
        "distance_group_long.csv",
        "distance_pairwise_long.csv",
    ] {
        assert!(res_path.join(output).is_file(), "Missing {}", output);
    }
    let long = std::fs::read_to_string(res_path.join("distance_pairwise_long.csv")).unwrap();
    assert!(long.starts_with("locus,sample_1,sample_2,group_1,group_2,compared_sites,distance"));
    tmp_dir.close().unwrap();
}

#[test]
fn test_align_distance_group_auto_fails_without_pattern() {
    initiate_cmd!(
        cmd,
        "align",
        "distance",
        "crates/segul-core/tests/files/alignments/",
        tmp_dir
    );
    cmd.arg("--group-auto");
    cmd.assert().failure();
    tmp_dir.close().unwrap();
}

#[test]
fn test_align_distance_group_delimiter_requires_fields() {
    initiate_cmd!(
        cmd,
        "align",
        "distance",
        "crates/segul-core/tests/files/alignments/",
        tmp_dir
    );
    cmd.arg("--group-delimiter").arg("_");
    cmd.assert().failure();
    tmp_dir.close().unwrap();
}
