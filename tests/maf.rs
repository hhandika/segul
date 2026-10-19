mod utils;

use std::env;
use std::path::PathBuf;

const SIMPLE_MAF: &str = "crates/segul-core/tests/files/maf/simple.maf";
const SIMPLE_REGIONS: &str = "crates/segul-core/tests/files/bed/simple_regions.bed";

fn fixture(path: &str) -> PathBuf {
    env::current_dir().unwrap().join(path)
}

#[test]
fn test_maf_convert_with_bed() {
    let tmp_dir = utils::create_tmp_dir().unwrap();
    let mut cmd = utils::segul(tmp_dir.path());
    cmd.args(["maf", "convert", "-i"])
        .arg(fixture(SIMPLE_MAF))
        .arg("--bed")
        .arg(fixture(SIMPLE_REGIONS));
    cmd.assert().success();
    let output = tmp_dir.path().join("Maf-Convert");
    for name in [
        "locA-chr7-27578828-27578866.fas",
        "clipA-chr7-27578830-27578834.fas",
        "gapIns-chr7-27578848-27578852.fas",
        "locMerge-chr7-27699739-27707234.fas",
        "maf_convert.csv",
    ] {
        assert!(output.join(name).is_file(), "Missing {}", name);
    }
    tmp_dir.close().unwrap();
}

#[test]
fn test_maf_convert_by_block_from_dir() {
    let tmp_dir = utils::create_tmp_dir().unwrap();
    let mut cmd = utils::segul(tmp_dir.path());
    cmd.args(["maf", "convert", "-d"])
        .arg(fixture("crates/segul-core/tests/files/maf"))
        .args(["-F", "nexus", "--min-score", "-1e9"]);
    cmd.assert().success();
    let output = tmp_dir.path().join("Maf-Convert");
    assert!(output.join("chr7-27578828-27578866.nex").is_file());
    tmp_dir.close().unwrap();
}

#[test]
fn test_maf_convert_rejects_group_by_with_bed() {
    let tmp_dir = utils::create_tmp_dir().unwrap();
    let mut cmd = utils::segul(tmp_dir.path());
    cmd.args(["maf", "convert", "-i"])
        .arg(fixture(SIMPLE_MAF))
        .arg("--bed")
        .arg(fixture(SIMPLE_REGIONS))
        .args(["--group-by", "chrom"]);
    cmd.assert().failure();
    tmp_dir.close().unwrap();
}

#[test]
fn test_old_maf_command_is_removed() {
    let tmp_dir = utils::create_tmp_dir().unwrap();
    let mut cmd = utils::segul(tmp_dir.path());
    cmd.args(["maf", "maf", "-i"]).arg(fixture(SIMPLE_MAF));
    cmd.assert().failure();
    tmp_dir.close().unwrap();
}

#[test]
fn test_maf_summary() {
    let tmp_dir = utils::create_tmp_dir().unwrap();
    let mut cmd = utils::segul(tmp_dir.path());
    cmd.args(["maf", "summary", "-i"]).arg(fixture(SIMPLE_MAF));
    cmd.assert().success();
    let output = tmp_dir.path().join("Maf-Summary");
    assert!(output.join("maf_summary_per_file.csv").is_file());
    assert!(output.join("maf_summary_per_species.csv").is_file());
    tmp_dir.close().unwrap();
}
