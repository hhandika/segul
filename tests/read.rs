mod utils;

use std::env;
use std::fs;
use std::path::{Path, PathBuf};

const QC_DIR: &str = "crates/segul-core/tests/files/raw_qc";
const OUTPUT: &str = "Read-Summary";

fn fixture(dir: &str) -> PathBuf {
    env::current_dir().unwrap().join(dir)
}

macro_rules! initiate_raw_cmd {
    ($cmd: ident, $tmp_dir: ident, $dir: expr_2021) => {
        let $tmp_dir = utils::create_tmp_dir().unwrap();
        let mut $cmd = utils::segul($tmp_dir.path());
        $cmd.arg("read")
            .arg("summary")
            .arg("-d")
            .arg(fixture($dir))
            .args(["-o", OUTPUT]);
    };
}

fn header(path: &Path) -> String {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .next()
        .unwrap()
        .to_string()
}

fn zip_entries(path: &Path) -> Vec<String> {
    let file = fs::File::open(path).unwrap();
    let archive = zip::ZipArchive::new(file).unwrap();
    archive
        .file_names()
        .map(|n| n.unwrap().to_string())
        .collect()
}

#[test]
fn test_read_cmd() {
    let tmp_dir = utils::create_tmp_dir().unwrap();
    let mut cmd = utils::segul(tmp_dir.path());
    cmd.args(["read", "summary", "--help"]).assert().success();
}

#[test]
fn test_read_summary_minimal() {
    initiate_raw_cmd!(cmd, tmp_dir, QC_DIR);
    cmd.args(["--mode", "minimal"]).assert().success();
    let summary = tmp_dir.path().join(OUTPUT).join("minimal-read-summary.csv");
    let text = fs::read_to_string(summary).unwrap();
    // Four files with six reads each, sorted by path.
    assert_eq!(text.lines().count(), 5);
    assert!(text.lines().skip(1).all(|l| l.ends_with(",6")));
}

#[test]
fn test_read_summary_default() {
    initiate_raw_cmd!(cmd, tmp_dir, QC_DIR);
    cmd.assert().success();
    let output = tmp_dir.path().join(OUTPUT);
    let header = header(&output.join("default-read-summary.csv"));
    assert!(header.starts_with("file_path,file_name,read_count,base_count"));
    assert!(header.ends_with("mean_read_QScore,low_quality_reads,phred_offset"));
    assert!(!output.join("read-report.json").exists());
}

#[test]
fn test_read_summary_complete() {
    initiate_raw_cmd!(cmd, tmp_dir, QC_DIR);
    cmd.args(["--mode", "complete", "--prefix", "qc"])
        .assert()
        .success();
    let output = tmp_dir.path().join(OUTPUT);
    let header = header(&output.join("qc_default-read-summary.csv"));
    assert!(header.ends_with("polyG_reads_pct,overrepresented_count"));
    let positions = zip_entries(&output.join("qc_read-pos-summary.zip"));
    assert!(positions.contains(&"qc_adapter.csv".to_string()));
    assert!(positions.contains(&"qc_bgzf.fastq.csv".to_string()));
    let distributions = zip_entries(&output.join("qc_read-dist-summary.zip"));
    assert!(distributions.contains(&"qc_adapter/adapters.csv".to_string()));
    assert!(distributions.contains(&"qc_phred64/duplication.csv".to_string()));
    let report = fs::read_to_string(output.join("qc_read-report.json")).unwrap();
    assert!(report.contains("\"phred_offset\":64"));
    assert!(report.contains("Illumina Universal"));
}

#[test]
fn test_read_summary_phred_override() {
    initiate_raw_cmd!(cmd, tmp_dir, QC_DIR);
    // The Phred+33 files have scores below '@', which Phred+64 cannot encode.
    cmd.args(["--phred", "64"]).assert().failure();
}

#[test]
fn test_read_summary_malformed() {
    initiate_raw_cmd!(cmd, tmp_dir, "crates/segul-core/tests/files/raw_invalid");
    cmd.assert().failure();
}

#[test]
fn test_read_summary_dup_sketch_size() {
    for size in ["0", "2"] {
        initiate_raw_cmd!(cmd, tmp_dir, QC_DIR);
        cmd.args(["--mode", "complete", "--dup-sketch-size", size])
            .assert()
            .success();
        let report =
            fs::read_to_string(tmp_dir.path().join(OUTPUT).join("read-report.json")).unwrap();
        // Six reads with four distinct sequences: a sketch of two is an estimate.
        let exact = if size == "0" { "true" } else { "false" };
        assert!(report.contains(&format!("\"exact\":{}", exact)));
    }
}
