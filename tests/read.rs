mod utils;

use std::env;
use std::fs;
use std::path::{Path, PathBuf};

use segul_core::stats::fastq::qc::QcProfile;
use segul_core::stats::fastq::report::ReadReport;

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
fn test_read_summary_basic() {
    for mode in ["basic", "default"] {
        initiate_raw_cmd!(cmd, tmp_dir, QC_DIR);
        cmd.args(["--mode", mode]).assert().success();
        let output = tmp_dir.path().join(OUTPUT);
        let header = header(&output.join("read-summary.csv"));
        assert!(header.starts_with("file_path,file_name,read_count,base_count"));
        assert!(header.contains("effective_mean_QScore,mean_error_rate"));
        assert!(header.contains("mean_complexity,low_complexity_reads_pct"));
        assert!(
            header.ends_with("platform,chastity_failed_reads,chastity_failed_pct,run_group_count")
        );
        assert!(!output.join("read-dist-summary.zip").exists());
        assert!(output.join("report.json").is_file());
    }
}

#[test]
fn test_read_summary_default_is_complete() {
    initiate_raw_cmd!(cmd, tmp_dir, QC_DIR);
    cmd.args(["--prefix", "qc"]).assert().success();
    let output = tmp_dir.path().join(OUTPUT);
    let header = header(&output.join("qc_read-summary.csv"));
    assert!(
        header.ends_with("QScore_tail_drop,effective_Q20_span,effective_Q30_span,head_bias_bp")
    );
    let positions = zip_entries(&output.join("qc_read-pos-summary.zip"));
    assert!(positions.contains(&"qc_adapter.csv".to_string()));
    assert!(positions.contains(&"qc_bgzf.csv".to_string()));
    let distributions = zip_entries(&output.join("qc_read-dist-summary.zip"));
    assert!(distributions.contains(&"qc_adapter/adapters.csv".to_string()));
    assert!(distributions.contains(&"qc_adapter/relative-pos.csv".to_string()));
    assert!(distributions.contains(&"qc_adapter/length-quality.csv".to_string()));
    assert!(distributions.contains(&"qc_adapter/expected-errors.csv".to_string()));
    assert!(distributions.contains(&"qc_adapter/complexity.csv".to_string()));
    assert!(distributions.contains(&"qc_phred64/duplication.csv".to_string()));
    let report = fs::read_to_string(output.join("qc_report.json")).unwrap();
    let parsed = ReadReport::from_json(&report).unwrap();
    assert!(parsed.samples.iter().any(|s| s.metrics.phred_offset == 64));
    assert!(report.contains("Illumina Universal"));
    assert!(
        parsed
            .samples
            .iter()
            .any(|s| s.metrics.effective_mean_qscore > 0.0)
    );
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
fn test_read_summary_qc_profile_and_pairs() {
    initiate_raw_cmd!(cmd, tmp_dir, "crates/segul-core/tests/files/raw");
    cmd.args(["--qc-profile", "amplicon"]).assert().success();
    let output = tmp_dir.path().join(OUTPUT);
    let report = fs::read_to_string(output.join("report.json")).unwrap();
    let parsed = ReadReport::from_json(&report).unwrap();
    assert_eq!(parsed.qc_profile, QcProfile::Amplicon);
    assert!(
        parsed
            .pairs
            .iter()
            .any(|p| p.pair_name == "read" && p.in_sync)
    );
    assert!(parsed.cohort.is_some());
    let distributions = zip_entries(&output.join("read-dist-summary.zip"));
    assert!(distributions.contains(&"pairs.csv".to_string()));
}

#[test]
fn test_read_summary_dup_size() {
    for size in ["0", "2"] {
        initiate_raw_cmd!(cmd, tmp_dir, QC_DIR);
        cmd.args(["--mode", "complete", "--dup-size", size])
            .assert()
            .success();
        let report = fs::read_to_string(tmp_dir.path().join(OUTPUT).join("report.json")).unwrap();
        let parsed = ReadReport::from_json(&report).unwrap();
        // Six reads with four distinct sequences: a sketch of two is an estimate.
        let exact = size == "0";
        assert!(
            parsed
                .samples
                .iter()
                .all(|s| s.duplication.as_ref().is_some_and(|d| d.exact == exact))
        );
    }
}
