mod utils;

use std::env;
use std::fs;
use std::path::PathBuf;

use segul_core::report::RunReport;
use segul_core::report::verify::VerifyReport;

fn fixture_dir(name: &str) -> PathBuf {
    env::current_dir().unwrap().join("tests/files").join(name)
}

#[test]
fn test_cli_checksum_verify_and_run_replay() {
    let tmp = utils::create_tmp_dir().unwrap();
    let out1 = tmp.path().join("summary_1");
    let out2 = tmp.path().join("summary_2");

    // 1. Run align summary with --checksum all.
    utils::segul(tmp.path())
        .args(["--checksum", "all", "align", "summary", "-d"])
        .arg(fixture_dir("alignments"))
        .args(["-f", "nexus", "-o"])
        .arg(&out1)
        .assert()
        .success();

    let report_path = out1.join("report.json");
    assert!(report_path.is_file());
    let r1 = RunReport::from_path(&report_path).unwrap();
    assert_eq!(r1.schema_version, 2);
    assert!(r1.inputs.primary.iter().all(|f| f.sha256.is_some()));
    assert!(r1.outputs.files.iter().all(|f| f.sha256.is_some()));

    // 2. Verify passes cleanly (both text and --json).
    let json_out = utils::segul(tmp.path())
        .args(["report", "verify"])
        .arg(&report_path)
        .arg("--json")
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let parsed_verify = VerifyReport::from_json(std::str::from_utf8(&json_out).unwrap()).unwrap();
    assert!(parsed_verify.passed);
    assert_eq!(parsed_verify.mismatched, 0);

    // 3. Replay from report.json via `segul run`.
    utils::segul(tmp.path())
        .args(["run"])
        .arg(&report_path)
        .args(["-o"])
        .arg(&out2)
        .args(["--verify-inputs", "--checksum", "all"])
        .assert()
        .success();

    let r2 = RunReport::from_path(&out2.join("report.json")).unwrap();
    assert_eq!(r1.parameters, r2.parameters);
    assert_eq!(r1.statistics, r2.statistics);

    // 4. Tamper an output file -> `segul report verify` exits non-zero.
    let locus_csv = out1.join("locus_summary.csv");
    fs::write(&locus_csv, b"tampered\n").unwrap();
    utils::segul(tmp.path())
        .args(["report", "verify"])
        .arg(&report_path)
        .assert()
        .failure();
}

#[test]
fn test_cli_chained_workflow_and_methods() {
    let tmp = utils::create_tmp_dir().unwrap();
    let step1_out = tmp.path().join("filtered_50p");
    let step2_dir = tmp.path().join("concat");
    let wf_out = tmp.path().join("wf_out");
    let wf_file = tmp.path().join("segul-workflow.json");
    let methods_file = tmp.path().join("methods.md");

    // Step 1: filter
    utils::segul(tmp.path())
        .args(["--checksum", "outputs", "align", "filter", "-d"])
        .arg(fixture_dir("gappy"))
        .args(["-f", "nexus", "--percent", "0.5", "-o", "filtered"])
        .assert()
        .success();
    assert!(step1_out.join("report.json").is_file());

    // Step 2: concat from step 1 outputs
    utils::segul(tmp.path())
        .args(["--checksum", "outputs", "align", "concat", "-d"])
        .arg(&step1_out)
        .args([
            "-f", "nexus", "-F", "nexus", "-p", "nexus", "--prefix", "combined", "-o",
        ])
        .arg(&step2_dir)
        .assert()
        .success();

    let step2_report = step2_dir.join("combined_report.json");
    let r2 = RunReport::from_path(&step2_report).unwrap();
    assert_eq!(r2.run.parent_reports.len(), 1);

    // Export workflow and replay via `segul run`.
    utils::segul(tmp.path())
        .args(["report", "workflow"])
        .arg(&step2_report)
        .args(["--name", "cli-pipeline", "-o"])
        .arg(&wf_file)
        .assert()
        .success();

    utils::segul(tmp.path())
        .args(["run"])
        .arg(&wf_file)
        .args(["-o"])
        .arg(&wf_out)
        .assert()
        .success();
    assert!(wf_out.join("01_align_filter").join("report.json").is_file());
    assert!(
        wf_out
            .join("02_align_concat")
            .join("combined_report.json")
            .is_file()
    );

    // Export methods in Markdown and plain text.
    utils::segul(tmp.path())
        .args(["report", "methods"])
        .arg(&step2_report)
        .args(["--format", "markdown", "-o"])
        .arg(&methods_file)
        .assert()
        .success();
    let md = fs::read_to_string(&methods_file).unwrap();
    assert!(md.contains("SEGUL"));
    assert!(md.contains("1. **"));
    assert!(md.contains("2. **"));
}
