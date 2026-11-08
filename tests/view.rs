mod utils;

#[test]
fn test_view_help() {
    utils::segul(utils::create_tmp_dir().unwrap().path())
        .arg("view")
        .arg("align")
        .arg("--help")
        .assert()
        .success()
        .stdout(predicates::str::contains("--theme"));
}

#[test]
fn test_view_missing_file() {
    utils::segul(utils::create_tmp_dir().unwrap().path())
        .arg("view")
        .arg("align")
        .arg("-i")
        .arg("missing.fas")
        .assert()
        .failure()
        .stderr(predicates::str::contains("missing.fas"));
}

#[test]
fn test_view_summary_help() {
    utils::segul(utils::create_tmp_dir().unwrap().path())
        .arg("view")
        .arg("summary")
        .arg("--help")
        .assert()
        .success()
        .stdout(predicates::str::contains("alignment-report.json"))
        .stdout(predicates::str::contains("distance-report.json"));
}

#[test]
fn test_view_summary_missing_file() {
    utils::segul(utils::create_tmp_dir().unwrap().path())
        .arg("view")
        .arg("summary")
        .arg("-i")
        .arg("missing-report.json")
        .assert()
        .failure()
        .stderr(predicates::str::contains("missing-report.json"));
}

#[test]
fn test_view_summary_not_a_report() {
    utils::segul(utils::create_tmp_dir().unwrap().path())
        .arg("view")
        .arg("summary")
        .arg("-i")
        .arg(
            std::env::current_dir()
                .unwrap()
                .join("crates/segul-core/tests/files/plot/locus_summary.csv"),
        )
        .assert()
        .failure()
        .stderr(predicates::str::contains("not a segul summary report"));
}
