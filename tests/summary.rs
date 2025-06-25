mod utils;

use std::env;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::Path;

use glob::glob;
use predicates::Predicate;

macro_rules! generate_test {
    ($tmp_dir: ident, $loc_count: expr_2021, $taxon_count: expr_2021, $loc_cols: expr_2021, $taxon_cols: expr_2021) => {
        let pred = predicates::path::is_dir();
        let output_dir = $tmp_dir.path().join("Align-Summary");
        let locus_path = output_dir.join("locus_summary.csv");
        let taxon_path = output_dir.join("taxon_summary.csv");

        // Check column counts
        let locus_cols = count_header_cols(&locus_path);
        let taxon_cols = count_header_cols(&taxon_path);

        // Check taxon row counts
        let loci = parse_csv(&locus_path);
        let taxon = parse_csv(&taxon_path);

        assert!(pred.eval(&output_dir));
        assert_eq!($loc_count, loci.len());
        assert_eq!($taxon_count, taxon.len());
        assert_eq!($loc_cols, locus_cols); // DNA chars = 18
        assert_eq!($taxon_cols, taxon_cols); // Cols = 8

        $tmp_dir.close().unwrap();
    };
}

macro_rules! generate_locus_test {
    ($tmp_dir:ident, $fcount: expr_2021, $cols: expr_2021) => {
        let pred = predicates::path::is_dir();
        let output_dir = $tmp_dir.path().join("Align-Summary");
        let files = glob(&format!("{}/*.csv", output_dir.display()))
            .expect("Failed globbing files")
            .filter_map(|ok| ok.ok())
            .collect::<Vec<_>>();
        let cols = count_header_cols(&files[0]);
        assert!(pred.eval(&output_dir));
        assert_eq!($fcount, files.len());
        assert_eq!($cols, cols);
        $tmp_dir.close().unwrap();
    };
}

#[test]
fn test_summary() {
    initiate_cmd!(cmd, "align", "summary", "tests/files/long-aln", tmp_dir);
    cmd.assert().success();
    let locus_counts = 4;
    let taxon_counts = 4;
    let locus_cols = 34;
    let taxon_cols = 26;
    generate_test!(tmp_dir, locus_counts, taxon_counts, locus_cols, taxon_cols);
}

#[test]
fn test_summary_aa() {
    initiate_cmd!(cmd, "align", "summary", "tests/files/concat-aa", tmp_dir);
    cmd.arg("--datatype").arg("aa").assert().success();
    let locus_counts = 4;
    let taxon_counts = 3;
    let locus_cols = 44; // AA chars = 31
    let taxon_cols = 36; // Cols = 5
    generate_test!(tmp_dir, locus_counts, taxon_counts, locus_cols, taxon_cols);
}

#[test]
fn test_locus_summary() {
    initiate_cmd!(cmd, "align", "summary", "tests/files/long-aln", tmp_dir);
    cmd.arg("--per-locus").assert().success();
    let fcount = 4;
    let cols = 25;
    generate_locus_test!(tmp_dir, fcount, cols);
}

#[test]
fn test_locus_summary_aa() {
    initiate_cmd!(cmd, "align", "summary", "tests/files/concat-aa", tmp_dir);
    cmd.arg("--datatype")
        .arg("aa")
        .arg("--per-locus")
        .assert()
        .success();
    let fcount = 4;
    let cols = 35;
    generate_locus_test!(tmp_dir, fcount, cols);
}

fn parse_csv(fpath: &Path) -> Vec<String> {
    let file = File::open(fpath).expect("Unable to open file");
    let buff = BufReader::new(file);
    let mut result = Vec::new();
    buff.lines()
        .filter_map(|ok| ok.ok())
        .skip(1)
        .for_each(|line| {
            let parts: Vec<&str> = line.split(',').map(|e| e.trim()).collect();
            result.push(parts[0].to_string());
        });

    result
}

fn count_header_cols(fpath: &Path) -> usize {
    let file = File::open(fpath).expect("Unable to open file");
    let buff = BufReader::new(file);
    let mut result = 0;
    buff.lines()
        .filter_map(|ok| ok.ok())
        .take(1)
        .for_each(|line| {
            let parts: Vec<&str> = line.split(',').map(|e| e.trim()).collect();
            result = parts.len();
        });

    result
}
