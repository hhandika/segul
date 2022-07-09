mod utils;

use std::env;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::Path;

use glob::glob;
use predicates::Predicate;

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

#[test]
fn test_summary() {
    initiate_cmd!(cmd, "summary", "tests/files/long-aln", tmp_dir);
    cmd.assert().success();
    let pred = predicates::path::is_dir();
    let output_dir = tmp_dir.path().join("SEGUL-Summary");
    let locus_path = output_dir.join("locus_summary.csv");
    let taxon_path = output_dir.join("taxon_summary.csv");

    // Check column counts
    let locus_cols = count_header_cols(&locus_path);
    let taxon_cols = count_header_cols(&taxon_path);

    // Check taxon row counts
    let loci = parse_csv(&locus_path);
    let taxon = parse_csv(&taxon_path);

    assert!(pred.eval(&output_dir));
    assert_eq!(4, loci.len());
    assert_eq!(4, taxon.len());
    assert_eq!(locus_cols, 34); // DNA chars = 18
    assert_eq!(taxon_cols, 26); // Cols = 8

    tmp_dir.close().unwrap();
}

#[test]
fn test_locus_summary() {
    initiate_cmd!(cmd, "summary", "tests/files/long-aln", tmp_dir);
    cmd.arg("--per-locus").assert().success();
    let pred = predicates::path::is_dir();
    let output_dir = tmp_dir.path().join("SEGUL-Summary");
    let files = glob(&format!("{}/*.csv", output_dir.display()))
        .expect("Failed globbing files")
        .filter_map(|ok| ok.ok())
        .collect::<Vec<_>>();
    let cols = count_header_cols(&files[0]);
    assert!(pred.eval(&output_dir));
    assert_eq!(4, files.len());
    assert_eq!(cols, 19);
    tmp_dir.close().unwrap();
}

#[test]
fn test_summary_aa() {
    initiate_cmd!(cmd, "summary", "tests/files/concat-aa", tmp_dir);
    cmd.arg("--datatype").arg("aa").assert().success();
    let pred = predicates::path::is_dir();
    let output_dir = tmp_dir.path().join("SEGUL-Summary");
    let locus_path = output_dir.join("locus_summary.csv");
    let taxon_path = output_dir.join("taxon_summary.csv");

    // Check column counts
    let locus_cols = count_header_cols(&locus_path);
    let taxon_cols = count_header_cols(&taxon_path);

    // Check taxon row counts
    let loci = parse_csv(&locus_path);
    let taxon = parse_csv(&taxon_path);

    assert!(pred.eval(&output_dir));
    assert_eq!(4, loci.len());
    assert_eq!(3, taxon.len());
    assert_eq!(locus_cols, 44); // AA chars = 31
    assert_eq!(taxon_cols, 36); // Cols = 5

    tmp_dir.close().unwrap();
}
