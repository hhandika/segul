mod utils;

use std::env;
use std::fs;

use segul_core::parser::tree::{Rootedness, TreeFmt, parse_trees_file};
use segul_core::report::{ReportStatistics, RunReport, TaskKind};

#[test]
fn test_tree_summary_defaults_and_beast_burnin_ca_nexus() {
    let tmp = utils::create_tmp_dir().unwrap();
    let input = env::current_dir()
        .unwrap()
        .join("tests/files/trees/posterior_beast.trees");
    let out = tmp.path().join("beast_mcc");

    // Default method is mcc, default output-format is newick.
    utils::segul(tmp.path())
        .args(["tree", "summary", "-i"])
        .arg(&input)
        .args([
            "--burnin",
            "20%",
            "--branch-lengths",
            "ca",
            "--checksum",
            "all",
            "-o",
        ])
        .arg(&out)
        .assert()
        .success();

    let nwk_path = out.join("tree_summary.nwk");
    assert!(nwk_path.is_file());
    assert!(out.join("tree_summary_clades.csv").is_file());
    assert!(out.join("tree_summary_topologies.csv").is_file());
    let rpath = out.join("report.json");
    assert!(rpath.is_file());

    let trees = parse_trees_file(&nwk_path, TreeFmt::Newick).unwrap();
    assert_eq!(trees.len(), 1);
    assert_eq!(trees[0].tip_count(), 4);
    for id in trees[0].preorder() {
        if !trees[0].is_root(id) {
            assert!(trees[0].branch_length(id).unwrap() >= 0.0);
        }
    }
    let nwk_text = fs::read_to_string(&nwk_path).unwrap();
    assert!(nwk_text.contains("height_95%_HPD="));
    assert!(nwk_text.contains("rate="));
    assert!(nwk_text.contains("loc=north"));

    let report = RunReport::from_path(&rpath).unwrap();
    assert_eq!(report.kind, TaskKind::TreeSummary);
    let Some(ReportStatistics::TreeSummary(stats)) = report.typed_statistics().unwrap() else {
        panic!("Expected TreeSummary statistics");
    };
    assert_eq!(stats.trees_read, 10);
    assert_eq!(stats.trees_burned_in, 2);
    assert_eq!(stats.trees_summarized, 8);
    assert_eq!(stats.unique_topologies, 3);
    assert_eq!(stats.method, "mcc");
    assert!(stats.rooted);

    // Verify, replay, and methods on the report.
    utils::segul(tmp.path())
        .args(["report", "verify"])
        .arg(&rpath)
        .assert()
        .success();

    let replay_out = tmp.path().join("beast_replay");
    utils::segul(tmp.path())
        .args(["run"])
        .arg(&rpath)
        .args(["--verify-inputs", "-o"])
        .arg(&replay_out)
        .assert()
        .success();
    assert_eq!(
        fs::read_to_string(&nwk_path).unwrap(),
        fs::read_to_string(replay_out.join("tree_summary.nwk")).unwrap()
    );

    utils::segul(tmp.path())
        .args(["report", "methods"])
        .arg(&rpath)
        .assert()
        .success()
        .stdout(predicates::str::contains("Tree summary"));
}

#[test]
fn test_tree_summary_unrooted_methods_and_all_output_formats() {
    let tmp = utils::create_tmp_dir().unwrap();
    let input = env::current_dir()
        .unwrap()
        .join("tests/files/trees/bootstrap_unrooted.nwk");

    let expected_internal = [
        ("strict", 1usize),
        ("mrc", 2usize),
        ("mre", 3usize),
        ("mcc", 3usize),
        ("map", 3usize),
    ];

    for (method, exp_nodes) in expected_internal {
        let out = tmp.path().join(format!("out_{method}"));
        utils::segul(tmp.path())
            .args(["tree", "summary", "-i"])
            .arg(&input)
            .args(["-m", method, "--support-format", "percent", "-o"])
            .arg(&out)
            .assert()
            .success();

        let trees = parse_trees_file(&out.join("tree_summary.nwk"), TreeFmt::Newick).unwrap();
        let internal = trees[0]
            .preorder()
            .filter(|&id| !trees[0].is_leaf(id) && !trees[0].is_root(id))
            .count();
        assert_eq!(internal, exp_nodes, "method={method}");
    }

    // Check NEXUS, PhyloXML, NeXML, and gzip + outgroup rooting
    let nex_out = tmp.path().join("out_nex");
    utils::segul(tmp.path())
        .args(["tree", "summary", "-i"])
        .arg(&input)
        .args([
            "-m",
            "mrc",
            "--outgroup",
            "A,B",
            "-F",
            "nexus",
            "--gzip",
            "--prefix",
            "boot",
            "-o",
        ])
        .arg(&nex_out)
        .assert()
        .success();
    let nex_gz = nex_out.join("boot.nex.gz");
    assert!(nex_gz.is_file());
    let nex_trees = parse_trees_file(&nex_gz, TreeFmt::Auto).unwrap();
    assert_eq!(nex_trees[0].metadata.rootedness, Rootedness::Rooted);
    assert_eq!(nex_trees[0].child_count(nex_trees[0].root()), 2);

    for (fmt, ext) in [("phyloxml", "phyloxml"), ("nexml", "nexml")] {
        let fmt_out = tmp.path().join(format!("out_{fmt}"));
        utils::segul(tmp.path())
            .args(["tree", "summary", "-i"])
            .arg(&input)
            .args(["-F", fmt, "-o"])
            .arg(&fmt_out)
            .assert()
            .success();
        let parsed =
            parse_trees_file(&fmt_out.join(format!("tree_summary.{ext}")), TreeFmt::Auto).unwrap();
        assert_eq!(parsed[0].tip_count(), 6);
    }
}

#[test]
fn test_tree_summary_partial_taxa_target_replay_and_errors() {
    let tmp = utils::create_tmp_dir().unwrap();
    let gene_trees = env::current_dir()
        .unwrap()
        .join("tests/files/trees/gene_trees_partial.nwk");
    let target = env::current_dir()
        .unwrap()
        .join("tests/files/trees/species_target.nwk");

    // Partial-taxa without flag must error
    utils::segul(tmp.path())
        .args(["tree", "summary", "-i"])
        .arg(&gene_trees)
        .args(["-o"])
        .arg(tmp.path().join("err_partial"))
        .assert()
        .failure();

    // Target method without --target must error
    utils::segul(tmp.path())
        .args(["tree", "summary", "-i"])
        .arg(&gene_trees)
        .args(["-m", "target", "-o"])
        .arg(tmp.path().join("err_no_target"))
        .assert()
        .failure();

    // Target + --allow-partial-taxa computes exact decisive gCF support and replays cleanly
    let out = tmp.path().join("partial_target");
    utils::segul(tmp.path())
        .args(["tree", "summary", "-i"])
        .arg(&gene_trees)
        .args(["-m", "target", "--target"])
        .arg(&target)
        .args(["--unrooted", "--allow-partial-taxa", "-o"])
        .arg(&out)
        .assert()
        .success();

    let clades_csv = fs::read_to_string(out.join("tree_summary_clades.csv")).unwrap();
    assert!(clades_csv.contains(",2,3,"));

    let replay_out = tmp.path().join("partial_replay");
    utils::segul(tmp.path())
        .args(["run"])
        .arg(out.join("report.json"))
        .args(["-o"])
        .arg(&replay_out)
        .assert()
        .success();
    assert_eq!(
        fs::read_to_string(out.join("tree_summary.nwk")).unwrap(),
        fs::read_to_string(replay_out.join("tree_summary.nwk")).unwrap()
    );
}
