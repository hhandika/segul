mod utils;

use std::env;
use std::fs;

use segul_core::parser::tree::{Rootedness, TreeFmt, parse_trees_file};
use segul_core::report::{RunReport, TaskKind};

#[test]
fn test_tree_estimate_distance_methods_and_formats() {
    let tmp = utils::create_tmp_dir().unwrap();
    let input = env::current_dir()
        .unwrap()
        .join("tests/files/tree_estimate/dna_6taxa.fas");

    for method in ["nj", "bionj", "upgma"] {
        let out_dir = tmp.path().join(format!("out_{method}"));
        let mut cmd = utils::segul(tmp.path());
        cmd.args(["tree", "infer", "-i"])
            .arg(&input)
            .args(["-m", method, "--model", "k2p", "-o"])
            .arg(&out_dir);
        cmd.assert().success();

        let tree_file = out_dir.join("dna_6taxa.nwk");
        assert!(tree_file.is_file());
        assert!(out_dir.join("tree_infer_summary.csv").is_file());
        assert!(out_dir.join("report.json").is_file());

        let trees = parse_trees_file(&tree_file, TreeFmt::Newick).unwrap();
        assert_eq!(trees.len(), 1);
        assert_eq!(trees[0].tip_count(), 6);

        let report = RunReport::from_path(&out_dir.join("report.json")).unwrap();
        assert_eq!(report.kind, TaskKind::TreeInfer);
    }

    // Check NEXUS, PhyloXML, NeXML, and gzip outputs
    let nex_dir = tmp.path().join("out_nex");
    utils::segul(tmp.path())
        .args(["tree", "infer", "-i"])
        .arg(&input)
        .args([
            "-m",
            "bionj",
            "-F",
            "nexus",
            "--no-translate",
            "--gzip",
            "-o",
        ])
        .arg(&nex_dir)
        .assert()
        .success();
    let nex_gz = nex_dir.join("dna_6taxa.nex.gz");
    assert!(nex_gz.is_file());
    let nex_trees = parse_trees_file(&nex_gz, TreeFmt::Auto).unwrap();
    assert_eq!(nex_trees.len(), 1);
    assert_eq!(nex_trees[0].tip_count(), 6);

    let pxml_dir = tmp.path().join("out_pxml");
    utils::segul(tmp.path())
        .args(["tree", "infer", "-i"])
        .arg(&input)
        .args(["-m", "nj", "-F", "phyloxml", "-o"])
        .arg(&pxml_dir)
        .assert()
        .success();
    let pxml_trees =
        parse_trees_file(&pxml_dir.join("dna_6taxa.phyloxml"), TreeFmt::PhyloXml).unwrap();
    assert_eq!(pxml_trees[0].tip_count(), 6);

    let nxml_dir = tmp.path().join("out_nxml");
    utils::segul(tmp.path())
        .args(["tree", "infer", "-i"])
        .arg(&input)
        .args(["-m", "nj", "-F", "nexml", "-o"])
        .arg(&nxml_dir)
        .assert()
        .success();
    let nxml_trees = parse_trees_file(&nxml_dir.join("dna_6taxa.nexml"), TreeFmt::NeXml).unwrap();
    assert_eq!(nxml_trees[0].tip_count(), 6);
}

#[test]
fn test_tree_estimate_parsimony_aa_and_bootstrap_determinism() {
    let tmp = utils::create_tmp_dir().unwrap();
    let dna_input = env::current_dir()
        .unwrap()
        .join("tests/files/tree_estimate/dna_6taxa.fas");
    let aa_input = env::current_dir()
        .unwrap()
        .join("tests/files/tree_estimate/aa_5taxa.fas");

    let mp1 = tmp.path().join("mp1");
    let mp2 = tmp.path().join("mp2");

    for out_dir in [&mp1, &mp2] {
        utils::segul(tmp.path())
            .args(["tree", "infer", "-i"])
            .arg(&dna_input)
            .args([
                "-m",
                "mp",
                "--search",
                "exact",
                "--mp-branch-lengths",
                "acctran",
                "-b",
                "20",
                "--seed",
                "123",
                "--outgroup",
                "T5,T6",
                "-o",
            ])
            .arg(out_dir)
            .assert()
            .success();
    }

    let nwk1 = fs::read_to_string(mp1.join("dna_6taxa.nwk")).unwrap();
    let nwk2 = fs::read_to_string(mp2.join("dna_6taxa.nwk")).unwrap();
    assert_eq!(nwk1, nwk2);

    let trees = parse_trees_file(&mp1.join("dna_6taxa.nwk"), TreeFmt::Newick).unwrap();
    assert_eq!(trees[0].tip_count(), 6);
    let supports: Vec<f64> = trees[0]
        .preorder()
        .filter(|&id| !trees[0].is_leaf(id) && !trees[0].is_root(id))
        .filter_map(|id| trees[0].support_value(id))
        .collect();
    assert!(!supports.is_empty());

    // Amino acid parsimony and Poisson distance
    let aa_mp = tmp.path().join("aa_mp");
    utils::segul(tmp.path())
        .args(["tree", "infer", "-i"])
        .arg(&aa_input)
        .args(["--datatype", "aa", "-m", "mp", "--midpoint", "-o"])
        .arg(&aa_mp)
        .assert()
        .success();
    let aa_trees = parse_trees_file(&aa_mp.join("aa_5taxa.nwk"), TreeFmt::Newick).unwrap();
    assert_eq!(aa_trees[0].tip_count(), 5);

    let aa_dist = tmp.path().join("aa_dist");
    utils::segul(tmp.path())
        .args(["tree", "infer", "-i"])
        .arg(&aa_input)
        .args(["--datatype", "aa", "-m", "nj", "--model", "poisson", "-o"])
        .arg(&aa_dist)
        .assert()
        .success();
    assert!(aa_dist.join("aa_5taxa.nwk").is_file());
}

#[test]
fn test_tree_estimate_batch_single_file_verify_replay_and_methods() {
    let tmp = utils::create_tmp_dir().unwrap();
    let dir = env::current_dir()
        .unwrap()
        .join("tests/files/tree_estimate");
    let out1 = tmp.path().join("batch_single");
    let out_replay = tmp.path().join("replayed");

    // Batch directory of DNA files (filter only dna_* or pass explicit files)
    let f1 = dir.join("dna_6taxa.fas");
    let f2 = dir.join("dna_locus2.nex");

    utils::segul(tmp.path())
        .args(["tree", "infer", "-i"])
        .arg(&f1)
        .args([
            "--checksum",
            "all",
            "-m",
            "bionj",
            "--model",
            "jc69",
            "--midpoint",
            "--prefix",
            "myrun",
            "--single-file",
            "-o",
        ])
        .arg(&out1)
        .assert()
        .success();

    let combined = out1.join("myrun.nwk");
    assert!(combined.is_file());
    assert!(out1.join("myrun_summary.csv").is_file());
    let rpath = out1.join("myrun_report.json");
    assert!(rpath.is_file());

    let trees = parse_trees_file(&combined, TreeFmt::Newick).unwrap();
    assert_eq!(trees.len(), 1);
    assert_eq!(trees[0].metadata.rootedness, Rootedness::Unknown); // Newick doesn't encode [&R]

    // Verify report checksums
    utils::segul(tmp.path())
        .args(["report", "verify"])
        .arg(&rpath)
        .assert()
        .success();

    // Replay report
    utils::segul(tmp.path())
        .args(["run"])
        .arg(&rpath)
        .args(["--verify-inputs", "-o"])
        .arg(&out_replay)
        .assert()
        .success();
    assert!(out_replay.join("myrun.nwk").is_file());
    assert_eq!(
        fs::read_to_string(&combined).unwrap(),
        fs::read_to_string(out_replay.join("myrun.nwk")).unwrap()
    );

    // Methods generation
    utils::segul(tmp.path())
        .args(["report", "methods"])
        .arg(&rpath)
        .assert()
        .success()
        .stdout(predicates::str::contains("Tree inference"));

    // Multi-file directory check using a temporary DNA-only folder
    let dna_dir = tmp.path().join("dna_only");
    fs::create_dir_all(&dna_dir).unwrap();
    fs::copy(&f1, dna_dir.join("dna_6taxa.fas")).unwrap();
    fs::copy(&f2, dna_dir.join("dna_locus2.nex")).unwrap();

    let batch_out = tmp.path().join("batch_multi");
    utils::segul(tmp.path())
        .args(["tree", "infer", "-d"])
        .arg(&dna_dir)
        .args(["-m", "upgma", "--single-file", "-F", "nexus", "-o"])
        .arg(&batch_out)
        .assert()
        .success();

    let nex_trees = parse_trees_file(&batch_out.join("trees.nex"), TreeFmt::Nexus).unwrap();
    assert_eq!(nex_trees.len(), 2);
}

#[test]
fn test_tree_estimate_error_cases() {
    let tmp = utils::create_tmp_dir().unwrap();
    let dna_input = env::current_dir()
        .unwrap()
        .join("tests/files/tree_estimate/dna_6taxa.fas");

    // Invalid model for DNA
    utils::segul(tmp.path())
        .args(["tree", "infer", "-i"])
        .arg(&dna_input)
        .args(["-m", "nj", "--model", "poisson", "-o"])
        .arg(tmp.path().join("err1"))
        .assert()
        .failure();

    // Unknown outgroup taxon
    utils::segul(tmp.path())
        .args(["tree", "infer", "-i"])
        .arg(&dna_input)
        .args(["-m", "nj", "--outgroup", "NO_SUCH_TAXON", "-o"])
        .arg(tmp.path().join("err2"))
        .assert()
        .failure();

    // Non-monophyletic outgroup
    utils::segul(tmp.path())
        .args(["tree", "infer", "-i"])
        .arg(&dna_input)
        .args(["-m", "nj", "--outgroup", "T1,T4", "-o"])
        .arg(tmp.path().join("err3"))
        .assert()
        .failure();
}
