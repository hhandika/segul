mod utils;

use std::env;
use std::path::Path;

use segul_core::helper::types::DataType;
use segul_core::parser::nexus::Nexus;
use segul_core::parser::tree::{TreeFmt, parse_trees_file};

#[test]
fn test_tree_convert_newick_to_nexus_and_back() {
    let tmp = utils::create_tmp_dir().unwrap();
    let input = env::current_dir()
        .unwrap()
        .join("tests/files/trees/simple.nwk");
    let nex_out = tmp.path().join("nex_out");

    let mut cmd = utils::segul(tmp.path());
    cmd.args(["tree", "convert", "-i"])
        .arg(&input)
        .args(["-f", "newick", "-F", "nexus", "-o"])
        .arg(&nex_out);
    cmd.assert().success();

    let nex_file = nex_out.join("simple.nex");
    assert!(nex_file.is_file());
    assert!(nex_out.join("report.json").is_file());
    let nex_text = std::fs::read_to_string(&nex_file).unwrap();
    assert!(nex_text.contains("TRANSLATE"));

    let trees = parse_trees_file(&nex_file, TreeFmt::Nexus).unwrap();
    assert_eq!(trees.len(), 1);
    assert_eq!(trees[0].tip_count(), 5);

    // Convert NEXUS back to Newick
    let nwk_out = tmp.path().join("nwk_out");
    let mut cmd2 = utils::segul(tmp.path());
    cmd2.args(["tree", "convert", "-i"])
        .arg(&nex_file)
        .args(["-F", "newick", "-o"])
        .arg(&nwk_out);
    cmd2.assert().success();

    let nwk_file = nwk_out.join("simple.nwk");
    assert!(nwk_file.is_file());
    let reparsed = parse_trees_file(&nwk_file, TreeFmt::Newick).unwrap();
    assert_eq!(reparsed.len(), 1);
    assert_eq!(
        reparsed[0]
            .leaves()
            .filter_map(|id| reparsed[0].label(id))
            .collect::<Vec<_>>(),
        vec!["Taxon_A", "Taxon_B", "Taxon_C", "Taxon_D", "Taxon_E"]
    );
}

#[test]
fn test_tree_convert_nexus_no_translate_and_strip_annotations() {
    let tmp = utils::create_tmp_dir().unwrap();
    let input = env::current_dir()
        .unwrap()
        .join("tests/files/trees/beast_annotated.trees");
    let out_dir = tmp.path().join("out");

    let mut cmd = utils::segul(tmp.path());
    cmd.args(["tree", "convert", "-i"])
        .arg(&input)
        .args(["-F", "nexus", "--no-translate", "--strip-annotations", "-o"])
        .arg(&out_dir);
    cmd.assert().success();

    let nex_file = out_dir.join("beast_annotated.nex");
    let text = std::fs::read_to_string(&nex_file).unwrap();
    assert!(!text.contains("TRANSLATE"));
    assert!(!text.contains("height_95%_HPD"));

    let trees = parse_trees_file(&nex_file, TreeFmt::Nexus).unwrap();
    assert_eq!(trees.len(), 1);
    assert_eq!(trees[0].tip_count(), 3);
}

#[test]
fn test_tree_convert_enewick_and_gzip() {
    let tmp = utils::create_tmp_dir().unwrap();
    let input = env::current_dir()
        .unwrap()
        .join("tests/files/trees/extended.enwk");
    let out_dir = tmp.path().join("out");

    let mut cmd = utils::segul(tmp.path());
    cmd.args(["tree", "convert", "-i"])
        .arg(&input)
        .args(["-F", "enewick", "--gzip", "-o"])
        .arg(&out_dir);
    cmd.assert().success();

    let gz_file = out_dir.join("extended.enwk.gz");
    assert!(gz_file.is_file());
    let trees = parse_trees_file(&gz_file, TreeFmt::Auto).unwrap();
    assert_eq!(trees.len(), 1);
    assert!(trees[0].is_network());
    assert_eq!(trees[0].reticulation_groups()["#H1"].len(), 2);
}

#[test]
fn test_tree_convert_dir_batch_and_mixed_nexus_isolation() {
    let tmp = utils::create_tmp_dir().unwrap();
    let tree_dir = env::current_dir().unwrap().join("tests/files/trees");
    let out_dir = tmp.path().join("batch_out");

    let mut cmd = utils::segul(tmp.path());
    cmd.args(["tree", "convert", "-d"])
        .arg(&tree_dir)
        .args(["-F", "newick", "-o"])
        .arg(&out_dir);
    cmd.assert().success();

    assert!(out_dir.join("simple.nwk").is_file());
    assert!(out_dir.join("simple_trees.nwk").is_file());
    assert!(out_dir.join("multi.nwk").is_file());
    assert!(out_dir.join("edge_cases.nwk").is_file());
    assert!(out_dir.join("beast_annotated.nwk").is_file());
    assert!(out_dir.join("taxa_no_translate.nwk").is_file());
    assert!(out_dir.join("mixed_data_and_trees.nwk").is_file());
    assert!(out_dir.join("sample.nwk").is_file());
    assert!(out_dir.join("sample_nexml.nwk").is_file());
    assert!(out_dir.join("report.json").is_file());

    // Confirm existing sequence NEXUS parser still parses mixed_data_and_trees.nex cleanly.
    let mixed = Path::new("tests/files/trees/mixed_data_and_trees.nex");
    let mut seq_nex = Nexus::new(mixed, &DataType::Dna);
    seq_nex.parse();
    assert_eq!(seq_nex.matrix.len(), 3);
    assert_eq!(seq_nex.header.nchar, 6);
    assert_eq!(seq_nex.matrix.get("TaxonB"), Some(b"TGCATG".as_slice()));
}

#[test]
fn test_tree_convert_phyloxml_and_nexml() {
    let tmp = utils::create_tmp_dir().unwrap();
    let pxml_in = env::current_dir()
        .unwrap()
        .join("tests/files/trees/sample.phyloxml");
    let nxml_in = env::current_dir()
        .unwrap()
        .join("tests/files/trees/sample_nexml.nexml");

    let nxml_out = tmp.path().join("to_nexml");
    let mut cmd1 = utils::segul(tmp.path());
    cmd1.args(["tree", "convert", "-i"])
        .arg(&pxml_in)
        .args(["-F", "nexml", "--gzip", "-o"])
        .arg(&nxml_out);
    cmd1.assert().success();

    let nxml_gz = nxml_out.join("sample.nexml.gz");
    assert!(nxml_gz.is_file());
    let from_nxml_gz = parse_trees_file(&nxml_gz, TreeFmt::Auto).unwrap();
    assert_eq!(from_nxml_gz.len(), 1);
    assert_eq!(
        from_nxml_gz[0]
            .leaves()
            .filter_map(|id| from_nxml_gz[0].label(id))
            .collect::<Vec<_>>(),
        vec!["Homo & sapiens", "Pan troglodytes", "Gorilla"]
    );

    let pxml_out = tmp.path().join("to_phyloxml");
    let mut cmd2 = utils::segul(tmp.path());
    cmd2.args(["tree", "convert", "-i"])
        .arg(&nxml_in)
        .args(["-F", "phyloxml", "--strip-annotations", "-o"])
        .arg(&pxml_out);
    cmd2.assert().success();

    let pxml_file = pxml_out.join("sample_nexml.phyloxml");
    assert!(pxml_file.is_file());
    let from_pxml = parse_trees_file(&pxml_file, TreeFmt::PhyloXml).unwrap();
    assert_eq!(from_pxml.len(), 1);
    assert_eq!(
        from_pxml[0]
            .leaves()
            .filter_map(|id| from_pxml[0].label(id))
            .collect::<Vec<_>>(),
        vec!["Taxon_A", "Taxon_B", "Taxon_C"]
    );
}
