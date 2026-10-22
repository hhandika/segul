mod utils;

use std::env;
use std::path::PathBuf;

const FILES: &str = "crates/segul-core/tests/files";

fn fixture(name: &str) -> PathBuf {
    env::current_dir().unwrap().join(FILES).join(name)
}

macro_rules! test_summary_plot {
    ($name: ident, [$($arg: expr_2021),+], $dir: expr_2021, $output: expr_2021, [$($chart: expr_2021),+]) => {
        #[test]
        fn $name() {
            let tmp_dir = utils::create_tmp_dir().unwrap();
            let mut cmd = utils::segul(tmp_dir.path());
            cmd.args([$($arg),+])
                .arg("-d")
                .arg(fixture($dir))
                .args(["-o", $output, "--plot"]);
            cmd.assert().success();
            let plots = tmp_dir.path().join($output).join("plots");
            $(assert!(plots.join($chart).is_file(), "Missing {}", $chart);)+
            tmp_dir.close().unwrap();
        }
    };
}

test_summary_plot!(
    test_align_summary_plot,
    ["align", "summary", "-f", "nexus", "--plot-format", "svg"],
    "alignments",
    "Align-Summary",
    ["completeness.svg", "pis-vs-length.svg", "taxon-loci.svg"]
);

test_summary_plot!(
    test_align_distance_plot,
    [
        "align",
        "distance",
        "-f",
        "nexus",
        "--long",
        "--cvd",
        "protanopia"
    ],
    "gappy",
    "Align-Distance",
    ["locus-mean.png", "matrix_gene_2.png", "locus-box.png"]
);

test_summary_plot!(
    test_read_summary_plot,
    ["read", "summary", "--mode", "complete"],
    "raw",
    "Read-Summary",
    [
        "read-count.png",
        "base-composition.png",
        "pos-qscore_read_1.png",
        "pos-composition_read_1.png"
    ]
);

test_summary_plot!(
    test_contig_summary_plot,
    ["contig", "summary"],
    "contigs",
    "Contig-Summary",
    ["nstats.png", "count-vs-n50.png"]
);

test_summary_plot!(
    test_maf_summary_plot,
    ["maf", "summary", "--prefix", "run"],
    "maf",
    "Maf-Summary",
    ["run_species-coverage.png", "run_file-coverage.png"]
);

#[test]
fn test_summary_without_plot() {
    let tmp_dir = utils::create_tmp_dir().unwrap();
    let mut cmd = utils::segul(tmp_dir.path());
    cmd.args(["contig", "summary", "-d"])
        .arg(fixture("contigs"))
        .args(["-o", "Contig-Summary"]);
    cmd.assert().success();
    assert!(!tmp_dir.path().join("Contig-Summary").join("plots").exists());
    tmp_dir.close().unwrap();
}

#[test]
fn test_plot_options_require_plot() {
    let tmp_dir = utils::create_tmp_dir().unwrap();
    let mut cmd = utils::segul(tmp_dir.path());
    cmd.args(["contig", "summary", "-d"])
        .arg(fixture("contigs"))
        .args(["--plot-format", "svg"]);
    cmd.assert().failure();
    tmp_dir.close().unwrap();
}

#[test]
fn test_plot_conflicts_with_per_locus() {
    let tmp_dir = utils::create_tmp_dir().unwrap();
    let mut cmd = utils::segul(tmp_dir.path());
    cmd.args(["align", "summary", "-f", "nexus", "-d"])
        .arg(fixture("alignments"))
        .args(["--per-locus", "--plot"]);
    cmd.assert().failure();
    tmp_dir.close().unwrap();
}

#[test]
fn test_plot_rejects_pdf() {
    let tmp_dir = utils::create_tmp_dir().unwrap();
    let mut cmd = utils::segul(tmp_dir.path());
    cmd.args(["contig", "summary", "-d"])
        .arg(fixture("contigs"))
        .args(["--plot", "--plot-format", "pdf"]);
    cmd.assert().failure();
    tmp_dir.close().unwrap();
}
