mod utils;

use tempdir::TempDir;

macro_rules! test_subcommand {
    ($func: ident, $subcommand: expr) => {
        #[test]
        fn $func() {
            let tmp_dir = TempDir::new(utils::DIR).unwrap();
            utils::segul(tmp_dir.path())
                .arg($subcommand)
                .arg("-h")
                .assert()
                .success();
            tmp_dir.close().unwrap();
        }
    };
}

#[test]
fn test_version() {
    use clap::crate_version;
    let version = crate_version!();
    utils::segul(utils::create_tmp_dir().unwrap().path())
        .arg("-V")
        .assert()
        .success()
        .stdout(predicates::str::contains(version));
}

test_subcommand! {test_concat_subcommand, "concat"}
test_subcommand! {test_convert_subcommand, "convert"}
test_subcommand! {test_extract_subcommand, "extract"}
test_subcommand! {test_filter_subcommand, "filter"}
test_subcommand! {test_id_subcommand, "id"}
test_subcommand! {test_partition_subcommand, "partition"}
test_subcommand! {test_remove_subcommand, "remove"}
test_subcommand! {test_rename_subcommand, "rename"}
test_subcommand! {test_split_subcommand, "split"}
test_subcommand! {test_summary_subcommand, "summary"}
test_subcommand! {test_translate_subcommand, "translate"}
