mod utils;

#[test]
fn test_segul() {
    utils::segul().arg("-V").assert().success();
}
