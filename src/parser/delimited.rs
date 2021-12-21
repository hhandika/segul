use std::ffi::OsStr;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::Path;

pub fn parse_delimited_text(fpath: &Path) -> Vec<(String, String)> {
    let file = File::open(fpath).expect("Unable to open file");
    let buff = BufReader::new(file);
    let mut result = Vec::new();
    let ext: &str = fpath
        .extension()
        .and_then(OsStr::to_str)
        .expect("Failed parsing extension");
    assert!(ext == "tsv" || ext == "csv");

    buff.lines()
        .filter_map(|ok| ok.ok())
        .skip(1)
        .for_each(|line| {
            let parts: Vec<&str> = match ext {
                "csv" => line.split(',').map(|e| e.trim()).collect(),
                "tsv" => line.split_whitespace().map(|e| e.trim()).collect(),
                _ => panic!("Unsupported file extension"),
            };
            assert_eq!(
                parts.len(),
                2,
                "Failed parsing delimited file. Expected 2 columns, found {}",
                parts.len()
            );
            result.push((parts[0].to_string(), parts[1].to_string()));
        });

    result
}

#[cfg(test)]
mod test {
    use super::*;

    macro_rules! test_parse_delimited_text {
        ($name:ident, $fpath:expr, $expected_or:expr,$expected_dest:expr) => {
            #[test]
            fn $name() {
                let result = parse_delimited_text($fpath);
                assert_eq!(result.len(), 1);
                assert_eq!(result[0].0, $expected_or);
                assert_eq!(result[0].1, $expected_dest);
            }
        };
    }

    test_parse_delimited_text!(
        test_parse_csv,
        Path::new("test_files/delimited.csv"),
        "sequence_old",
        "sequence_new"
    );

    test_parse_delimited_text!(
        test_parse_tsv,
        Path::new("test_files/delimited.tsv"),
        "sequence_old",
        "sequence_new"
    );
}
