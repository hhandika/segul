use std::collections::BTreeMap;

fn index_sites(mat: BTreeMap<String, String>) -> BTreeMap<usize, String> {
    let mut index: BTreeMap<usize, String> = BTreeMap::new();
    mat.values().for_each(|seq| {
        seq.chars().enumerate().for_each(|(idx, dna)| {
            if index.contains_key(&idx) {
                if let Some(value) = index.get_mut(&idx) {
                    match dna {
                        '-' | 'N' | '?' => (),
                        _ => value.push(dna),
                    }
                }
            } else {
                match dna {
                    '-' | 'N' | '?' => (),
                    _ => {
                        index.insert(idx, dna.to_string());
                    }
                }
            }
        })
    });

    index
}

fn get_parsimony_informative(matrix: &BTreeMap<usize, String>) -> usize {
    let mut parsim: usize = 0;
    matrix.values().for_each(|site| {
        let n_pattern = get_pattern(&site);
        println!("Pattern: {}", n_pattern);
        if n_pattern > 1 {
            parsim += 1
        }
    });
    parsim
}

fn get_pattern(site: &str) -> usize {
    let mut uniques: Vec<char> = site.chars().collect();
    uniques.sort();
    uniques.dedup();
    let mut pattern = 0;
    uniques.iter().for_each(|c| {
        let n_pattern: String = site.matches(&c.to_string()).collect();
        if n_pattern.len() >= 2 {
            pattern += 1;
        }
    });

    pattern
}

#[cfg(test)]
mod test {
    use super::*;

    fn get_matrix(id: &[&str], seq: &[&str]) -> BTreeMap<String, String> {
        let mut matrix = BTreeMap::new();
        id.iter().zip(seq.iter()).for_each(|(i, s)| {
            matrix.insert(i.to_string(), s.to_string());
        });

        matrix
    }

    #[test]
    fn pattern_count_test() {
        let site = "AATT";
        let site_1 = "AATTGG";
        let pattern = get_pattern(&site);
        assert_eq!(2, pattern);
        assert_eq!(3, get_pattern(site_1));
    }

    #[test]
    fn pattern_count_all_test() {
        let site = "AAAA";
        let pattern = get_pattern(&site);
        assert_eq!(1, pattern);
    }

    #[test]
    fn count_parsimony_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT", "ATTA", "ATGC", "ATGA"];
        let mat = get_matrix(&id, &seq);
        let dna = index_sites(mat);
        assert_eq!(1, get_parsimony_informative(&dna));
    }
}
