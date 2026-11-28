//! Distance-based phylogenetic tree estimation: NJ, BIONJ, and UPGMA.

use anyhow::{Result, bail};

use crate::parser::tree::{NodeId, Rootedness, Tree};

use super::matrix::DistanceMatrix;

/// Estimate an unrooted tree using Saitou & Nei (1987) Neighbor-Joining.
pub fn estimate_nj(dm: &DistanceMatrix) -> Result<Tree> {
    estimate_nj_inner(dm, false)
}

/// Estimate an unrooted tree using Gascuel (1997) variance-weighted BIONJ.
pub fn estimate_bionj(dm: &DistanceMatrix) -> Result<Tree> {
    estimate_nj_inner(dm, true)
}

/// Estimate a rooted ultrametric tree using Sokal & Michener (1958) UPGMA.
pub fn estimate_upgma(dm: &DistanceMatrix) -> Result<Tree> {
    let n = dm.n;
    if n < 3 {
        bail!("At least 3 taxa are required for UPGMA (found {})", n);
    }

    let max_clusters = 2 * n - 1;
    let mut d = vec![0.0f64; max_clusters * max_clusters];
    for i in 0..n {
        for j in 0..n {
            d[i * max_clusters + j] = dm.get(i, j);
        }
    }

    let label_cap = dm.taxa.iter().map(String::len).sum::<usize>();
    let mut tree = Tree::with_capacity(max_clusters, label_cap);
    let mut cluster_nodes: Vec<Option<NodeId>> = vec![None; max_clusters];
    let mut sizes = vec![0usize; max_clusters];
    let mut heights = vec![0.0f64; max_clusters];
    let mut active: Vec<usize> = Vec::with_capacity(n);

    for (i, taxon) in dm.taxa.iter().enumerate() {
        let node = tree.add_node();
        tree.set_label(node, taxon, false);
        cluster_nodes[i] = Some(node);
        sizes[i] = 1;
        active.push(i);
    }

    let mut next_cluster = n;
    while active.len() > 1 {
        let mut best_ai = 0usize;
        let mut best_aj = 1usize;
        let mut best_dist = f64::INFINITY;

        for (ai, &ci) in active.iter().enumerate() {
            for (aj, &cj) in active.iter().enumerate().skip(ai + 1) {
                let dij = d[ci * max_clusters + cj];
                if dij < best_dist {
                    best_dist = dij;
                    best_ai = ai;
                    best_aj = aj;
                }
            }
        }

        let ci = active[best_ai];
        let cj = active[best_aj];
        let h_u = (0.5 * best_dist).max(heights[ci]).max(heights[cj]);
        let vi = (h_u - heights[ci]).max(0.0);
        let vj = (h_u - heights[cj]).max(0.0);

        let node_i = cluster_nodes[ci].expect("Missing UPGMA cluster node");
        let node_j = cluster_nodes[cj].expect("Missing UPGMA cluster node");
        let node_u = tree.add_node();
        tree.set_branch_length(node_i, vi);
        tree.set_branch_length(node_j, vj);
        tree.attach_child(node_u, node_i, None);
        tree.attach_child(node_u, node_j, Some(node_i));

        let cu = next_cluster;
        next_cluster += 1;
        cluster_nodes[cu] = Some(node_u);
        sizes[cu] = sizes[ci] + sizes[cj];
        heights[cu] = h_u;

        active.remove(best_aj);
        active.remove(best_ai);

        let si = sizes[ci] as f64;
        let sj = sizes[cj] as f64;
        let su = (sizes[ci] + sizes[cj]) as f64;
        for &ck in &active {
            let duk = (si * d[ci * max_clusters + ck] + sj * d[cj * max_clusters + ck]) / su;
            d[cu * max_clusters + ck] = duk;
            d[ck * max_clusters + cu] = duk;
        }
        active.push(cu);
    }

    let root = cluster_nodes[active[0]].expect("Missing UPGMA root node");
    tree.set_root(root);
    tree.metadata.rootedness = Rootedness::Rooted;
    Ok(tree)
}

fn estimate_nj_inner(dm: &DistanceMatrix, bionj: bool) -> Result<Tree> {
    let n = dm.n;
    if n < 3 {
        bail!("At least 3 taxa are required for NJ (found {})", n);
    }

    let max_clusters = 2 * n - 2;
    let mut d = vec![0.0f64; max_clusters * max_clusters];
    let mut v = if bionj {
        vec![0.0f64; max_clusters * max_clusters]
    } else {
        Vec::new()
    };

    for i in 0..n {
        for j in 0..n {
            let val = dm.get(i, j);
            d[i * max_clusters + j] = val;
            if bionj {
                v[i * max_clusters + j] = val;
            }
        }
    }

    let label_cap = dm.taxa.iter().map(String::len).sum::<usize>();
    let mut tree = Tree::with_capacity(max_clusters, label_cap);
    let mut cluster_nodes: Vec<Option<NodeId>> = vec![None; max_clusters];
    let mut active: Vec<usize> = Vec::with_capacity(n);

    for (i, taxon) in dm.taxa.iter().enumerate() {
        let node = tree.add_node();
        tree.set_label(node, taxon, false);
        cluster_nodes[i] = Some(node);
        active.push(i);
    }

    let mut row_sums = vec![0.0f64; n];
    let mut next_cluster = n;

    while active.len() > 3 {
        let m = active.len();
        let r = m as f64;

        row_sums[..m].fill(0.0);
        for ai in 0..m {
            let ci = active[ai];
            let mut sum = 0.0f64;
            for aj in 0..m {
                if ai != aj {
                    sum += d[ci * max_clusters + active[aj]];
                }
            }
            row_sums[ai] = sum;
        }

        let mut best_ai = 0usize;
        let mut best_aj = 1usize;
        let mut best_q = f64::INFINITY;

        for ai in 0..m {
            let ci = active[ai];
            let si = row_sums[ai];
            for aj in (ai + 1)..m {
                let cj = active[aj];
                let q = (r - 2.0) * d[ci * max_clusters + cj] - si - row_sums[aj];
                if q < best_q {
                    best_q = q;
                    best_ai = ai;
                    best_aj = aj;
                }
            }
        }

        let ci = active[best_ai];
        let cj = active[best_aj];
        let dij = d[ci * max_clusters + cj];
        let mut vi = 0.5 * dij + (row_sums[best_ai] - row_sums[best_aj]) / (2.0 * (r - 2.0));
        let mut vj = dij - vi;

        if vi < 0.0 {
            vj = (vj + vi).max(0.0);
            vi = 0.0;
        } else if vj < 0.0 {
            vi = (vi + vj).max(0.0);
            vj = 0.0;
        }

        let node_i = cluster_nodes[ci].expect("Missing NJ cluster node");
        let node_j = cluster_nodes[cj].expect("Missing NJ cluster node");
        let node_u = tree.add_node();
        tree.set_branch_length(node_i, vi);
        tree.set_branch_length(node_j, vj);
        tree.attach_child(node_u, node_i, None);
        tree.attach_child(node_u, node_j, Some(node_i));

        let cu = next_cluster;
        next_cluster += 1;
        cluster_nodes[cu] = Some(node_u);

        active.remove(best_aj);
        active.remove(best_ai);

        if bionj {
            let vij = v[ci * max_clusters + cj];
            let lambda = if vij <= 0.0 {
                0.5
            } else {
                let diff_sum: f64 = active
                    .iter()
                    .map(|&ck| v[cj * max_clusters + ck] - v[ci * max_clusters + ck])
                    .sum();
                (0.5 + diff_sum / (2.0 * (r - 2.0) * vij)).clamp(0.0, 1.0)
            };
            for &ck in &active {
                let dik = d[ci * max_clusters + ck];
                let djk = d[cj * max_clusters + ck];
                let duk = lambda * (dik - vi) + (1.0 - lambda) * (djk - vj);
                d[cu * max_clusters + ck] = duk;
                d[ck * max_clusters + cu] = duk;

                let vik = v[ci * max_clusters + ck];
                let vjk = v[cj * max_clusters + ck];
                let vuk = lambda * vik + (1.0 - lambda) * vjk - lambda * (1.0 - lambda) * vij;
                v[cu * max_clusters + ck] = vuk;
                v[ck * max_clusters + cu] = vuk;
            }
        } else {
            for &ck in &active {
                let duk = 0.5 * (d[ci * max_clusters + ck] + d[cj * max_clusters + ck] - dij);
                d[cu * max_clusters + ck] = duk;
                d[ck * max_clusters + cu] = duk;
            }
        }

        active.push(cu);
    }

    let a = active[0];
    let b = active[1];
    let c = active[2];
    let dab = d[a * max_clusters + b];
    let dac = d[a * max_clusters + c];
    let dbc = d[b * max_clusters + c];

    let va = (0.5 * (dab + dac - dbc)).max(0.0);
    let vb = (0.5 * (dab + dbc - dac)).max(0.0);
    let vc = (0.5 * (dac + dbc - dab)).max(0.0);

    let root = tree.add_node();
    tree.set_root(root);
    let na = cluster_nodes[a].expect("Missing NJ cluster node");
    let nb = cluster_nodes[b].expect("Missing NJ cluster node");
    let nc = cluster_nodes[c].expect("Missing NJ cluster node");
    tree.set_branch_length(na, va);
    tree.set_branch_length(nb, vb);
    tree.set_branch_length(nc, vc);
    tree.attach_child(root, na, None);
    tree.attach_child(root, nb, Some(na));
    tree.attach_child(root, nc, Some(nb));
    tree.metadata.rootedness = Rootedness::Unrooted;

    Ok(tree)
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use assert_approx_eq::assert_approx_eq;

    use super::*;

    fn tip_to_tip_distances(tree: &Tree) -> HashMap<(String, String), f64> {
        let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); tree.len()];
        for id in tree.preorder() {
            if let Some(parent) = tree.parent(id) {
                let len = tree.branch_length(id).unwrap_or(0.0);
                adj[id.index()].push((parent.index(), len));
                adj[parent.index()].push((id.index(), len));
            }
        }

        let leaves: Vec<(usize, String)> = tree
            .leaves()
            .map(|id| (id.index(), tree.label(id).unwrap().to_string()))
            .collect();

        let mut out = HashMap::new();
        for (src_idx, src_name) in &leaves {
            let mut dist = vec![f64::NAN; tree.len()];
            dist[*src_idx] = 0.0;
            let mut stack = vec![*src_idx];
            while let Some(u) = stack.pop() {
                for &(v, w) in &adj[u] {
                    if dist[v].is_nan() {
                        dist[v] = dist[u] + w;
                        stack.push(v);
                    }
                }
            }
            for (dst_idx, dst_name) in &leaves {
                out.insert((src_name.clone(), dst_name.clone()), dist[*dst_idx]);
            }
        }
        out
    }

    fn additive_5_taxon_matrix() -> DistanceMatrix {
        // Tree: ((A:0.1,B:0.2):0.3,C:0.4,(D:0.15,E:0.25):0.35)
        let taxa = ["A", "B", "C", "D", "E"]
            .iter()
            .map(|s| s.to_string())
            .collect::<Vec<_>>();
        let vals = vec![
            0.0, 0.3, 0.8, 0.9, 1.0, // A
            0.3, 0.0, 0.9, 1.0, 1.1, // B
            0.8, 0.9, 0.0, 0.9, 1.0, // C
            0.9, 1.0, 0.9, 0.0, 0.4, // D
            1.0, 1.1, 1.0, 0.4, 0.0, // E
        ];
        DistanceMatrix {
            taxa,
            n: 5,
            values: vals,
            undefined_pairs: 0,
        }
    }

    #[test]
    fn test_nj_and_bionj_exact_additive_recovery() {
        let dm = additive_5_taxon_matrix();
        for tree in [estimate_nj(&dm).unwrap(), estimate_bionj(&dm).unwrap()] {
            assert_eq!(tree.tip_count(), 5);
            assert_eq!(tree.metadata.rootedness, Rootedness::Unrooted);
            let recovered = tip_to_tip_distances(&tree);
            for i in 0..dm.n {
                for j in 0..dm.n {
                    let actual = recovered[&(dm.taxa[i].clone(), dm.taxa[j].clone())];
                    assert_approx_eq!(actual, dm.get(i, j), 1e-9);
                }
            }
        }
    }

    #[test]
    fn test_upgma_ultrametric_recovery() {
        // Ultrametric tree: (((A:1,B:1):1,C:2):2,(D:3,E:3):1) -> heights 1, 2, 3, 4
        let taxa = ["A", "B", "C", "D", "E"]
            .iter()
            .map(|s| s.to_string())
            .collect::<Vec<_>>();
        let vals = vec![
            0.0, 2.0, 4.0, 8.0, 8.0, // A
            2.0, 0.0, 4.0, 8.0, 8.0, // B
            4.0, 4.0, 0.0, 8.0, 8.0, // C
            8.0, 8.0, 8.0, 0.0, 6.0, // D
            8.0, 8.0, 8.0, 6.0, 0.0, // E
        ];
        let dm = DistanceMatrix {
            taxa,
            n: 5,
            values: vals,
            undefined_pairs: 0,
        };
        let tree = estimate_upgma(&dm).unwrap();
        assert_eq!(tree.tip_count(), 5);
        assert_eq!(tree.metadata.rootedness, Rootedness::Rooted);

        let recovered = tip_to_tip_distances(&tree);
        for i in 0..dm.n {
            for j in 0..dm.n {
                let actual = recovered[&(dm.taxa[i].clone(), dm.taxa[j].clone())];
                assert_approx_eq!(actual, dm.get(i, j), 1e-9);
            }
        }

        // Verify equal root-to-tip distance (height = 4.0) for every leaf
        for leaf in tree.leaves() {
            let mut height = 0.0f64;
            let mut cur = leaf;
            while let Some(p) = tree.parent(cur) {
                height += tree.branch_length(cur).unwrap();
                cur = p;
            }
            assert_approx_eq!(height, 4.0, 1e-9);
        }
    }
}
