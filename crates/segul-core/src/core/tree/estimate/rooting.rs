//! Tree rooting operations: outgroup and midpoint rooting.

use ahash::AHashSet;
use anyhow::{Result, bail};

use crate::parser::tree::{NodeId, Rootedness, Tree};

/// Rooting strategy applied after tree estimation.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub enum RootingMethod {
    /// Leave rooting as produced by the estimator (unrooted for NJ/BIONJ/MP, rooted for UPGMA).
    #[default]
    None,
    /// Root at the midpoint of the longest tip-to-tip path.
    Midpoint,
    /// Root on the branch separating the monophyletic outgroup taxa from the ingroup.
    Outgroup(Vec<String>),
}

/// Reroot `tree` on the edge separating `outgroup` from the remaining taxa.
pub fn reroot_outgroup(tree: &Tree, outgroup: &[String]) -> Result<Tree> {
    if outgroup.is_empty() {
        bail!("Outgroup taxon list must not be empty");
    }
    if tree.tip_count() < 3 {
        bail!("Cannot reroot a tree with fewer than 3 tips");
    }

    let graph = TreeGraph::from_tree(tree);
    let mut leaf_by_name = ahash::AHashMap::with_capacity(graph.leaf_count);
    for (idx, lbl) in graph.labels.iter().enumerate() {
        if graph.is_leaf[idx] {
            let name = lbl
                .as_deref()
                .ok_or_else(|| anyhow::anyhow!("Cannot reroot tree containing an unlabeled tip"))?;
            leaf_by_name.insert(name, idx);
        }
    }

    let mut out_set = AHashSet::with_capacity(outgroup.len());
    for name in outgroup {
        let trimmed = name.trim();
        let &idx = leaf_by_name.get(trimmed).ok_or_else(|| {
            anyhow::anyhow!("Outgroup taxon '{}' was not found in the tree", trimmed)
        })?;
        out_set.insert(idx);
    }

    if out_set.len() >= graph.leaf_count {
        bail!(
            "Outgroup cannot include all {} taxa in the tree",
            graph.leaf_count
        );
    }

    let target_edge = find_bipartition_edge(&graph, &out_set).ok_or_else(|| {
        anyhow::anyhow!(
            "Specified outgroup taxa ({}) do not form a monophyletic clade in the tree",
            outgroup.join(", ")
        )
    })?;

    let (u, v, len) = target_edge;
    let half = len.map(|l| l * 0.5);
    Ok(rebuild_rooted_tree(&graph, tree, u, v, half, half))
}

/// Reroot `tree` at the midpoint of its longest tip-to-tip path (tree diameter).
pub fn reroot_midpoint(tree: &Tree) -> Result<Tree> {
    if tree.tip_count() < 2 {
        bail!("Cannot midpoint-root a tree with fewer than 2 tips");
    }

    let graph = TreeGraph::from_tree(tree);
    for (u, nbrs) in graph.adj.iter().enumerate() {
        for &(v, len) in nbrs {
            if u < v {
                let l = len.ok_or_else(|| {
                    anyhow::anyhow!(
                        "Midpoint rooting requires branch lengths on all edges \
                        (missing branch length on tree)"
                    )
                })?;
                if !l.is_finite() || l < 0.0 {
                    bail!("Midpoint rooting requires finite non-negative branch lengths");
                }
            }
        }
    }

    let first_leaf = graph
        .is_leaf
        .iter()
        .position(|&is_l| is_l)
        .ok_or_else(|| anyhow::anyhow!("Tree has no leaf nodes"))?;

    let (end_a, _, _) = farthest_leaf(&graph, first_leaf);
    let (end_b, diameter, parents) = farthest_leaf(&graph, end_a);
    if diameter <= 0.0 {
        bail!("Cannot midpoint-root a tree with zero total branch length");
    }

    // Reconstruct path from end_a to end_b
    let mut path = Vec::new();
    let mut cur = end_b;
    while cur != usize::MAX {
        path.push(cur);
        cur = parents[cur];
    }
    path.reverse();

    let target = diameter * 0.5;
    let mut acc = 0.0f64;
    for (idx, win) in path.windows(2).enumerate() {
        let u = win[0];
        let v = win[1];
        let edge_len = graph.edge_length(u, v).unwrap_or(0.0);
        if acc + edge_len >= target || idx + 2 == path.len() {
            let d_u = (target - acc).clamp(0.0, edge_len);
            let d_v = (edge_len - d_u).max(0.0);
            return Ok(rebuild_rooted_tree(
                &graph,
                tree,
                u,
                v,
                Some(d_u),
                Some(d_v),
            ));
        }
        acc += edge_len;
    }

    unreachable!("Midpoint must lie along diameter path")
}

struct TreeGraph {
    adj: Vec<Vec<(usize, Option<f64>)>>,
    labels: Vec<Option<String>>,
    quoted: Vec<bool>,
    is_leaf: Vec<bool>,
    leaf_count: usize,
}

impl TreeGraph {
    fn from_tree(tree: &Tree) -> Self {
        let n = tree.len();
        let mut adj = vec![Vec::new(); n];
        let mut labels = vec![None; n];
        let mut quoted = vec![false; n];
        let mut is_leaf = vec![false; n];
        let mut leaf_count = 0usize;

        for id in tree.preorder() {
            let idx = id.index();
            labels[idx] = tree.label(id).map(str::to_string);
            quoted[idx] = tree.is_quoted_label(id);
            let leaf = tree.is_leaf(id);
            is_leaf[idx] = leaf;
            if leaf {
                leaf_count += 1;
            }
            if let Some(parent) = tree.parent(id) {
                let len = tree.branch_length(id);
                adj[idx].push((parent.index(), len));
                adj[parent.index()].push((idx, len));
            }
        }

        // Collapse degree-2 root if present so the undirected graph has no redundant node.
        let root_idx = tree.root().index();
        if adj[root_idx].len() == 2 {
            let (c1, l1) = adj[root_idx][0];
            let (c2, l2) = adj[root_idx][1];
            let combined = match (l1, l2) {
                (Some(a), Some(b)) => Some(a + b),
                (Some(a), None) => Some(a),
                (None, Some(b)) => Some(b),
                (None, None) => None,
            };
            adj[root_idx].clear();
            for item in &mut adj[c1] {
                if item.0 == root_idx {
                    *item = (c2, combined);
                }
            }
            for item in &mut adj[c2] {
                if item.0 == root_idx {
                    *item = (c1, combined);
                }
            }
        }

        Self {
            adj,
            labels,
            quoted,
            is_leaf,
            leaf_count,
        }
    }

    fn edge_length(&self, u: usize, v: usize) -> Option<f64> {
        self.adj[u]
            .iter()
            .find(|&&(nbr, _)| nbr == v)
            .and_then(|&(_, len)| len)
    }
}

fn find_bipartition_edge(
    graph: &TreeGraph,
    out_set: &AHashSet<usize>,
) -> Option<(usize, usize, Option<f64>)> {
    let start = graph.is_leaf.iter().position(|&l| l)?;
    let target_size = out_set.len();
    let total_leaves = graph.leaf_count;

    // Postorder traversal collecting (total_leaves_in_subtree, outgroup_leaves_in_subtree).
    let n = graph.adj.len();
    let mut parent = vec![usize::MAX; n];
    let mut order = Vec::with_capacity(n);
    let mut stack = vec![(start, usize::MAX)];
    while let Some((u, p)) = stack.pop() {
        parent[u] = p;
        order.push(u);
        for &(v, _) in &graph.adj[u] {
            if v != p {
                stack.push((v, u));
            }
        }
    }

    let mut sub_total = vec![0usize; n];
    let mut sub_out = vec![0usize; n];

    for &u in order.iter().rev() {
        if graph.is_leaf[u] {
            sub_total[u] = 1;
            sub_out[u] = usize::from(out_set.contains(&u));
        } else {
            let mut tot = 0usize;
            let mut out = 0usize;
            for &(v, _) in &graph.adj[u] {
                if v != parent[u] {
                    tot += sub_total[v];
                    out += sub_out[v];
                }
            }
            sub_total[u] = tot;
            sub_out[u] = out;
        }

        let p = parent[u];
        if p != usize::MAX {
            let is_exact_out = sub_out[u] == target_size && sub_total[u] == target_size;
            let is_exact_in = sub_out[u] == 0 && (total_leaves - sub_total[u]) == target_size;
            if is_exact_out || is_exact_in {
                let len = graph.edge_length(u, p);
                return Some((u, p, len));
            }
        }
    }
    None
}

fn farthest_leaf(graph: &TreeGraph, start: usize) -> (usize, f64, Vec<usize>) {
    let n = graph.adj.len();
    let mut dist = vec![f64::NAN; n];
    let mut parent = vec![usize::MAX; n];
    dist[start] = 0.0;

    let mut stack = vec![start];
    let mut best_leaf = start;
    let mut best_dist = -1.0f64;

    while let Some(u) = stack.pop() {
        if graph.is_leaf[u] && dist[u] > best_dist {
            best_dist = dist[u];
            best_leaf = u;
        }
        for &(v, len) in &graph.adj[u] {
            if dist[v].is_nan() {
                dist[v] = dist[u] + len.unwrap_or(0.0);
                parent[v] = u;
                stack.push(v);
            }
        }
    }

    (best_leaf, best_dist, parent)
}

fn rebuild_rooted_tree(
    graph: &TreeGraph,
    orig: &Tree,
    edge_u: usize,
    edge_v: usize,
    len_u: Option<f64>,
    len_v: Option<f64>,
) -> Tree {
    let label_cap: usize = graph.labels.iter().flatten().map(String::len).sum();
    let mut out = Tree::with_capacity(graph.adj.len() + 1, label_cap);
    out.metadata = orig.metadata.clone();
    out.metadata.rootedness = Rootedness::Rooted;

    let root = out.add_node();
    out.set_root(root);

    let mut mapped: Vec<Option<NodeId>> = vec![None; graph.adj.len()];
    let mut attach_branch =
        |start: usize, avoid: usize, bl: Option<f64>, prev: Option<NodeId>| -> NodeId {
            let start_id = out.add_node();
            mapped[start] = Some(start_id);
            if let Some(lbl) = graph.labels[start].as_deref()
                && graph.is_leaf[start]
            {
                out.set_label(start_id, lbl, graph.quoted[start]);
            }
            if let Some(l) = bl {
                out.set_branch_length(start_id, l);
            }
            out.attach_child(root, start_id, prev);

            let mut stack = vec![(start, avoid)];
            while let Some((u, p)) = stack.pop() {
                let parent_id = mapped[u].expect("Mapped graph node");
                let mut child_prev: Option<NodeId> = None;
                for &(v, edge_len) in &graph.adj[u] {
                    if v == p {
                        continue;
                    }
                    let child_id = out.add_node();
                    mapped[v] = Some(child_id);
                    if let Some(lbl) = graph.labels[v].as_deref()
                        && graph.is_leaf[v]
                    {
                        out.set_label(child_id, lbl, graph.quoted[v]);
                    }
                    if let Some(l) = edge_len {
                        out.set_branch_length(child_id, l);
                    }
                    out.attach_child(parent_id, child_id, child_prev);
                    child_prev = Some(child_id);
                    stack.push((v, u));
                }
            }
            start_id
        };

    let child_u = attach_branch(edge_u, edge_v, len_u, None);
    let _child_v = attach_branch(edge_v, edge_u, len_v, Some(child_u));

    out
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;

    use crate::parser::tree::parse_tree_str;

    use super::*;

    #[test]
    fn test_outgroup_rooting_single_and_clade() {
        let tree = parse_tree_str("((A:0.1,B:0.2):0.3,C:0.4,(D:0.1,E:0.2):0.5);").unwrap();

        let rooted_c = reroot_outgroup(&tree, &["C".into()]).unwrap();
        assert_eq!(rooted_c.metadata.rootedness, Rootedness::Rooted);
        assert_eq!(rooted_c.tip_count(), 5);
        assert_eq!(rooted_c.child_count(rooted_c.root()), 2);

        let rooted_de = reroot_outgroup(&tree, &["D".into(), "E".into()]).unwrap();
        assert_eq!(rooted_de.metadata.rootedness, Rootedness::Rooted);
        assert_eq!(rooted_de.tip_count(), 5);

        // Non-monophyletic outgroup must error
        let err = reroot_outgroup(&tree, &["A".into(), "D".into()]).unwrap_err();
        assert!(err.to_string().contains("monophyletic"));

        // Unknown outgroup must error
        let err_missing = reroot_outgroup(&tree, &["UNKNOWN".into()]).unwrap_err();
        assert!(err_missing.to_string().contains("not found"));
    }

    #[test]
    fn test_midpoint_rooting_balances_diameter() {
        let tree = parse_tree_str("((A:1.0,B:1.0):1.0,C:2.0,(D:4.0,E:1.0):2.0);").unwrap();
        // Diameter is A/B(1+1) + (2+4)D = 8.0, or C(2) + (2+4)D = 8.0. Max dist from root = 4.0!
        let rooted = reroot_midpoint(&tree).unwrap();
        assert_eq!(rooted.metadata.rootedness, Rootedness::Rooted);
        assert_eq!(rooted.tip_count(), 5);

        let root_children: Vec<NodeId> = rooted.children(rooted.root()).collect();
        assert_eq!(root_children.len(), 2);

        fn max_depth(t: &Tree, id: NodeId) -> f64 {
            let bl = t.branch_length(id).unwrap_or(0.0);
            if t.is_leaf(id) {
                bl
            } else {
                bl + t
                    .children(id)
                    .map(|c| max_depth(t, c))
                    .fold(0.0f64, f64::max)
            }
        }

        let left = max_depth(&rooted, root_children[0]);
        let right = max_depth(&rooted, root_children[1]);
        assert_approx_eq!(left, 4.0, 1e-9);
        assert_approx_eq!(right, 4.0, 1e-9);
    }
}
