//! Tree rooting operations: outgroup and midpoint rooting.

use ahash::AHashSet;
use anyhow::{Result, bail};

use crate::parser::tree::{Annotation, AnnotationTarget, NodeId, Rootedness, Tree};

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
        for &(v, len, _) in nbrs {
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

#[derive(Clone, Default)]
struct EdgeMeta {
    label: Option<String>,
    quoted: bool,
    node_annotations: Vec<Annotation>,
    branch_annotations: Vec<Annotation>,
}

impl EdgeMeta {
    fn is_empty(&self) -> bool {
        self.label.is_none()
            && self.node_annotations.is_empty()
            && self.branch_annotations.is_empty()
    }

    fn merge_with(self, other: Self) -> Self {
        if self.is_empty() {
            return other;
        }
        if other.is_empty() {
            return self;
        }
        let (label, quoted) = if self.label.is_some() {
            (self.label, self.quoted)
        } else {
            (other.label, other.quoted)
        };
        let node_annotations = if !self.node_annotations.is_empty() {
            self.node_annotations
        } else {
            other.node_annotations
        };
        let branch_annotations = if !self.branch_annotations.is_empty() {
            self.branch_annotations
        } else {
            other.branch_annotations
        };
        Self {
            label,
            quoted,
            node_annotations,
            branch_annotations,
        }
    }
}

struct TreeGraph {
    adj: Vec<Vec<(usize, Option<f64>, EdgeMeta)>>,
    labels: Vec<Option<String>>,
    quoted: Vec<bool>,
    leaf_node_annotations: Vec<Vec<Annotation>>,
    root_node_annotations: Vec<Annotation>,
    is_leaf: Vec<bool>,
    leaf_count: usize,
}

impl TreeGraph {
    fn from_tree(tree: &Tree) -> Self {
        let n = tree.len();
        let mut adj = vec![Vec::new(); n];
        let mut labels = vec![None; n];
        let mut quoted = vec![false; n];
        let mut leaf_node_annotations = vec![Vec::new(); n];
        let mut is_leaf = vec![false; n];
        let mut leaf_count = 0usize;
        let root_node_annotations = tree.node_annotations(tree.root()).to_vec();

        for id in tree.preorder() {
            let idx = id.index();
            labels[idx] = tree.label(id).map(str::to_string);
            quoted[idx] = tree.is_quoted_label(id);
            let leaf = tree.is_leaf(id);
            is_leaf[idx] = leaf;
            if leaf {
                leaf_count += 1;
                leaf_node_annotations[idx] = tree.node_annotations(id).to_vec();
            }
            if let Some(parent) = tree.parent(id) {
                let len = tree.branch_length(id);
                let edge_meta = if leaf {
                    EdgeMeta {
                        label: None,
                        quoted: false,
                        node_annotations: Vec::new(),
                        branch_annotations: tree.branch_annotations(id).to_vec(),
                    }
                } else {
                    EdgeMeta {
                        label: labels[idx].clone(),
                        quoted: quoted[idx],
                        node_annotations: tree.node_annotations(id).to_vec(),
                        branch_annotations: tree.branch_annotations(id).to_vec(),
                    }
                };
                adj[idx].push((parent.index(), len, edge_meta.clone()));
                adj[parent.index()].push((idx, len, edge_meta));
            }
        }

        // Collapse degree-2 root if present so the undirected graph has no redundant node.
        let root_idx = tree.root().index();
        if adj[root_idx].len() == 2 {
            let (c1, l1, m1) = adj[root_idx][0].clone();
            let (c2, l2, m2) = adj[root_idx][1].clone();
            let combined = match (l1, l2) {
                (Some(a), Some(b)) => Some(a + b),
                (Some(a), None) => Some(a),
                (None, Some(b)) => Some(b),
                (None, None) => None,
            };
            let combined_meta = m1.merge_with(m2);
            adj[root_idx].clear();
            for item in &mut adj[c1] {
                if item.0 == root_idx {
                    *item = (c2, combined, combined_meta.clone());
                }
            }
            for item in &mut adj[c2] {
                if item.0 == root_idx {
                    *item = (c1, combined, combined_meta.clone());
                }
            }
        }

        Self {
            adj,
            labels,
            quoted,
            leaf_node_annotations,
            root_node_annotations,
            is_leaf,
            leaf_count,
        }
    }

    fn edge_length(&self, u: usize, v: usize) -> Option<f64> {
        self.adj[u]
            .iter()
            .find(|(nbr, _, _)| *nbr == v)
            .and_then(|(_, len, _)| *len)
    }

    fn edge_meta(&self, u: usize, v: usize) -> Option<&EdgeMeta> {
        self.adj[u]
            .iter()
            .find(|(nbr, _, _)| *nbr == v)
            .map(|(_, _, meta)| meta)
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
        for &(v, _, _) in &graph.adj[u] {
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
            for &(v, _, _) in &graph.adj[u] {
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
        for &(v, len, _) in &graph.adj[u] {
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
    for annot in &graph.root_node_annotations {
        out.add_annotation(root, AnnotationTarget::Node, annot.clone());
    }

    let apply_node_and_edge = |out_tree: &mut Tree,
                               node_id: NodeId,
                               graph_node: usize,
                               edge_meta: Option<&EdgeMeta>,
                               bl: Option<f64>| {
        if graph.is_leaf[graph_node] {
            if let Some(lbl) = graph.labels[graph_node].as_deref() {
                out_tree.set_label(node_id, lbl, graph.quoted[graph_node]);
            }
            for annot in &graph.leaf_node_annotations[graph_node] {
                out_tree.add_annotation(node_id, AnnotationTarget::Node, annot.clone());
            }
            if let Some(meta) = edge_meta {
                for annot in &meta.branch_annotations {
                    out_tree.add_annotation(node_id, AnnotationTarget::Branch, annot.clone());
                }
            }
        } else if let Some(meta) = edge_meta {
            if let Some(lbl) = meta.label.as_deref() {
                out_tree.set_label(node_id, lbl, meta.quoted);
            }
            for annot in &meta.node_annotations {
                out_tree.add_annotation(node_id, AnnotationTarget::Node, annot.clone());
            }
            for annot in &meta.branch_annotations {
                out_tree.add_annotation(node_id, AnnotationTarget::Branch, annot.clone());
            }
        }
        if let Some(l) = bl {
            out_tree.set_branch_length(node_id, l);
        }
    };

    let split_edge_meta = graph.edge_meta(edge_u, edge_v);
    let mut mapped: Vec<Option<NodeId>> = vec![None; graph.adj.len()];
    let mut attach_branch =
        |start: usize, avoid: usize, bl: Option<f64>, prev: Option<NodeId>| -> NodeId {
            let start_id = out.add_node();
            mapped[start] = Some(start_id);
            apply_node_and_edge(&mut out, start_id, start, split_edge_meta, bl);
            out.attach_child(root, start_id, prev);

            let mut stack = vec![(start, avoid)];
            while let Some((u, p)) = stack.pop() {
                let parent_id = mapped[u].expect("Mapped graph node");
                let mut child_prev: Option<NodeId> = None;
                for (v, edge_len, edge_meta) in &graph.adj[u] {
                    if *v == p {
                        continue;
                    }
                    let child_id = out.add_node();
                    mapped[*v] = Some(child_id);
                    apply_node_and_edge(&mut out, child_id, *v, Some(edge_meta), *edge_len);
                    out.attach_child(parent_id, child_id, child_prev);
                    child_prev = Some(child_id);
                    stack.push((*v, u));
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

    #[test]
    fn test_reroot_preserves_internal_supports_and_annotations() {
        let tree =
            parse_tree_str("((A:0.1,B:0.2)95[&posterior=0.95]:0.3,C:0.4,(D:0.1,E:0.2)88:0.5);")
                .unwrap();

        let rooted_c = reroot_outgroup(&tree, &["C".into()]).unwrap();
        let out_nwk = crate::writer::tree::to_newick_string(&rooted_c, &Default::default());
        assert!(out_nwk.contains("95[&posterior=0.95]"));
        assert!(out_nwk.contains("88"));

        let mid = reroot_midpoint(&tree).unwrap();
        let mid_nwk = crate::writer::tree::to_newick_string(&mid, &Default::default());
        assert!(mid_nwk.contains("95[&posterior=0.95]"));
        assert!(mid_nwk.contains("88"));
    }
}
