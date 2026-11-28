//! Lay out trees for interactive viewers.
//!
//! The layout is computed once per tree in linear time. Viewers then ask for
//! the rows in view, so drawing costs depend on the screen size, not the
//! tree size. Each tip, or each folded clade, takes one row. Multi-tree
//! files are read one tree at a time.
//!
//! Values drawn on top of the branches, such as support or annotations, come
//! from [`TreeLayer`]s, so viewers can add new kinds of layers without
//! changing the layout.
//!
//! # Example
//!
//! ```
//! use segul_core::core::tree::view::{Segment, TreeScale, TreeView};
//! use segul_core::parser::tree::parse_tree_str;
//!
//! let tree = parse_tree_str("((A:1,B:2)95:1,C:3);").unwrap();
//! let mut view = TreeView::new(tree);
//! assert_eq!(view.row_count(), 3);
//! assert_eq!(view.label(view.row_node(0)), Some("A"));
//! // Tip B sits three units from the root.
//! assert_eq!(view.x(view.row_node(1), TreeScale::Length), 3.0);
//!
//! // Fold the clade of A and B into one row.
//! let clade = view.tree().parent(view.row_node(0)).unwrap();
//! view.toggle_fold(clade);
//! assert_eq!(view.row_count(), 2);
//! assert!(view
//!     .segments(0, 2, TreeScale::Length)
//!     .iter()
//!     .any(|s| matches!(s, Segment::Branch { row: 1, .. })));
//! ```

use std::fmt::Write as _;
use std::path::Path;
use std::sync::OnceLock;

use ahash::AHashSet;
use anyhow::{Context, Result};

use crate::core::tree::estimate::rooting::{reroot_midpoint_mapped, reroot_on_branch_mapped};
use crate::parser::tree::{
    AnnotationValue, NodeId, ParseOptions, Rootedness, Tree, TreeFmt, TreeReader,
};

const HIDDEN: u32 = u32::MAX;

/// How branch positions are measured.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TreeScale {
    /// Distance from the root along branch lengths.
    #[default]
    Length,
    /// Equal steps per node, with every tip aligned on the right.
    Cladogram,
}

/// A piece of the tree that crosses the rows in view.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Segment {
    /// Horizontal branch on `row` from the parent (`from`) to the node (`to`).
    /// `y` is the exact height in rows, midway between the node's first and
    /// last child; `row` is `y` rounded, for drawing on whole rows.
    Branch {
        node: NodeId,
        row: usize,
        y: f64,
        from: f64,
        to: f64,
    },
    /// Vertical line joining the children of `node`, from its first child's
    /// row (`top`) to its last child's row (`bottom`), or exactly from
    /// `top_y` to `bottom_y`. Rows may lie outside the requested range.
    Connector {
        node: NodeId,
        x: f64,
        top: usize,
        bottom: usize,
        top_y: f64,
        bottom_y: f64,
    },
}

/// A tree laid out in rows, with foldable clades.
pub struct TreeView {
    tree: Tree,
    has_lengths: bool,
    length_x: Vec<f64>,
    cladogram_x: Vec<f64>,
    length_range: (f64, f64),
    cladogram_max: f64,
    tips_below: Vec<u32>,
    folded: Vec<bool>,
    // Row to the tip or folded clade shown there.
    rows: Vec<NodeId>,
    // Node to its row, or `HIDDEN` inside a folded clade.
    node_rows: Vec<u32>,
    // Node to its exact height in rows; tips and folded clades sit on whole
    // rows, internal nodes midway between their first and last child.
    node_ys: Vec<f64>,
    // Lowercase tip names in tip order, built on the first search.
    search_names: OnceLock<Vec<(NodeId, String)>>,
    // The tree as read, kept once the view is rerooted.
    original: Option<Tree>,
}

impl TreeView {
    /// Lay out `tree`, with every clade unfolded.
    pub fn new(tree: Tree) -> Self {
        Self::with_original(tree, None)
    }

    fn with_original(tree: Tree, original: Option<Tree>) -> Self {
        let len = tree.len();
        let preorder: Vec<NodeId> = tree.preorder().collect();
        let has_lengths = preorder
            .iter()
            .any(|&id| !tree.is_root(id) && tree.branch_length(id).is_some());
        let mut length_x = vec![0.0; len];
        for &id in &preorder {
            if let Some(parent) = tree.parent(id) {
                let step = if has_lengths {
                    tree.branch_length(id).unwrap_or(0.0)
                } else {
                    1.0
                };
                length_x[id.index()] = length_x[parent.index()] + step;
            }
        }
        // Children come after parents in preorder, so walking it backwards
        // finishes every subtree before its parent.
        let mut heights = vec![0u32; len];
        let mut tips_below = vec![0u32; len];
        for &id in preorder.iter().rev() {
            if tree.is_leaf(id) {
                tips_below[id.index()] = 1;
            }
            if let Some(parent) = tree.parent(id) {
                let (p, i) = (parent.index(), id.index());
                heights[p] = heights[p].max(heights[i] + 1);
                tips_below[p] += tips_below[i];
            }
        }
        let cladogram_max = preorder
            .first()
            .map(|root| heights[root.index()] as f64)
            .unwrap_or(0.0);
        let cladogram_x = heights.iter().map(|&h| cladogram_max - h as f64).collect();
        let length_range = length_x
            .iter()
            .fold((0.0f64, 0.0f64), |(lo, hi), &x| (lo.min(x), hi.max(x)));
        let mut view = Self {
            tree,
            has_lengths,
            length_x,
            cladogram_x,
            length_range,
            cladogram_max,
            tips_below,
            folded: vec![false; len],
            rows: Vec::new(),
            node_rows: vec![HIDDEN; len],
            node_ys: vec![0.0; len],
            search_names: OnceLock::new(),
            original,
        };
        view.layout_rows();
        view
    }

    /// The tree being shown.
    pub fn tree(&self) -> &Tree {
        &self.tree
    }

    /// Label of `id`, if any.
    pub fn label(&self, id: NodeId) -> Option<&str> {
        self.tree.label(id)
    }

    /// Tips in the whole tree.
    pub fn tip_count(&self) -> usize {
        self.tree.tip_count()
    }

    /// Tips in the subtree of `id`.
    pub fn tips_below(&self, id: NodeId) -> usize {
        self.tips_below[id.index()] as usize
    }

    /// Whether any branch has a length. Without lengths every branch counts
    /// as one unit.
    pub fn has_lengths(&self) -> bool {
        self.has_lengths
    }

    /// Whether the file marks the tree as rooted or unrooted.
    pub fn rootedness(&self) -> Rootedness {
        self.tree.metadata.rootedness
    }

    /// Whether the view is rooted somewhere other than in the file.
    pub fn is_rerooted(&self) -> bool {
        self.original.is_some()
    }

    /// Root on the middle of the branch above `id`, keeping support values
    /// and annotations with their branches. Clades are unfolded. Returns the
    /// new id of `id`.
    pub fn root_on(&mut self, id: NodeId) -> Result<NodeId> {
        let (tree, mapped) = reroot_on_branch_mapped(&self.tree, id)?;
        self.replace_tree(tree);
        Ok(mapped
            .get(id.index())
            .copied()
            .flatten()
            .unwrap_or_else(|| self.tree.root()))
    }

    /// Root at the middle of the longest path between two tips. Needs branch
    /// lengths. Returns the new root.
    pub fn root_at_midpoint(&mut self) -> Result<NodeId> {
        let (tree, _) = reroot_midpoint_mapped(&self.tree)?;
        self.replace_tree(tree);
        Ok(self.tree.root())
    }

    /// Go back to the root from the file. Returns whether anything changed.
    pub fn restore_root(&mut self) -> bool {
        match self.original.take() {
            Some(original) => {
                *self = Self::new(original);
                true
            }
            None => false,
        }
    }

    // Lays out a rerooted tree, keeping the tree from the file.
    fn replace_tree(&mut self, tree: Tree) {
        let original = self
            .original
            .take()
            .unwrap_or_else(|| std::mem::take(&mut self.tree));
        *self = Self::with_original(tree, Some(original));
    }

    /// Rows in the current layout: one per tip or folded clade.
    pub fn row_count(&self) -> usize {
        self.rows.len()
    }

    /// Tip or folded clade on `row`.
    pub fn row_node(&self, row: usize) -> NodeId {
        self.rows[row]
    }

    /// Row of `id`, or `None` when it is inside a folded clade. Internal
    /// nodes sit midway between their first and last child.
    pub fn row_of(&self, id: NodeId) -> Option<usize> {
        match self.node_rows[id.index()] {
            HIDDEN => None,
            row => Some(row as usize),
        }
    }

    /// Exact height of `id` in rows, or `None` when it is inside a folded
    /// clade. Internal nodes sit midway between their first and last child,
    /// so they may fall between rows.
    pub fn y(&self, id: NodeId) -> Option<f64> {
        self.row_of(id).map(|_| self.node_ys[id.index()])
    }

    /// Horizontal position of `id` from the root.
    pub fn x(&self, id: NodeId, scale: TreeScale) -> f64 {
        match scale {
            TreeScale::Length => self.length_x[id.index()],
            TreeScale::Cladogram => self.cladogram_x[id.index()],
        }
    }

    /// Smallest and largest horizontal position. The smallest is below zero
    /// only with negative branch lengths.
    pub fn x_range(&self, scale: TreeScale) -> (f64, f64) {
        match scale {
            TreeScale::Length => self.length_range,
            TreeScale::Cladogram => (0.0, self.cladogram_max),
        }
    }

    /// Whether `id` is folded into one row.
    pub fn is_folded(&self, id: NodeId) -> bool {
        self.folded[id.index()]
    }

    /// Fold or unfold the clade of `id`. Tips cannot be folded. Returns
    /// whether anything changed.
    pub fn toggle_fold(&mut self, id: NodeId) -> bool {
        if self.tree.is_leaf(id) {
            return false;
        }
        self.folded[id.index()] = !self.folded[id.index()];
        self.layout_rows();
        true
    }

    /// Unfold every clade.
    pub fn unfold_all(&mut self) {
        self.folded.iter_mut().for_each(|f| *f = false);
        self.layout_rows();
    }

    /// Unfold the clades hiding `id`, so it gets a row.
    pub fn reveal(&mut self, id: NodeId) {
        let mut changed = false;
        let mut current = self.tree.parent(id);
        while let Some(node) = current {
            if self.folded[node.index()] {
                self.folded[node.index()] = false;
                changed = true;
            }
            current = self.tree.parent(node);
        }
        if changed {
            self.layout_rows();
        }
    }

    /// Branches and connectors that cross rows `first..first + count`, with
    /// positions on `scale`.
    ///
    /// Walks up from each row only until it reaches a node already seen, so
    /// the cost grows with the rows asked for plus the tree depth.
    pub fn segments(&self, first: usize, count: usize, scale: TreeScale) -> Vec<Segment> {
        let end = first.saturating_add(count).min(self.rows.len());
        let mut seen = AHashSet::new();
        let mut segments = Vec::new();
        for row in first..end {
            let mut current = Some(self.rows[row]);
            while let Some(id) = current {
                if !seen.insert(id) {
                    break;
                }
                let x = self.x(id, scale);
                let node_row = self.node_rows[id.index()] as usize;
                if (first..end).contains(&node_row) {
                    let from = self.tree.parent(id).map_or(x, |p| self.x(p, scale));
                    segments.push(Segment::Branch {
                        node: id,
                        row: node_row,
                        y: self.node_ys[id.index()],
                        from,
                        to: x,
                    });
                }
                if let Some((top, bottom)) = self.child_span(id)
                    && top < end
                    && bottom >= first
                {
                    let (top_y, bottom_y) = self.children_ys(id);
                    segments.push(Segment::Connector {
                        node: id,
                        x,
                        top,
                        bottom,
                        top_y,
                        bottom_y,
                    });
                }
                current = self.tree.parent(id);
            }
        }
        segments
    }

    /// Tips whose names contain `query`, ignoring case, in row order.
    pub fn find_tips(&self, query: &str) -> Vec<NodeId> {
        if query.is_empty() {
            return Vec::new();
        }
        let query = query.to_lowercase();
        self.search_names()
            .iter()
            .filter(|(_, name)| name.contains(&query))
            .map(|(id, _)| *id)
            .collect()
    }

    /// Internal nodes with a label that is not a number, such as clade names.
    pub fn named_clades(&self) -> Vec<NodeId> {
        self.tree
            .preorder()
            .filter(|&id| !self.tree.is_leaf(id))
            .filter(|&id| {
                self.tree
                    .label(id)
                    .is_some_and(|l| !l.trim().is_empty() && l.trim().parse::<f64>().is_err())
            })
            .collect()
    }

    /// Mean tip position for `bins` runs of rows, scaled to 0 to 1, for a
    /// whole-tree overview. Folded clades count at their own position.
    pub fn overview(&self, bins: usize, scale: TreeScale) -> Vec<f64> {
        let rows = self.rows.len();
        if bins == 0 || rows == 0 {
            return Vec::new();
        }
        let (lo, hi) = self.x_range(scale);
        let span = (hi - lo).max(f64::MIN_POSITIVE);
        (0..bins)
            .map(|bin| {
                let start = bin * rows / bins;
                let end = ((bin + 1) * rows / bins).max(start + 1).min(rows);
                if start >= rows {
                    return 0.0;
                }
                let sum: f64 = self.rows[start..end]
                    .iter()
                    .map(|&id| (self.x(id, scale) - lo) / span)
                    .sum();
                sum / (end - start) as f64
            })
            .collect()
    }

    /// The subtree of `id` as a Newick string with labels and branch lengths.
    pub fn subtree_newick(&self, id: NodeId) -> String {
        enum Step {
            Open(NodeId),
            Comma,
            Close(NodeId),
        }
        let mut out = String::new();
        let mut stack = vec![Step::Open(id)];
        while let Some(step) = stack.pop() {
            match step {
                Step::Open(node) if self.tree.is_leaf(node) => self.push_node(&mut out, node, id),
                Step::Open(node) => {
                    out.push('(');
                    stack.push(Step::Close(node));
                    let children: Vec<NodeId> = self.tree.children(node).collect();
                    for (i, &child) in children.iter().enumerate().rev() {
                        stack.push(Step::Open(child));
                        if i > 0 {
                            stack.push(Step::Comma);
                        }
                    }
                }
                Step::Comma => out.push(','),
                Step::Close(node) => {
                    out.push(')');
                    self.push_node(&mut out, node, id);
                }
            }
        }
        out.push(';');
        out
    }

    fn push_node(&self, out: &mut String, node: NodeId, top: NodeId) {
        if let Some(label) = self.tree.label(node) {
            out.push_str(&quote_label(label));
        }
        if node != top
            && let Some(length) = self.tree.branch_length(node)
        {
            let _ = write!(out, ":{length}");
        }
    }

    // First and last child rows of a visible, unfolded internal node.
    fn child_span(&self, id: NodeId) -> Option<(usize, usize)> {
        if self.node_rows[id.index()] == HIDDEN {
            return None;
        }
        self.children_rows(id)
    }

    // First and last child rows, before the node itself has a row.
    fn children_rows(&self, id: NodeId) -> Option<(usize, usize)> {
        if self.folded[id.index()] {
            return None;
        }
        let first = self.tree.first_child(id)?;
        let last = self.tree.children(id).last()?;
        Some((
            self.node_rows[first.index()] as usize,
            self.node_rows[last.index()] as usize,
        ))
    }

    // Exact heights of the first and last child of an unfolded node.
    fn children_ys(&self, id: NodeId) -> (f64, f64) {
        let first = self.tree.first_child(id).expect("Connector has children");
        let last = self.tree.children(id).last().unwrap_or(first);
        (self.node_ys[first.index()], self.node_ys[last.index()])
    }

    fn search_names(&self) -> &[(NodeId, String)] {
        self.search_names.get_or_init(|| {
            self.tree
                .leaves()
                .map(|id| (id, self.tree.label(id).unwrap_or("").to_lowercase()))
                .collect()
        })
    }

    // Gives each tip or folded clade a row, then places internal nodes
    // midway between their first and last child. Heights are kept exact, so
    // nested clades stay centered; rows are the heights rounded.
    fn layout_rows(&mut self) {
        self.rows.clear();
        self.node_rows.iter_mut().for_each(|r| *r = HIDDEN);
        if self.tree.is_empty() {
            return;
        }
        let mut visited = Vec::new();
        let mut stack = vec![self.tree.root()];
        while let Some(id) = stack.pop() {
            visited.push(id);
            if self.tree.is_leaf(id) || self.folded[id.index()] {
                self.node_rows[id.index()] = self.rows.len() as u32;
                self.node_ys[id.index()] = self.rows.len() as f64;
                self.rows.push(id);
                continue;
            }
            let start = stack.len();
            stack.extend(self.tree.children(id));
            stack[start..].reverse();
        }
        for &id in visited.iter().rev() {
            if self.children_rows(id).is_some() {
                let (top, bottom) = self.children_ys(id);
                let y = (top + bottom) / 2.0;
                self.node_ys[id.index()] = y;
                self.node_rows[id.index()] = y.round() as u32;
            }
        }
    }
}

// Newick needs quotes around labels with spaces or punctuation.
fn quote_label(label: &str) -> String {
    let needs_quotes = label
        .chars()
        .any(|c| c.is_whitespace() || "()[]{}:;,'\"".contains(c));
    if needs_quotes {
        format!("'{}'", label.replace('\'', "''"))
    } else {
        label.to_string()
    }
}

/// Read tree number `index` (zero-based) from a tree file, skipping the
/// trees before it without building them.
pub fn read_tree_at(path: &Path, fmt: TreeFmt, index: usize) -> Result<Tree> {
    let mut reader = TreeReader::open(path, fmt, ParseOptions::default())?;
    if index > 0 {
        reader.skip_trees(index)?;
    }
    reader
        .next()
        .transpose()?
        .with_context(|| format!("Failed finding tree {} in {}", index + 1, path.display()))
}

/// Count the trees in a tree file without building them.
pub fn count_trees(path: &Path, fmt: TreeFmt) -> Result<usize> {
    TreeReader::open(path, fmt, ParseOptions::default())?.skip_trees(usize::MAX)
}

/// A value a layer shows for a node.
#[derive(Debug, Clone, PartialEq)]
pub enum LayerValue {
    Number(f64),
    Text(String),
}

impl std::fmt::Display for LayerValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            // Whole numbers without decimals, others with up to four digits.
            Self::Number(v) if v.fract() == 0.0 && v.abs() < 1e15 => write!(f, "{v:.0}"),
            Self::Number(v) => {
                let text = format!("{v:.4}");
                f.write_str(text.trim_end_matches('0').trim_end_matches('.'))
            }
            Self::Text(s) => f.write_str(s),
        }
    }
}

/// Per-node values drawn over the tree, such as support or annotations.
///
/// Implement this to add a new layer; viewers decide how to draw the values
/// (as labels or branch colors).
pub trait TreeLayer: Send + Sync {
    /// Short name shown in the viewer.
    fn name(&self) -> &str;

    /// Value for `id`, if it has one.
    fn value(&self, tree: &Tree, id: NodeId) -> Option<LayerValue>;
}

/// Branch support from labels or annotations.
pub struct SupportLayer;

impl TreeLayer for SupportLayer {
    fn name(&self) -> &str {
        "Support"
    }

    fn value(&self, tree: &Tree, id: NodeId) -> Option<LayerValue> {
        if tree.is_leaf(id) {
            return None;
        }
        tree.support_value(id).map(LayerValue::Number)
    }
}

/// Length of the branch leading to each node.
pub struct BranchLengthLayer;

impl TreeLayer for BranchLengthLayer {
    fn name(&self) -> &str {
        "Lengths"
    }

    fn value(&self, tree: &Tree, id: NodeId) -> Option<LayerValue> {
        if tree.is_root(id) {
            return None;
        }
        tree.branch_length(id).map(LayerValue::Number)
    }
}

/// One key from BEAST, FigTree, MrBayes, or NHX annotations.
pub struct AnnotationLayer {
    key: String,
}

impl AnnotationLayer {
    pub fn new(key: &str) -> Self {
        Self {
            key: key.to_string(),
        }
    }
}

impl TreeLayer for AnnotationLayer {
    fn name(&self) -> &str {
        &self.key
    }

    fn value(&self, tree: &Tree, id: NodeId) -> Option<LayerValue> {
        tree.node_annotations(id)
            .iter()
            .chain(tree.branch_annotations(id))
            .flat_map(|a| a.entries.iter())
            .find(|(k, _)| *k == self.key)
            .map(|(_, v)| match v {
                AnnotationValue::Float(x) => LayerValue::Number(*x),
                AnnotationValue::Int(x) => LayerValue::Number(*x as f64),
                AnnotationValue::Str(s) => s
                    .trim()
                    .parse::<f64>()
                    .ok()
                    .filter(|x| x.is_finite())
                    .map_or_else(|| LayerValue::Text(s.clone()), LayerValue::Number),
                other => LayerValue::Text(other.to_string()),
            })
    }
}

/// Annotation keys used anywhere in `tree`, in the order first seen.
pub fn annotation_keys(tree: &Tree) -> Vec<String> {
    let mut seen = AHashSet::new();
    let mut keys = Vec::new();
    for id in tree.preorder() {
        let annotations = tree.node_annotations(id).iter();
        for annotation in annotations.chain(tree.branch_annotations(id)) {
            for (key, _) in &annotation.entries {
                if seen.insert(key.as_str()) {
                    keys.push(key.clone());
                }
            }
        }
    }
    keys
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::parser::tree::parse_tree_str;

    const TREES: &str = "../../tests/files/trees";

    fn view(newick: &str) -> TreeView {
        TreeView::new(parse_tree_str(newick).expect("Failed parsing tree"))
    }

    fn row_labels(view: &TreeView) -> Vec<String> {
        (0..view.row_count())
            .map(|r| view.label(view.row_node(r)).unwrap_or("").to_string())
            .collect()
    }

    #[test]
    fn test_layout() {
        let path = Path::new(TREES).join("simple.nwk");
        let tree = read_tree_at(&path, TreeFmt::Auto, 0).expect("Failed reading tree");
        let view = TreeView::new(tree);
        assert_eq!(view.tip_count(), 5);
        assert!(view.has_lengths());
        assert_eq!(
            row_labels(&view),
            ["Taxon_A", "Taxon_B", "Taxon_C", "Taxon_D", "Taxon_E"]
        );
        let b = view.row_node(1);
        assert!((view.x(b, TreeScale::Length) - 0.9).abs() < 1e-9);
        assert_eq!(view.x(b, TreeScale::Cladogram), 2.0);
        assert!((view.x_range(TreeScale::Length).1 - 0.9).abs() < 1e-9);
        // The clade of C, D, and E sits on the row of its middle child.
        let clade = view.tree().parent(view.row_node(2)).unwrap();
        assert_eq!(view.row_of(clade), Some(3));
        assert_eq!(view.tips_below(clade), 3);
        // The root sits midway between its clades at 0.5 and 3.
        let root = view.tree().root();
        assert_eq!(view.y(root), Some(1.75));
        assert_eq!(view.row_of(root), Some(2));
    }

    #[test]
    fn test_no_lengths() {
        let view = view("((A,B),C);");
        assert!(!view.has_lengths());
        assert_eq!(view.x(view.row_node(0), TreeScale::Length), 2.0);
        assert_eq!(view.x(view.row_node(2), TreeScale::Length), 1.0);
        assert_eq!(view.x(view.row_node(2), TreeScale::Cladogram), 2.0);
    }

    #[test]
    fn test_fold() {
        let mut view = view("((A,B)ab,(C,D,E)cde);");
        let ab = view.tree().parent(view.row_node(0)).unwrap();
        assert!(view.toggle_fold(ab));
        assert_eq!(row_labels(&view), ["ab", "C", "D", "E"]);
        assert!(view.is_folded(ab));
        assert_eq!(view.row_of(view.find_tips("a")[0]), None);
        assert!(!view.toggle_fold(view.row_node(1)));

        let d = view.find_tips("d")[0];
        let cde = view.tree().parent(d).unwrap();
        view.toggle_fold(cde);
        assert_eq!(view.row_count(), 2);
        view.reveal(d);
        assert_eq!(view.row_of(d), Some(2));
        assert!(view.is_folded(ab));
        view.unfold_all();
        assert_eq!(view.row_count(), 5);
    }

    #[test]
    fn test_segments_match_full_layout() {
        let view = view("(((A:1,B:1):1,C:2):1,((D:1,E:1):1,(F:1,(G:1,H:1):1):1):1);");
        let full = view.segments(0, view.row_count(), TreeScale::Length);
        // Every node has one branch.
        assert_eq!(
            full.iter()
                .filter(|s| matches!(s, Segment::Branch { .. }))
                .count(),
            view.tree().len()
        );
        // A slice holds exactly the full layout's pieces that cross it.
        for first in 0..view.row_count() {
            for count in 1..=3 {
                let end = first + count;
                let mut part = view.segments(first, count, TreeScale::Length);
                let mut expected: Vec<Segment> = full
                    .iter()
                    .copied()
                    .filter(|s| match *s {
                        Segment::Branch { row, .. } => row >= first && row < end,
                        Segment::Connector { top, bottom, .. } => top < end && bottom >= first,
                    })
                    .collect();
                let key = |s: &Segment| match *s {
                    Segment::Branch { node, .. } => (node, 0),
                    Segment::Connector { node, .. } => (node, 1),
                };
                part.sort_by_key(key);
                expected.sort_by_key(key);
                assert_eq!(part, expected, "rows {first}..{end}");
            }
        }
    }

    #[test]
    fn test_deep_tree() {
        let path = Path::new(TREES).join("caterpillar.nwk");
        let tree = read_tree_at(&path, TreeFmt::Auto, 0).expect("Failed reading tree");
        let view = TreeView::new(tree);
        assert_eq!(view.tip_count(), 2000);
        let segments = view.segments(1990, 10, TreeScale::Cladogram);
        assert!(segments.len() > 10);
        assert!(view.subtree_newick(view.tree().root()).ends_with(';'));
    }

    #[test]
    fn test_multi_tree_file() {
        let path = Path::new(TREES).join("multi.nwk");
        assert_eq!(count_trees(&path, TreeFmt::Auto).unwrap(), 3);
        let third = TreeView::new(read_tree_at(&path, TreeFmt::Auto, 2).unwrap());
        assert_eq!(row_labels(&third), ["TaxonE", "TaxonF"]);
        assert!(read_tree_at(&path, TreeFmt::Auto, 3).is_err());
    }

    #[test]
    fn test_find_and_clades() {
        let view = view("((Homo_sapiens,Pan)Hominini,Gorilla)95;");
        assert_eq!(view.find_tips("HOMO").len(), 1);
        assert!(view.find_tips("").is_empty());
        let clades = view.named_clades();
        assert_eq!(clades.len(), 1);
        assert_eq!(view.label(clades[0]), Some("Hominini"));
    }

    #[test]
    fn test_overview() {
        let view = view("((A:1,B:1):1,C:1);");
        let bins = view.overview(3, TreeScale::Length);
        assert_eq!(bins, vec![1.0, 1.0, 0.5]);
        assert!(view.overview(0, TreeScale::Length).is_empty());
    }

    #[test]
    fn test_subtree_newick() {
        let view = view("((A:1,'B c':2)90:0.5,C:3);");
        let clade = view.tree().parent(view.row_node(0)).unwrap();
        assert_eq!(view.subtree_newick(clade), "(A:1,'B c':2)90;");
    }

    #[test]
    fn test_layers() {
        let path = Path::new(TREES).join("beast_annotated.trees");
        let tree = read_tree_at(&path, TreeFmt::Auto, 0).expect("Failed reading tree");
        let keys = annotation_keys(&tree);
        assert!(keys.contains(&String::from("posterior")));
        assert!(keys.contains(&String::from("rate")));
        let view = TreeView::new(tree);
        let tree = view.tree();
        let clade = tree.parent(view.row_node(0)).unwrap();
        assert_eq!(
            SupportLayer.value(tree, clade),
            Some(LayerValue::Number(0.99))
        );
        assert_eq!(SupportLayer.value(tree, view.row_node(0)), None);
        let rate = AnnotationLayer::new("rate");
        assert_eq!(
            rate.value(tree, view.row_node(0)),
            Some(LayerValue::Number(0.002))
        );
        assert_eq!(
            BranchLengthLayer.value(tree, view.row_node(2)),
            Some(LayerValue::Number(0.5))
        );
        assert_eq!(LayerValue::Number(95.0).to_string(), "95");
        assert_eq!(LayerValue::Number(0.031_24).to_string(), "0.0312");
    }

    #[test]
    fn test_reroot() {
        let mut view = view("((A:1,B:2)90:0.5,C:1,(D:1,E:1)80:0.4);");
        assert!(!view.is_rerooted());
        let ab = view.tree().parent(view.row_node(0)).unwrap();
        view.toggle_fold(ab);
        let c = view.find_tips("c")[0];
        let new_c = view.root_on(c).unwrap();
        assert!(view.is_rerooted());
        assert_eq!(view.label(new_c), Some("C"));
        assert_eq!(view.tree().parent(new_c), Some(view.tree().root()));
        assert_eq!(view.row_count(), 5);
        assert_eq!(view.rootedness(), Rootedness::Rooted);
        // Rerooting twice still restores the file's root.
        view.root_on(view.find_tips("d")[0]).unwrap();
        assert!(view.restore_root());
        assert!(!view.is_rerooted());
        assert_eq!(row_labels(&view), ["A", "B", "C", "D", "E"]);
        assert!(!view.restore_root());
        let root = view.root_at_midpoint().unwrap();
        assert_eq!(root, view.tree().root());
        assert!(view.root_on(view.tree().root()).is_err());
        assert!(
            super::TreeView::new(parse_tree_str("((A,B),C);").unwrap())
                .root_at_midpoint()
                .is_err()
        );
    }
}
