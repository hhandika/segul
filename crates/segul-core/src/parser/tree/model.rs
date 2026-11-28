//! Cache-compact arena-backed phylogenetic tree and network data model.

use std::fmt;
use std::sync::Arc;

use ahash::AHashMap;
use indexmap::IndexSet;

const NONE: u32 = u32::MAX;

const FLAG_INTERNED_LABEL: u32 = 1 << 0;
const FLAG_QUOTED_LABEL: u32 = 1 << 1;
const FLAG_HAS_NODE_ANNOT: u32 = 1 << 2;
const FLAG_HAS_BRANCH_ANNOT: u32 = 1 << 3;
const FLAG_IS_HYBRID: u32 = 1 << 4;

/// Compact 32-bit node identifier into a [`Tree`] arena.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct NodeId(u32);

impl NodeId {
    /// Return the zero-based arena index of this node.
    #[inline]
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

/// Compact 32-bit identifier into a shared [`TaxonTable`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct TaxonId(u32);

impl TaxonId {
    /// Return the zero-based index of this taxon.
    #[inline]
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

/// 32-byte Left-Child Right-Sibling (LCRS) node stored in a flat arena.
///
/// Two nodes fit in a single 64-byte cache line with zero per-node heap
/// allocations.
#[derive(Debug, Clone, Copy)]
pub struct Node {
    parent: u32,
    first_child: u32,
    next_sibling: u32,
    label_start: u32,
    label_len: u32,
    flags: u32,
    branch_length: f64,
}

impl Node {
    #[inline]
    fn empty() -> Self {
        Self {
            parent: NONE,
            first_child: NONE,
            next_sibling: NONE,
            label_start: 0,
            label_len: 0,
            flags: 0,
            branch_length: f64::NAN,
        }
    }

    /// Return true if this node has no children.
    #[inline]
    pub fn is_leaf(&self) -> bool {
        self.first_child == NONE
    }

    /// Return true if this node has no parent.
    #[inline]
    pub fn is_root(&self) -> bool {
        self.parent == NONE
    }

    /// Return the branch length leading into this node, if present.
    #[inline]
    pub fn branch_length(&self) -> Option<f64> {
        if self.branch_length.is_nan() {
            None
        } else {
            Some(self.branch_length)
        }
    }
}

/// Ordered deduplicated table of taxon names shared across multi-tree files.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TaxonTable {
    taxa: IndexSet<String>,
}

impl TaxonTable {
    /// Create an empty taxon table.
    pub fn new() -> Self {
        Self {
            taxa: IndexSet::new(),
        }
    }

    /// Create an empty taxon table pre-allocated for `capacity` taxa.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            taxa: IndexSet::with_capacity(capacity),
        }
    }

    /// Insert `name` if not already present and return its [`TaxonId`].
    pub fn intern(&mut self, name: &str) -> TaxonId {
        if let Some(idx) = self.taxa.get_index_of(name) {
            TaxonId(u32::try_from(idx).expect("TaxonTable index overflow"))
        } else {
            let (idx, _) = self.taxa.insert_full(name.to_string());
            TaxonId(u32::try_from(idx).expect("TaxonTable index overflow"))
        }
    }

    /// Look up the taxon name for `id`.
    #[inline]
    pub fn get(&self, id: TaxonId) -> Option<&str> {
        self.taxa.get_index(id.index()).map(String::as_str)
    }

    /// Look up the [`TaxonId`] for `name`, if present.
    #[inline]
    pub fn get_index_of(&self, name: &str) -> Option<TaxonId> {
        self.taxa.get_index_of(name).map(|idx| TaxonId(idx as u32))
    }

    /// Number of unique taxa in the table.
    #[inline]
    pub fn len(&self) -> usize {
        self.taxa.len()
    }

    /// Return true when the table contains no taxa.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.taxa.is_empty()
    }

    /// Iterate over taxon names in insertion order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = &str> {
        self.taxa.iter().map(String::as_str)
    }
}

/// Target to which a bracket annotation is attached.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AnnotationTarget {
    /// Comment placed before `:` (attached to the node itself).
    Node,
    /// Comment placed after `:` (attached to the subtending branch/edge).
    Branch,
}

/// Typed value parsed from BEAST, FigTree, MrBayes, or NHX metadata comments.
#[derive(Debug, Clone, PartialEq)]
pub enum AnnotationValue {
    Str(String),
    Float(f64),
    Int(i64),
    Bool(bool),
    Array(Vec<AnnotationValue>),
}

impl fmt::Display for AnnotationValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Str(s) => f.write_str(s),
            Self::Float(v) => write!(f, "{v}"),
            Self::Int(v) => write!(f, "{v}"),
            Self::Bool(v) => write!(f, "{v}"),
            Self::Array(items) => {
                f.write_str("{")?;
                for (i, item) in items.iter().enumerate() {
                    if i > 0 {
                        f.write_str(",")?;
                    }
                    write!(f, "{item}")?;
                }
                f.write_str("}")
            }
        }
    }
}

/// Dialect of a parsed bracket annotation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AnnotationFormat {
    /// BEAST / FigTree / MrBayes `[&key=value,...]` metadata.
    Beast,
    /// New Hampshire Extended `[&&NHX:key=value:...]` metadata.
    Nhx,
    /// Unstructured comment or branch support such as `[95]`.
    Raw,
}

/// Parsed bracket comment attached to a node or branch.
#[derive(Debug, Clone, PartialEq)]
pub struct Annotation {
    pub format: AnnotationFormat,
    pub entries: Vec<(String, AnnotationValue)>,
    pub raw: String,
}

/// Extended Newick (eNewick) reticulation metadata for a hybrid node.
#[derive(Debug, Clone, PartialEq)]
pub struct HybridInfo {
    /// Hybrid identifier tag including `#`, e.g. `"#H1"`.
    pub tag: String,
    /// Optional second colon field (`:len:support:gamma`).
    pub edge_support: Option<f64>,
    /// Optional inheritance probability (`gamma`) in the third colon field.
    pub inheritance_prob: Option<f64>,
}

/// Sparse side-table storing node/branch annotations and eNewick hybrid metadata.
#[derive(Debug, Clone, Default)]
pub struct AnnotationStore {
    node_annotations: AHashMap<NodeId, Vec<Annotation>>,
    branch_annotations: AHashMap<NodeId, Vec<Annotation>>,
    hybrids: AHashMap<NodeId, HybridInfo>,
}

/// Rootedness state declared by a tree command or `[&R]` / `[&U]` tag.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Rootedness {
    #[default]
    Unknown,
    Rooted,
    Unrooted,
}

/// Tree-level metadata parsed from NEXUS headers or leading Newick comments.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TreeMetadata {
    pub name: Option<String>,
    pub rootedness: Rootedness,
    pub weight: Option<f64>,
    pub is_default: bool,
    pub comments: Vec<String>,
}

/// Arena-backed phylogenetic tree or Extended Newick network.
#[derive(Debug, Clone)]
pub struct Tree {
    nodes: Vec<Node>,
    root: NodeId,
    labels: String,
    taxa: Option<Arc<TaxonTable>>,
    annotations: Option<Box<AnnotationStore>>,
    pub metadata: TreeMetadata,
    leaf_count: usize,
}

impl Default for Tree {
    fn default() -> Self {
        Self::new()
    }
}

impl Tree {
    /// Create an empty tree with no nodes.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            root: NodeId(0),
            labels: String::new(),
            taxa: None,
            annotations: None,
            metadata: TreeMetadata::default(),
            leaf_count: 0,
        }
    }

    /// Create an empty tree pre-allocated for `node_cap` nodes and `label_cap` label bytes.
    pub fn with_capacity(node_cap: usize, label_cap: usize) -> Self {
        Self {
            nodes: Vec::with_capacity(node_cap),
            root: NodeId(0),
            labels: String::with_capacity(label_cap),
            taxa: None,
            annotations: None,
            metadata: TreeMetadata::default(),
            leaf_count: 0,
        }
    }

    /// Set the shared [`TaxonTable`] used for interned leaf labels.
    pub fn set_taxa(&mut self, taxa: Arc<TaxonTable>) {
        self.taxa = Some(taxa);
    }

    /// Return the shared [`TaxonTable`], if this tree uses interned taxa.
    #[inline]
    pub fn taxa(&self) -> Option<&Arc<TaxonTable>> {
        self.taxa.as_ref()
    }

    /// Allocate a new disconnected node in the arena and return its [`NodeId`].
    #[inline]
    pub(crate) fn add_node(&mut self) -> NodeId {
        let raw = u32::try_from(self.nodes.len()).expect("Tree node count exceeds u32::MAX");
        assert!(raw != NONE, "Tree node count reached sentinel limit");
        self.nodes.push(Node::empty());
        self.leaf_count += 1;
        NodeId(raw)
    }

    /// Set the root node of the tree.
    #[inline]
    pub(crate) fn set_root(&mut self, root: NodeId) {
        self.root = root;
    }

    /// Attach `child` as the last child of `parent` in O(1) when `prev_sibling` is known.
    #[inline]
    pub(crate) fn attach_child(
        &mut self,
        parent: NodeId,
        child: NodeId,
        prev_sibling: Option<NodeId>,
    ) {
        self.nodes[child.index()].parent = parent.0;
        if let Some(prev) = prev_sibling {
            self.nodes[prev.index()].next_sibling = child.0;
            return;
        }
        let parent_first = self.nodes[parent.index()].first_child;
        if parent_first == NONE {
            self.nodes[parent.index()].first_child = child.0;
            self.leaf_count -= 1;
        } else {
            let mut cur = parent_first;
            while self.nodes[cur as usize].next_sibling != NONE {
                cur = self.nodes[cur as usize].next_sibling;
            }
            self.nodes[cur as usize].next_sibling = child.0;
        }
    }

    /// Store an inline label in the tree's contiguous label buffer.
    #[inline]
    pub(crate) fn set_label(&mut self, id: NodeId, label: &str, quoted: bool) {
        let start = u32::try_from(self.labels.len()).expect("Label buffer exceeds u32::MAX");
        let len = u32::try_from(label.len()).expect("Label length exceeds u32::MAX");
        self.labels.push_str(label);
        let node = &mut self.nodes[id.index()];
        node.label_start = start;
        node.label_len = len;
        node.flags &= !FLAG_INTERNED_LABEL;
        if quoted {
            node.flags |= FLAG_QUOTED_LABEL;
        }
    }

    /// Store a reference to an interned taxon in the shared [`TaxonTable`].
    #[inline]
    pub(crate) fn set_interned_label(&mut self, id: NodeId, taxon_id: TaxonId, quoted: bool) {
        let node = &mut self.nodes[id.index()];
        node.label_start = taxon_id.0;
        node.label_len = 0;
        node.flags |= FLAG_INTERNED_LABEL;
        if quoted {
            node.flags |= FLAG_QUOTED_LABEL;
        }
    }

    /// Set the branch length leading into `id`.
    #[inline]
    pub(crate) fn set_branch_length(&mut self, id: NodeId, length: f64) {
        self.nodes[id.index()].branch_length = length;
    }

    /// Attach a parsed annotation to a node or its subtending branch.
    pub(crate) fn add_annotation(
        &mut self,
        id: NodeId,
        target: AnnotationTarget,
        annotation: Annotation,
    ) {
        let store = self
            .annotations
            .get_or_insert_with(|| Box::new(AnnotationStore::default()));
        match target {
            AnnotationTarget::Node => {
                self.nodes[id.index()].flags |= FLAG_HAS_NODE_ANNOT;
                store
                    .node_annotations
                    .entry(id)
                    .or_default()
                    .push(annotation);
            }
            AnnotationTarget::Branch => {
                self.nodes[id.index()].flags |= FLAG_HAS_BRANCH_ANNOT;
                store
                    .branch_annotations
                    .entry(id)
                    .or_default()
                    .push(annotation);
            }
        }
    }

    /// Record Extended Newick hybrid metadata for `id`.
    pub(crate) fn set_hybrid_info(&mut self, id: NodeId, info: HybridInfo) {
        self.nodes[id.index()].flags |= FLAG_IS_HYBRID;
        let store = self
            .annotations
            .get_or_insert_with(|| Box::new(AnnotationStore::default()));
        store.hybrids.insert(id, info);
    }

    /// Return a mutable reference to `id`'s [`HybridInfo`] if present.
    pub(crate) fn hybrid_info_mut(&mut self, id: NodeId) -> Option<&mut HybridInfo> {
        self.annotations.as_mut()?.hybrids.get_mut(&id)
    }

    /// Root node identifier.
    #[inline]
    pub fn root(&self) -> NodeId {
        self.root
    }

    /// Node at arena index `index`, or `None` past the last node. Use it to
    /// turn an index from [`NodeId::index`] back into a [`NodeId`].
    #[inline]
    pub fn node_id(&self, index: usize) -> Option<NodeId> {
        (index < self.nodes.len()).then_some(NodeId(index as u32))
    }

    /// Total number of nodes (tips + internal nodes) in the tree.
    #[inline]
    pub fn len(&self) -> usize {
        self.nodes.len()
    }

    /// Return true when the tree has no nodes.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty()
    }

    /// Number of leaf (tip) nodes in the tree.
    #[inline]
    pub fn tip_count(&self) -> usize {
        self.leaf_count
    }

    /// Return the parent of `id`, or `None` for the root.
    #[inline]
    pub fn parent(&self, id: NodeId) -> Option<NodeId> {
        let raw = self.nodes.get(id.index())?.parent;
        if raw == NONE { None } else { Some(NodeId(raw)) }
    }

    /// Return the first child of `id`, or `None` for a leaf.
    #[inline]
    pub fn first_child(&self, id: NodeId) -> Option<NodeId> {
        let raw = self.nodes.get(id.index())?.first_child;
        if raw == NONE { None } else { Some(NodeId(raw)) }
    }

    /// Return the next right sibling of `id`, or `None` if `id` is the last child.
    #[inline]
    pub fn next_sibling(&self, id: NodeId) -> Option<NodeId> {
        let raw = self.nodes.get(id.index())?.next_sibling;
        if raw == NONE { None } else { Some(NodeId(raw)) }
    }

    /// Iterate over direct children of `id` in left-to-right order.
    #[inline]
    pub fn children(&self, id: NodeId) -> ChildrenIter<'_> {
        ChildrenIter {
            tree: self,
            next: self.first_child(id),
        }
    }

    /// Return true if `id` is a leaf node.
    #[inline]
    pub fn is_leaf(&self, id: NodeId) -> bool {
        self.nodes
            .get(id.index())
            .is_some_and(|node| node.is_leaf())
    }

    /// Return true if `id` is the root node.
    #[inline]
    pub fn is_root(&self, id: NodeId) -> bool {
        self.nodes
            .get(id.index())
            .is_some_and(|node| node.is_root())
    }

    /// Count direct children of `id`.
    #[inline]
    pub fn child_count(&self, id: NodeId) -> usize {
        self.children(id).count()
    }

    /// Return the label for `id`, resolving interned taxa when present.
    #[inline]
    pub fn label(&self, id: NodeId) -> Option<&str> {
        let node = self.nodes.get(id.index())?;
        if (node.flags & FLAG_INTERNED_LABEL) != 0 {
            return self.taxa.as_ref()?.get(TaxonId(node.label_start));
        }
        if node.label_len == 0 {
            if (node.flags & FLAG_QUOTED_LABEL) != 0 {
                return Some("");
            }
            return None;
        }
        let start = node.label_start as usize;
        let end = start + node.label_len as usize;
        self.labels.get(start..end)
    }

    /// Return true if the node's label was explicitly quoted in the input.
    #[inline]
    pub fn is_quoted_label(&self, id: NodeId) -> bool {
        self.nodes
            .get(id.index())
            .is_some_and(|node| (node.flags & FLAG_QUOTED_LABEL) != 0)
    }

    /// Return the branch length leading into `id`, if present.
    #[inline]
    pub fn branch_length(&self, id: NodeId) -> Option<f64> {
        self.nodes.get(id.index())?.branch_length()
    }

    /// Return node-level annotations (`[...]` before `:`) for `id`.
    #[inline]
    pub fn node_annotations(&self, id: NodeId) -> &[Annotation] {
        let Some(node) = self.nodes.get(id.index()) else {
            return &[];
        };
        if (node.flags & FLAG_HAS_NODE_ANNOT) == 0 {
            return &[];
        }
        self.annotations
            .as_deref()
            .and_then(|s| s.node_annotations.get(&id))
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    /// Return branch-level annotations (`[...]` after `:`) for `id`.
    #[inline]
    pub fn branch_annotations(&self, id: NodeId) -> &[Annotation] {
        let Some(node) = self.nodes.get(id.index()) else {
            return &[];
        };
        if (node.flags & FLAG_HAS_BRANCH_ANNOT) == 0 {
            return &[];
        }
        self.annotations
            .as_deref()
            .and_then(|s| s.branch_annotations.get(&id))
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    /// Return Extended Newick hybrid metadata for `id`, if `id` is a hybrid node.
    #[inline]
    pub fn hybrid_info(&self, id: NodeId) -> Option<&HybridInfo> {
        let node = self.nodes.get(id.index())?;
        if (node.flags & FLAG_IS_HYBRID) == 0 {
            return None;
        }
        self.annotations.as_deref()?.hybrids.get(&id)
    }

    /// Return true if this tree contains Extended Newick reticulation nodes.
    #[inline]
    pub fn is_network(&self) -> bool {
        self.annotations
            .as_deref()
            .is_some_and(|s| !s.hybrids.is_empty())
    }

    /// Return all Extended Newick hybrid nodes in arena order.
    pub fn hybrid_nodes(&self) -> Vec<(NodeId, &HybridInfo)> {
        let Some(store) = self.annotations.as_deref() else {
            return Vec::new();
        };
        let mut items: Vec<(NodeId, &HybridInfo)> =
            store.hybrids.iter().map(|(&id, info)| (id, info)).collect();
        items.sort_unstable_by_key(|(id, _)| *id);
        items
    }

    /// Group Extended Newick hybrid nodes by their `#tag`.
    pub fn reticulation_groups(&self) -> AHashMap<&str, Vec<NodeId>> {
        let mut groups: AHashMap<&str, Vec<NodeId>> = AHashMap::new();
        for (id, info) in self.hybrid_nodes() {
            if !info.tag.is_empty() {
                groups.entry(info.tag.as_str()).or_default().push(id);
            }
        }
        groups
    }

    /// Extract a primary numeric branch support value for `id`, if present.
    ///
    /// Checks in order:
    /// 1. Extended Newick `:len:support` edge support
    /// 2. Raw numeric branch annotation (e.g. `:0.05[95]`)
    /// 3. BEAST `posterior`/`prob` or NHX `B` annotation on branch or node
    /// 4. Numeric internal node label (e.g. `)95:0.05`)
    pub fn support_value(&self, id: NodeId) -> Option<f64> {
        if let Some(sup) = self.hybrid_info(id).and_then(|h| h.edge_support) {
            return Some(sup);
        }
        for annot in self.branch_annotations(id) {
            if let Some(val) = annotation_support(annot) {
                return Some(val);
            }
        }
        for annot in self.node_annotations(id) {
            if let Some(val) = annotation_support(annot) {
                return Some(val);
            }
        }
        if !self.is_leaf(id)
            && let Some(lbl) = self.label(id)
        {
            return lbl.trim().parse::<f64>().ok().filter(|v| v.is_finite());
        }
        None
    }

    /// Parse `/`-delimited multi-support strings on internal nodes (e.g. `95.4/99`).
    pub fn support_values(&self, id: NodeId) -> Vec<f64> {
        if self.is_leaf(id) {
            return Vec::new();
        }
        let Some(lbl) = self.label(id) else {
            return Vec::new();
        };
        let parts: Vec<f64> = lbl
            .split('/')
            .filter_map(|p| p.trim().parse::<f64>().ok().filter(|v| v.is_finite()))
            .collect();
        if parts.is_empty() { Vec::new() } else { parts }
    }

    /// Stack-safe iterative preorder traversal (root before children, left-to-right).
    pub fn preorder(&self) -> PreorderIter<'_> {
        let mut stack = Vec::new();
        if !self.nodes.is_empty() {
            stack.push(self.root);
        }
        PreorderIter { tree: self, stack }
    }

    /// Stack-safe iterative postorder traversal (children left-to-right before parent).
    pub fn postorder(&self) -> PostorderIter<'_> {
        let mut stack = Vec::new();
        if !self.nodes.is_empty() {
            stack.push((self.root, false));
        }
        PostorderIter { tree: self, stack }
    }

    /// Stack-safe iterator over leaf nodes in Newick left-to-right order.
    pub fn leaves(&self) -> impl Iterator<Item = NodeId> + '_ {
        self.preorder().filter(|&id| self.is_leaf(id))
    }
}

fn annotation_support(annot: &Annotation) -> Option<f64> {
    match annot.format {
        AnnotationFormat::Raw => annot
            .raw
            .trim()
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite()),
        AnnotationFormat::Beast => {
            for (k, v) in &annot.entries {
                if k.eq_ignore_ascii_case("posterior")
                    || k.eq_ignore_ascii_case("prob")
                    || k.eq_ignore_ascii_case("confidence")
                    || k.eq_ignore_ascii_case("bootstrap")
                {
                    match v {
                        AnnotationValue::Float(x) if x.is_finite() => return Some(*x),
                        AnnotationValue::Int(x) => return Some(*x as f64),
                        AnnotationValue::Str(s) => {
                            if let Ok(x) = s.trim().parse::<f64>()
                                && x.is_finite()
                            {
                                return Some(x);
                            }
                        }
                        _ => {}
                    }
                }
            }
            None
        }
        AnnotationFormat::Nhx => {
            for (k, v) in &annot.entries {
                if k == "B" {
                    match v {
                        AnnotationValue::Float(x) if x.is_finite() => return Some(*x),
                        AnnotationValue::Int(x) => return Some(*x as f64),
                        AnnotationValue::Str(s) => {
                            if let Ok(x) = s.trim().parse::<f64>()
                                && x.is_finite()
                            {
                                return Some(x);
                            }
                        }
                        _ => {}
                    }
                }
            }
            None
        }
    }
}

/// Iterator over the direct children of a node.
pub struct ChildrenIter<'a> {
    tree: &'a Tree,
    next: Option<NodeId>,
}

impl Iterator for ChildrenIter<'_> {
    type Item = NodeId;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let current = self.next?;
        self.next = self.tree.next_sibling(current);
        Some(current)
    }
}

/// Iterative preorder node iterator.
pub struct PreorderIter<'a> {
    tree: &'a Tree,
    stack: Vec<NodeId>,
}

impl Iterator for PreorderIter<'_> {
    type Item = NodeId;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let id = self.stack.pop()?;
        let base = self.stack.len();
        for child in self.tree.children(id) {
            self.stack.push(child);
        }
        self.stack[base..].reverse();
        Some(id)
    }
}

/// Iterative postorder node iterator.
pub struct PostorderIter<'a> {
    tree: &'a Tree,
    stack: Vec<(NodeId, bool)>,
}

impl Iterator for PostorderIter<'_> {
    type Item = NodeId;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        while let Some((id, visited)) = self.stack.pop() {
            if visited || self.tree.is_leaf(id) {
                return Some(id);
            }
            self.stack.push((id, true));
            let base = self.stack.len();
            for child in self.tree.children(id) {
                self.stack.push((child, false));
            }
            self.stack[base..].reverse();
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_node_size_is_32_bytes() {
        assert_eq!(std::mem::size_of::<Node>(), 32);
    }

    #[test]
    fn test_lcrs_traversals_and_counts() {
        let mut tree = Tree::with_capacity(5, 16);
        let root = tree.add_node();
        tree.set_root(root);
        let internal = tree.add_node();
        let a = tree.add_node();
        let b = tree.add_node();
        let c = tree.add_node();

        tree.attach_child(root, internal, None);
        tree.attach_child(root, c, Some(internal));
        tree.attach_child(internal, a, None);
        tree.attach_child(internal, b, Some(a));

        tree.set_label(a, "A", false);
        tree.set_label(b, "B", false);
        tree.set_label(c, "C", false);
        tree.set_label(internal, "95/99", false);

        assert_eq!(tree.len(), 5);
        assert_eq!(tree.tip_count(), 3);
        assert_eq!(
            tree.preorder().collect::<Vec<_>>(),
            vec![root, internal, a, b, c]
        );
        assert_eq!(
            tree.postorder().collect::<Vec<_>>(),
            vec![a, b, internal, c, root]
        );
        assert_eq!(
            tree.leaves()
                .filter_map(|id| tree.label(id))
                .collect::<Vec<_>>(),
            vec!["A", "B", "C"]
        );
        assert_eq!(tree.support_values(internal), vec![95.0, 99.0]);
        assert_eq!(tree.node_id(internal.index()), Some(internal));
        assert_eq!(tree.node_id(5), None);
    }
}
