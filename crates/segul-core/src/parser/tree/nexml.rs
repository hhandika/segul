//! Streaming iterative NeXML parser with `<otus>` interning and `<characters>` skipping.

use std::io::BufRead;
use std::sync::Arc;

use ahash::AHashMap;
use anyhow::{Context, Result, bail};

use super::model::{
    Annotation, AnnotationFormat, AnnotationTarget, AnnotationValue, HybridInfo, NodeId,
    Rootedness, TaxonId, TaxonTable, Tree,
};
use super::newick::ParseOptions;
use super::xml_pull::{XmlEvent, XmlPullReader, attr_value};

#[derive(Debug, Clone)]
struct OtuMap {
    taxa: Arc<TaxonTable>,
    by_id: AHashMap<String, TaxonId>,
}

struct RawNode {
    id: String,
    label: Option<String>,
    otu: Option<String>,
    is_root: bool,
    annots: Vec<(String, AnnotationValue)>,
}

struct RawEdge {
    source: String,
    target: String,
    length: Option<f64>,
    annots: Vec<(String, AnnotationValue)>,
}

struct NextTreeHeader {
    is_network: bool,
    id: Option<String>,
    label: Option<String>,
    empty: bool,
}

/// Streaming NeXML reader yielding one [`Tree`] per `<tree>` or `<network>` element.
pub struct NeXmlReader<R: BufRead> {
    xml: XmlPullReader<R>,
    opts: ParseOptions,
    otus_maps: AHashMap<String, Arc<OtuMap>>,
    in_trees_block: bool,
    current_otus: Option<Arc<OtuMap>>,
    eof: bool,
}

impl<R: BufRead> NeXmlReader<R> {
    /// Create a new streaming NeXML reader.
    pub fn new(reader: R, opts: ParseOptions) -> Self {
        Self {
            xml: XmlPullReader::new(reader),
            opts,
            otus_maps: AHashMap::new(),
            in_trees_block: false,
            current_otus: None,
            eof: false,
        }
    }

    /// Skip up to `count` `<tree>` / `<network>` elements without constructing topology.
    pub fn skip_trees(&mut self, count: usize) -> Result<usize> {
        let mut skipped = 0usize;
        while skipped < count {
            let Some(header) = self.advance_to_next_tree()? else {
                break;
            };
            if !header.empty {
                self.xml.skip_element()?;
            }
            skipped += 1;
        }
        Ok(skipped)
    }

    /// Advance until the next `<tree>` or `<network>` inside a `<trees>` block.
    fn advance_to_next_tree(&mut self) -> Result<Option<NextTreeHeader>> {
        if self.eof {
            return Ok(None);
        }

        loop {
            match self.xml.next_event()? {
                XmlEvent::Eof => {
                    self.eof = true;
                    return Ok(None);
                }
                XmlEvent::Start { name, attrs, empty } => {
                    if !self.in_trees_block {
                        match name.as_str() {
                            "nexml" => {
                                if empty {
                                    self.eof = true;
                                    return Ok(None);
                                }
                            }
                            "otus" => {
                                if !empty {
                                    self.parse_otus_block(&attrs)?;
                                }
                            }
                            "trees" => {
                                if !empty {
                                    let otus_ref = attr_value(&attrs, "otus");
                                    self.current_otus = otus_ref
                                        .and_then(|id| self.otus_maps.get(id).cloned())
                                        .or_else(|| {
                                            if self.otus_maps.len() == 1 {
                                                self.otus_maps.values().next().cloned()
                                            } else {
                                                None
                                            }
                                        });
                                    self.in_trees_block = true;
                                }
                            }
                            _ => {
                                // Skip <characters>, <meta>, or other non-tree blocks cleanly.
                                if !empty {
                                    self.xml.skip_element()?;
                                }
                            }
                        }
                    } else {
                        match name.as_str() {
                            "tree" | "network" => {
                                return Ok(Some(NextTreeHeader {
                                    is_network: name == "network",
                                    id: attr_value(&attrs, "id").map(str::to_string),
                                    label: attr_value(&attrs, "label").map(str::to_string),
                                    empty,
                                }));
                            }
                            _ => {
                                if !empty {
                                    self.xml.skip_element()?;
                                }
                            }
                        }
                    }
                }
                XmlEvent::End(name) => {
                    if self.in_trees_block && name == "trees" {
                        self.in_trees_block = false;
                        self.current_otus = None;
                    }
                }
                XmlEvent::Text(_) => {}
            }
        }
    }

    fn parse_otus_block(&mut self, attrs: &[(String, String)]) -> Result<()> {
        let otus_id = attr_value(attrs, "id").unwrap_or("").to_string();
        let mut table = TaxonTable::new();
        let mut by_id = AHashMap::new();

        loop {
            match self.xml.next_event()? {
                XmlEvent::Start { name, attrs, empty } => {
                    if name == "otu" {
                        let id = attr_value(&attrs, "id")
                            .context("NeXML <otu> missing 'id' attribute")?
                            .to_string();
                        let label = attr_value(&attrs, "label")
                            .filter(|s| !s.is_empty())
                            .unwrap_or(&id);
                        let taxon_id = table.intern(label);
                        by_id.insert(id, taxon_id);
                        if !empty {
                            self.xml.skip_element()?;
                        }
                    } else if !empty {
                        self.xml.skip_element()?;
                    }
                }
                XmlEvent::End(name) if name == "otus" => break,
                XmlEvent::Eof => bail!("Unexpected EOF inside NeXML <otus>"),
                _ => {}
            }
        }

        let map = Arc::new(OtuMap {
            taxa: Arc::new(table),
            by_id,
        });
        self.otus_maps.insert(otus_id, map);
        Ok(())
    }

    fn read_next_tree(&mut self) -> Result<Option<Tree>> {
        let Some(header) = self.advance_to_next_tree()? else {
            return Ok(None);
        };
        if header.empty {
            bail!("Empty NeXML <tree> element has no nodes");
        }
        let is_network = header.is_network;
        let tree_id = header.id;
        let tree_label = header.label;

        let end_tag = if is_network { "network" } else { "tree" };
        let mut raw_nodes: Vec<RawNode> = Vec::new();
        let mut raw_edges: Vec<RawEdge> = Vec::new();
        let mut root_edge_len: Option<f64> = None;
        let mut root_edge_target: Option<String> = None;

        #[derive(Clone, Copy)]
        enum ActiveTarget {
            Node(usize),
            Edge(usize),
        }
        let mut active: Option<ActiveTarget> = None;

        loop {
            match self.xml.next_event()? {
                XmlEvent::Start { name, attrs, empty } => match name.as_str() {
                    "node" => {
                        let id = attr_value(&attrs, "id")
                            .context("NeXML <node> missing 'id' attribute")?
                            .to_string();
                        let label = attr_value(&attrs, "label").map(str::to_string);
                        let otu = attr_value(&attrs, "otu").map(str::to_string);
                        let is_root = attr_value(&attrs, "root")
                            .is_some_and(|v| v.eq_ignore_ascii_case("true") || v == "1");
                        let idx = raw_nodes.len();
                        raw_nodes.push(RawNode {
                            id,
                            label,
                            otu,
                            is_root,
                            annots: Vec::new(),
                        });
                        if !empty {
                            active = Some(ActiveTarget::Node(idx));
                        }
                    }
                    "edge" => {
                        let source = attr_value(&attrs, "source")
                            .context("NeXML <edge> missing 'source' attribute")?
                            .to_string();
                        let target = attr_value(&attrs, "target")
                            .context("NeXML <edge> missing 'target' attribute")?
                            .to_string();
                        let length = parse_optional_length(
                            attr_value(&attrs, "length"),
                            self.opts.allow_negative_branch_lengths,
                        )?;
                        let idx = raw_edges.len();
                        raw_edges.push(RawEdge {
                            source,
                            target,
                            length,
                            annots: Vec::new(),
                        });
                        if !empty {
                            active = Some(ActiveTarget::Edge(idx));
                        }
                    }
                    "rootedge" => {
                        root_edge_target = attr_value(&attrs, "target").map(str::to_string);
                        root_edge_len = parse_optional_length(
                            attr_value(&attrs, "length"),
                            self.opts.allow_negative_branch_lengths,
                        )?;
                        if !empty {
                            self.xml.skip_element()?;
                        }
                    }
                    "meta" => {
                        let entry = self.read_meta_entry(&attrs, empty)?;
                        if self.opts.parse_annotations
                            && let Some(kv) = entry
                            && let Some(target) = active
                        {
                            match target {
                                ActiveTarget::Node(idx) => raw_nodes[idx].annots.push(kv),
                                ActiveTarget::Edge(idx) => raw_edges[idx].annots.push(kv),
                            }
                        }
                    }
                    _ => {
                        if !empty {
                            self.xml.skip_element()?;
                        }
                    }
                },
                XmlEvent::End(name) => {
                    if name == "node" || name == "edge" {
                        active = None;
                    } else if name == end_tag {
                        break;
                    }
                }
                XmlEvent::Eof => bail!("Unexpected EOF inside NeXML <{}>", end_tag),
                XmlEvent::Text(_) => {}
            }
        }

        let mut tree = build_nexml_tree(
            &raw_nodes,
            &raw_edges,
            root_edge_target.as_deref(),
            root_edge_len,
            is_network,
            self.current_otus.as_deref(),
            &self.opts,
        )?;
        tree.metadata.name = tree_label.filter(|s| !s.is_empty()).or(tree_id);
        Ok(Some(tree))
    }

    fn read_meta_entry(
        &mut self,
        attrs: &[(String, String)],
        empty: bool,
    ) -> Result<Option<(String, AnnotationValue)>> {
        let key = attr_value(attrs, "property")
            .or_else(|| attr_value(attrs, "rel"))
            .unwrap_or("meta")
            .to_string();
        let content = attr_value(attrs, "content").map(str::to_string);
        let text = if !empty {
            let t = self.xml.read_element_text()?;
            let trimmed = t.trim();
            if trimmed.is_empty() {
                None
            } else {
                Some(trimmed.to_string())
            }
        } else {
            None
        };
        let Some(raw_val) = content.or(text) else {
            return Ok(None);
        };
        Ok(Some((key, parse_scalar_value(&raw_val))))
    }
}

impl<R: BufRead> Iterator for NeXmlReader<R> {
    type Item = Result<Tree>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.read_next_tree() {
            Ok(Some(tree)) => Some(Ok(tree)),
            Ok(None) => None,
            Err(err) => {
                self.eof = true;
                Some(Err(err))
            }
        }
    }
}

fn build_nexml_tree(
    raw_nodes: &[RawNode],
    raw_edges: &[RawEdge],
    root_edge_target: Option<&str>,
    root_edge_len: Option<f64>,
    is_network: bool,
    otus: Option<&OtuMap>,
    opts: &ParseOptions,
) -> Result<Tree> {
    if raw_nodes.is_empty() {
        bail!("NeXML tree contains no <node> elements");
    }

    let mut id_to_idx: AHashMap<&str, usize> = AHashMap::with_capacity(raw_nodes.len());
    for (idx, n) in raw_nodes.iter().enumerate() {
        if id_to_idx.insert(n.id.as_str(), idx).is_some() {
            bail!("Duplicate NeXML node id '{}'", n.id);
        }
    }

    let mut in_degree = vec![0usize; raw_nodes.len()];
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); raw_nodes.len()];
    for (e_idx, edge) in raw_edges.iter().enumerate() {
        let &src = id_to_idx
            .get(edge.source.as_str())
            .with_context(|| format!("NeXML edge references unknown source '{}'", edge.source))?;
        let &dst = id_to_idx
            .get(edge.target.as_str())
            .with_context(|| format!("NeXML edge references unknown target '{}'", edge.target))?;
        adj[src].push(e_idx);
        in_degree[dst] += 1;
    }

    let declared_root = raw_nodes.iter().position(|n| n.is_root);
    let root_idx = if let Some(idx) = declared_root {
        idx
    } else {
        let zero_in: Vec<usize> = in_degree
            .iter()
            .enumerate()
            .filter_map(|(i, &d)| if d == 0 { Some(i) } else { None })
            .collect();
        match zero_in.as_slice() {
            [only] => *only,
            [] => bail!("NeXML graph has no root node (possible cycle)"),
            _ => bail!(
                "NeXML tree has multiple root candidates ({})",
                zero_in.len()
            ),
        }
    };

    if !is_network && in_degree.iter().any(|&d| d > 1) {
        bail!("NeXML <tree> contains a node with multiple incoming edges");
    }

    let mut tree = Tree::with_capacity(raw_nodes.len(), raw_nodes.len() * 12);
    if let Some(map) = otus {
        tree.set_taxa(Arc::clone(&map.taxa));
    }
    if declared_root.is_some() {
        tree.metadata.rootedness = Rootedness::Rooted;
    }

    let mut hybrid_tags: AHashMap<usize, String> = AHashMap::new();
    let mut next_hybrid_id = 1usize;
    for (idx, &deg) in in_degree.iter().enumerate() {
        if deg > 1 {
            hybrid_tags.insert(idx, format!("#H{}", next_hybrid_id));
            next_hybrid_id += 1;
        }
    }

    let root_node = tree.add_node();
    tree.set_root(root_node);
    apply_raw_node(&mut tree, root_node, &raw_nodes[root_idx], None, otus, opts);
    if let Some(len) = root_edge_len
        && (root_edge_target.is_none() || root_edge_target == Some(raw_nodes[root_idx].id.as_str()))
    {
        tree.set_branch_length(root_node, len);
    }

    struct StackItem {
        raw_idx: usize,
        tree_node: NodeId,
        edge_cursor: usize,
        last_child: Option<NodeId>,
    }

    let mut expanded = vec![false; raw_nodes.len()];
    let mut on_stack = vec![false; raw_nodes.len()];
    expanded[root_idx] = true;
    on_stack[root_idx] = true;

    let mut stack = vec![StackItem {
        raw_idx: root_idx,
        tree_node: root_node,
        edge_cursor: 0,
        last_child: None,
    }];

    while let Some(frame) = stack.last_mut() {
        let cur_raw = frame.raw_idx;
        if frame.edge_cursor >= adj[cur_raw].len() {
            on_stack[cur_raw] = false;
            stack.pop();
            continue;
        }

        let e_idx = adj[cur_raw][frame.edge_cursor];
        frame.edge_cursor += 1;
        let edge = &raw_edges[e_idx];
        let target_raw = id_to_idx[edge.target.as_str()];

        if on_stack[target_raw] {
            bail!(
                "Cycle detected in NeXML graph at node '{}'",
                raw_nodes[target_raw].id
            );
        }

        let child_node = tree.add_node();
        tree.attach_child(frame.tree_node, child_node, frame.last_child);
        frame.last_child = Some(child_node);

        if let Some(len) = edge.length {
            tree.set_branch_length(child_node, len);
        }
        if opts.parse_annotations && !edge.annots.is_empty() {
            let raw_str = format_meta_raw(&edge.annots);
            tree.add_annotation(
                child_node,
                AnnotationTarget::Branch,
                Annotation {
                    format: AnnotationFormat::Beast,
                    entries: edge.annots.clone(),
                    raw: raw_str,
                },
            );
        }

        let hybrid_tag = hybrid_tags.get(&target_raw).map(String::as_str);
        apply_raw_node(
            &mut tree,
            child_node,
            &raw_nodes[target_raw],
            hybrid_tag,
            otus,
            opts,
        );

        if expanded[target_raw] {
            // Secondary reticulation reference: keep as a hybrid leaf node.
            continue;
        }

        expanded[target_raw] = true;
        on_stack[target_raw] = true;
        stack.push(StackItem {
            raw_idx: target_raw,
            tree_node: child_node,
            edge_cursor: 0,
            last_child: None,
        });
    }

    if !is_network && tree.len() != raw_nodes.len() {
        bail!(
            "NeXML <tree> is disconnected (visited {} of {} nodes)",
            tree.len(),
            raw_nodes.len()
        );
    }

    Ok(tree)
}

fn apply_raw_node(
    tree: &mut Tree,
    node: NodeId,
    raw: &RawNode,
    hybrid_tag: Option<&str>,
    otus: Option<&OtuMap>,
    opts: &ParseOptions,
) {
    if let Some(tag) = hybrid_tag {
        let base_label = raw
            .otu
            .as_deref()
            .and_then(|id| {
                otus.and_then(|m| m.by_id.get(id))
                    .and_then(|&tid| otus?.taxa.get(tid))
            })
            .or(raw.label.as_deref());
        let combined = match base_label {
            Some(lbl) if !lbl.is_empty() && !lbl.contains('#') => format!("{lbl}{tag}"),
            Some(lbl) if !lbl.is_empty() => lbl.to_string(),
            _ => tag.to_string(),
        };
        tree.set_label(node, &combined, false);
        if opts.parse_extended_newick {
            tree.set_hybrid_info(
                node,
                HybridInfo {
                    tag: tag.to_string(),
                    edge_support: None,
                    inheritance_prob: None,
                },
            );
        }
    } else if let Some(otu_id) = raw.otu.as_deref()
        && let Some(map) = otus
        && let Some(&taxon_id) = map.by_id.get(otu_id)
    {
        tree.set_interned_label(node, taxon_id, false);
    } else if let Some(lbl) = raw.label.as_deref()
        && !lbl.is_empty()
    {
        tree.set_label(node, lbl, false);
    }

    if opts.parse_annotations && !raw.annots.is_empty() {
        let raw_str = format_meta_raw(&raw.annots);
        tree.add_annotation(
            node,
            AnnotationTarget::Node,
            Annotation {
                format: AnnotationFormat::Beast,
                entries: raw.annots.clone(),
                raw: raw_str,
            },
        );
    }
}

fn format_meta_raw(entries: &[(String, AnnotationValue)]) -> String {
    let mut raw = String::from("&");
    for (i, (k, v)) in entries.iter().enumerate() {
        if i > 0 {
            raw.push(',');
        }
        raw.push_str(k);
        raw.push('=');
        raw.push_str(&v.to_string());
    }
    raw
}

fn parse_optional_length(raw: Option<&str>, allow_negative: bool) -> Result<Option<f64>> {
    let Some(s) = raw else {
        return Ok(None);
    };
    let trimmed = s.trim();
    if trimmed.is_empty() {
        return Ok(None);
    }
    let val: f64 = trimmed
        .parse()
        .with_context(|| format!("Invalid NeXML branch length '{}'", trimmed))?;
    if !val.is_finite() {
        bail!("Non-finite NeXML branch length '{}'", trimmed);
    }
    if !allow_negative && val < 0.0 {
        bail!("Negative NeXML branch length '{}'", trimmed);
    }
    Ok(Some(val))
}

fn parse_scalar_value(raw: &str) -> AnnotationValue {
    if raw.eq_ignore_ascii_case("true") {
        return AnnotationValue::Bool(true);
    }
    if raw.eq_ignore_ascii_case("false") {
        return AnnotationValue::Bool(false);
    }
    if let Ok(i) = raw.parse::<i64>() {
        return AnnotationValue::Int(i);
    }
    if let Ok(f) = raw.parse::<f64>()
        && f.is_finite()
    {
        return AnnotationValue::Float(f);
    }
    AnnotationValue::Str(raw.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_nexml_otus_skip_characters_and_trees() {
        let nexml = r#"<?xml version="1.0" encoding="UTF-8"?>
<nex:nexml xmlns:nex="http://www.nexml.org/2009" version="0.9">
  <otus id="otus1">
    <otu id="t1" label="Homo sapiens"/>
    <otu id="t2" label="Pan troglodytes"/>
    <otu id="t3" label="Gorilla gorilla"/>
  </otus>
  <characters id="chars1" otus="otus1">
    <matrix><row id="r1" otu="t1"><seq>ACGTACGT</seq></row></matrix>
  </characters>
  <trees id="trees1" otus="otus1">
    <tree id="tree1" label="species_tree">
      <rootedge target="n0" length="0.05"/>
      <node id="n0" root="true"/>
      <node id="n1"/>
      <node id="n2" otu="t1"/>
      <node id="n3" otu="t2"/>
      <node id="n4" otu="t3"/>
      <edge id="e1" source="n0" target="n1" length="0.2">
        <meta property="posterior" content="0.98"/>
      </edge>
      <edge id="e2" source="n1" target="n2" length="0.1"/>
      <edge id="e3" source="n1" target="n3" length="0.1"/>
      <edge id="e4" source="n0" target="n4" length="0.3"/>
    </tree>
  </trees>
</nex:nexml>"#;
        let trees: Vec<Tree> =
            NeXmlReader::new(Cursor::new(nexml.as_bytes()), ParseOptions::default())
                .collect::<Result<_>>()
                .unwrap();
        assert_eq!(trees.len(), 1);
        let t = &trees[0];
        assert_eq!(t.metadata.name.as_deref(), Some("species_tree"));
        assert_eq!(t.metadata.rootedness, Rootedness::Rooted);
        assert_eq!(t.branch_length(t.root()), Some(0.05));
        assert_eq!(
            t.leaves().filter_map(|id| t.label(id)).collect::<Vec<_>>(),
            vec!["Homo sapiens", "Pan troglodytes", "Gorilla gorilla"]
        );
        let n1 = t.first_child(t.root()).unwrap();
        assert_eq!(t.branch_length(n1), Some(0.2));
        assert_eq!(t.branch_annotations(n1).len(), 1);
    }

    #[test]
    fn test_nexml_network_and_skip_trees() {
        let nexml = r#"<nexml version="0.9">
  <otus id="o1"><otu id="a" label="A"/><otu id="b" label="B"/></otus>
  <trees id="ts" otus="o1">
    <tree id="t1"><node id="r" root="true"/><node id="x" otu="a"/><edge source="r" target="x" length="0.1"/></tree>
    <network id="net1" label="hyb">
      <node id="r" root="true"/>
      <node id="n1"/>
      <node id="h" label="Hyb"/>
      <node id="leaf_a" otu="a"/>
      <node id="leaf_b" otu="b"/>
      <edge source="r" target="n1" length="0.1"/>
      <edge source="n1" target="leaf_a" length="0.2"/>
      <edge source="n1" target="h" length="0.05"/>
      <edge source="r" target="h" length="0.15"/>
      <edge source="h" target="leaf_b" length="0.3"/>
    </network>
  </trees>
</nexml>"#;
        let mut reader = NeXmlReader::new(Cursor::new(nexml.as_bytes()), ParseOptions::default());
        assert_eq!(reader.skip_trees(1).unwrap(), 1);
        let net = reader.next().unwrap().unwrap();
        assert!(net.is_network());
        assert_eq!(net.reticulation_groups()["#H1"].len(), 2);
    }
}
