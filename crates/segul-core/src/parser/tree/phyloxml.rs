//! Streaming iterative PhyloXML parser.

use std::io::BufRead;

use anyhow::{Context, Result, bail};

use super::model::{
    Annotation, AnnotationFormat, AnnotationTarget, AnnotationValue, NodeId, Rootedness, Tree,
};
use super::newick::ParseOptions;
use super::xml_pull::{XmlEvent, XmlPullReader, attr_value};

type TaxonomyInfo = (Vec<(String, AnnotationValue)>, Option<String>);

/// Streaming PhyloXML reader yielding one [`Tree`] per `<phylogeny>` element.
pub struct PhyloXmlReader<R: BufRead> {
    xml: XmlPullReader<R>,
    opts: ParseOptions,
    eof: bool,
}

impl<R: BufRead> PhyloXmlReader<R> {
    /// Create a new streaming PhyloXML reader.
    pub fn new(reader: R, opts: ParseOptions) -> Self {
        Self {
            xml: XmlPullReader::new(reader),
            opts,
            eof: false,
        }
    }

    /// Skip up to `count` `<phylogeny>` trees without constructing topology.
    pub fn skip_trees(&mut self, count: usize) -> Result<usize> {
        let mut skipped = 0usize;
        while skipped < count {
            let Some((_, empty)) = self.advance_to_phylogeny()? else {
                break;
            };
            if !empty {
                self.xml.skip_element()?;
            }
            skipped += 1;
        }
        Ok(skipped)
    }

    fn advance_to_phylogeny(&mut self) -> Result<Option<(Option<String>, bool)>> {
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
                    if name == "phylogeny" {
                        let rooted = attr_value(&attrs, "rooted").map(str::to_string);
                        return Ok(Some((rooted, empty)));
                    }
                    if name != "phyloxml" && !empty {
                        self.xml.skip_element()?;
                    }
                }
                XmlEvent::End(_) | XmlEvent::Text(_) => {}
            }
        }
    }

    fn read_next_phylogeny(&mut self) -> Result<Option<Tree>> {
        let Some((rooted_attr, empty)) = self.advance_to_phylogeny()? else {
            return Ok(None);
        };
        if empty {
            bail!("PhyloXML <phylogeny> element has no <clade>");
        }

        let mut tree = Tree::new();
        if let Some(rooted) = rooted_attr.as_deref() {
            if rooted.eq_ignore_ascii_case("true") || rooted == "1" {
                tree.metadata.rootedness = Rootedness::Rooted;
            } else if rooted.eq_ignore_ascii_case("false") || rooted == "0" {
                tree.metadata.rootedness = Rootedness::Unrooted;
            }
        }

        struct CladeFrame {
            node: NodeId,
            last_child: Option<NodeId>,
            has_name: bool,
            fallback_label: Option<String>,
        }

        let mut stack: Vec<CladeFrame> = Vec::with_capacity(64);
        let mut has_root_clade = false;

        loop {
            match self.xml.next_event()? {
                XmlEvent::Start { name, attrs, empty } => {
                    if stack.is_empty() {
                        match name.as_str() {
                            "clade" => {
                                if has_root_clade {
                                    bail!(
                                        "PhyloXML <phylogeny> has multiple root <clade> elements"
                                    );
                                }
                                let node = tree.add_node();
                                tree.set_root(node);
                                has_root_clade = true;
                                if let Some(len_str) = attr_value(&attrs, "branch_length") {
                                    let len = parse_phylo_float(
                                        len_str,
                                        "branch_length",
                                        self.opts.allow_negative_branch_lengths,
                                    )?;
                                    tree.set_branch_length(node, len);
                                }
                                if !empty {
                                    stack.push(CladeFrame {
                                        node,
                                        last_child: None,
                                        has_name: false,
                                        fallback_label: None,
                                    });
                                }
                            }
                            "name" => {
                                if !empty {
                                    let text = self.xml.read_element_text()?;
                                    let trimmed = text.trim();
                                    if !trimmed.is_empty() {
                                        tree.metadata.name = Some(trimmed.to_string());
                                    }
                                }
                            }
                            "description" => {
                                if !empty {
                                    let text = self.xml.read_element_text()?;
                                    let trimmed = text.trim();
                                    if !trimmed.is_empty() {
                                        tree.metadata.comments.push(trimmed.to_string());
                                    }
                                }
                            }
                            _ => {
                                if !empty {
                                    self.xml.skip_element()?;
                                }
                            }
                        }
                    } else {
                        match name.as_str() {
                            "clade" => {
                                let child = tree.add_node();
                                let parent = stack.last_mut().expect("stack is non-empty");
                                tree.attach_child(parent.node, child, parent.last_child);
                                parent.last_child = Some(child);
                                if let Some(len_str) = attr_value(&attrs, "branch_length") {
                                    let len = parse_phylo_float(
                                        len_str,
                                        "branch_length",
                                        self.opts.allow_negative_branch_lengths,
                                    )?;
                                    tree.set_branch_length(child, len);
                                }
                                if !empty {
                                    stack.push(CladeFrame {
                                        node: child,
                                        last_child: None,
                                        has_name: false,
                                        fallback_label: None,
                                    });
                                }
                            }
                            "name" => {
                                if !empty {
                                    let text = self.xml.read_element_text()?;
                                    let trimmed = text.trim();
                                    if !trimmed.is_empty() {
                                        let cur = stack.last_mut().expect("stack is non-empty");
                                        tree.set_label(cur.node, trimmed, false);
                                        cur.has_name = true;
                                    }
                                }
                            }
                            "branch_length" => {
                                if !empty {
                                    let text = self.xml.read_element_text()?;
                                    let trimmed = text.trim();
                                    if !trimmed.is_empty() {
                                        let len = parse_phylo_float(
                                            trimmed,
                                            "branch_length",
                                            self.opts.allow_negative_branch_lengths,
                                        )?;
                                        let cur = stack.last().expect("stack is non-empty");
                                        tree.set_branch_length(cur.node, len);
                                    }
                                }
                            }
                            "confidence" => {
                                if !empty {
                                    let text = self.xml.read_element_text()?;
                                    let trimmed = text.trim();
                                    if self.opts.parse_annotations && !trimmed.is_empty() {
                                        let mut entries =
                                            vec![("confidence".to_string(), parse_scalar(trimmed))];
                                        if let Some(kind) = attr_value(&attrs, "type") {
                                            entries.push((
                                                "confidence_type".to_string(),
                                                AnnotationValue::Str(kind.to_string()),
                                            ));
                                        }
                                        let raw = format_entries_raw(&entries);
                                        let cur = stack.last().expect("stack is non-empty");
                                        tree.add_annotation(
                                            cur.node,
                                            AnnotationTarget::Node,
                                            Annotation {
                                                format: AnnotationFormat::Beast,
                                                entries,
                                                raw,
                                            },
                                        );
                                    }
                                }
                            }
                            "property" => {
                                if !empty {
                                    let text = self.xml.read_element_text()?;
                                    let trimmed = text.trim();
                                    if self.opts.parse_annotations && !trimmed.is_empty() {
                                        let key = attr_value(&attrs, "ref")
                                            .unwrap_or("property")
                                            .to_string();
                                        let target = if attr_value(&attrs, "applies_to")
                                            .is_some_and(|v| v == "parent_branch")
                                        {
                                            AnnotationTarget::Branch
                                        } else {
                                            AnnotationTarget::Node
                                        };
                                        let entries = vec![(key, parse_scalar(trimmed))];
                                        let raw = format_entries_raw(&entries);
                                        let cur = stack.last().expect("stack is non-empty");
                                        tree.add_annotation(
                                            cur.node,
                                            target,
                                            Annotation {
                                                format: AnnotationFormat::Beast,
                                                entries,
                                                raw,
                                            },
                                        );
                                    }
                                }
                            }
                            "taxonomy" => {
                                if !empty {
                                    let (entries, fallback) = self.read_taxonomy()?;
                                    let cur = stack.last_mut().expect("stack is non-empty");
                                    if cur.fallback_label.is_none() {
                                        cur.fallback_label = fallback;
                                    }
                                    if self.opts.parse_annotations && !entries.is_empty() {
                                        let raw = format_entries_raw(&entries);
                                        tree.add_annotation(
                                            cur.node,
                                            AnnotationTarget::Node,
                                            Annotation {
                                                format: AnnotationFormat::Beast,
                                                entries,
                                                raw,
                                            },
                                        );
                                    }
                                }
                            }
                            "events" => {
                                if !empty {
                                    let entries = self.read_events()?;
                                    if self.opts.parse_annotations && !entries.is_empty() {
                                        let raw = format_entries_raw(&entries);
                                        let cur = stack.last().expect("stack is non-empty");
                                        tree.add_annotation(
                                            cur.node,
                                            AnnotationTarget::Node,
                                            Annotation {
                                                format: AnnotationFormat::Beast,
                                                entries,
                                                raw,
                                            },
                                        );
                                    }
                                }
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
                    if name == "clade" {
                        let frame = stack.pop().context("Unmatched </clade> in PhyloXML")?;
                        if !frame.has_name
                            && tree.is_leaf(frame.node)
                            && let Some(lbl) = frame.fallback_label
                        {
                            tree.set_label(frame.node, &lbl, false);
                        }
                    } else if name == "phylogeny" {
                        if !stack.is_empty() {
                            bail!("Unclosed <clade> inside PhyloXML <phylogeny>");
                        }
                        if !has_root_clade {
                            bail!("PhyloXML <phylogeny> element has no <clade>");
                        }
                        return Ok(Some(tree));
                    }
                }
                XmlEvent::Eof => bail!("Unexpected EOF inside PhyloXML <phylogeny>"),
                XmlEvent::Text(_) => {}
            }
        }
    }

    fn read_taxonomy(&mut self) -> Result<TaxonomyInfo> {
        let mut entries = Vec::new();
        let mut sci_name: Option<String> = None;
        let mut code: Option<String> = None;

        loop {
            match self.xml.next_event()? {
                XmlEvent::Start { name, attrs, empty } => {
                    if empty {
                        continue;
                    }
                    match name.as_str() {
                        "scientific_name" => {
                            let val = self.xml.read_element_text()?;
                            let trimmed = val.trim();
                            if !trimmed.is_empty() {
                                sci_name = Some(trimmed.to_string());
                                entries.push((
                                    "taxonomy.scientific_name".to_string(),
                                    AnnotationValue::Str(trimmed.to_string()),
                                ));
                            }
                        }
                        "code" => {
                            let val = self.xml.read_element_text()?;
                            let trimmed = val.trim();
                            if !trimmed.is_empty() {
                                code = Some(trimmed.to_string());
                                entries.push((
                                    "taxonomy.code".to_string(),
                                    AnnotationValue::Str(trimmed.to_string()),
                                ));
                            }
                        }
                        "common_name" | "rank" => {
                            let val = self.xml.read_element_text()?;
                            let trimmed = val.trim();
                            if !trimmed.is_empty() {
                                entries.push((
                                    format!("taxonomy.{name}"),
                                    AnnotationValue::Str(trimmed.to_string()),
                                ));
                            }
                        }
                        "id" => {
                            let provider = attr_value(&attrs, "provider").map(str::to_string);
                            let val = self.xml.read_element_text()?;
                            let trimmed = val.trim();
                            if !trimmed.is_empty() {
                                entries.push(("taxonomy.id".to_string(), parse_scalar(trimmed)));
                            }
                            if let Some(prov) = provider {
                                entries.push((
                                    "taxonomy.id_provider".to_string(),
                                    AnnotationValue::Str(prov),
                                ));
                            }
                        }
                        _ => self.xml.skip_element()?,
                    }
                }
                XmlEvent::End(name) if name == "taxonomy" => break,
                XmlEvent::Eof => bail!("Unexpected EOF inside PhyloXML <taxonomy>"),
                _ => {}
            }
        }

        Ok((entries, sci_name.or(code)))
    }

    fn read_events(&mut self) -> Result<Vec<(String, AnnotationValue)>> {
        let mut entries = Vec::new();
        loop {
            match self.xml.next_event()? {
                XmlEvent::Start { name, empty, .. } => {
                    if empty {
                        continue;
                    }
                    match name.as_str() {
                        "speciations" | "duplications" | "losses" | "type" => {
                            let val = self.xml.read_element_text()?;
                            let trimmed = val.trim();
                            if !trimmed.is_empty() {
                                entries.push((format!("events.{name}"), parse_scalar(trimmed)));
                            }
                        }
                        _ => self.xml.skip_element()?,
                    }
                }
                XmlEvent::End(name) if name == "events" => break,
                XmlEvent::Eof => bail!("Unexpected EOF inside PhyloXML <events>"),
                _ => {}
            }
        }
        Ok(entries)
    }
}

impl<R: BufRead> Iterator for PhyloXmlReader<R> {
    type Item = Result<Tree>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.read_next_phylogeny() {
            Ok(Some(tree)) => Some(Ok(tree)),
            Ok(None) => None,
            Err(err) => {
                self.eof = true;
                Some(Err(err))
            }
        }
    }
}

fn parse_phylo_float(raw: &str, field: &str, allow_negative: bool) -> Result<f64> {
    let trimmed = raw.trim();
    let val: f64 = trimmed
        .parse()
        .with_context(|| format!("Invalid PhyloXML {} '{}'", field, trimmed))?;
    if !val.is_finite() {
        bail!("Non-finite PhyloXML {} '{}'", field, trimmed);
    }
    if !allow_negative && val < 0.0 {
        bail!("Negative PhyloXML {} '{}'", field, trimmed);
    }
    Ok(val)
}

fn parse_scalar(raw: &str) -> AnnotationValue {
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

fn format_entries_raw(entries: &[(String, AnnotationValue)]) -> String {
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_phyloxml_streaming_annotations_and_taxonomy_fallback() {
        let xml = r#"<?xml version="1.0" encoding="UTF-8"?>
<phyloxml xmlns="http://www.phyloxml.org">
  <phylogeny rooted="true">
    <name>tree_one</name>
    <description>Primate gene tree</description>
    <clade>
      <clade branch_length="0.1">
        <confidence type="bootstrap">98</confidence>
        <property ref="rate" applies_to="parent_branch">0.002</property>
        <clade branch_length="0.05">
          <name>Human &amp; Chimp</name>
        </clade>
        <clade>
          <branch_length>0.07</branch_length>
          <taxonomy>
            <scientific_name>Pan troglodytes</scientific_name>
          </taxonomy>
        </clade>
      </clade>
      <clade branch_length="0.25">
        <name>Gorilla</name>
      </clade>
    </clade>
  </phylogeny>
  <phylogeny rooted="false">
    <name>tree_two</name>
    <clade>
      <clade branch_length="0.3"><name>A</name></clade>
      <clade branch_length="0.4"><name>B</name></clade>
    </clade>
  </phylogeny>
</phyloxml>"#;
        let trees: Vec<Tree> =
            PhyloXmlReader::new(Cursor::new(xml.as_bytes()), ParseOptions::default())
                .collect::<Result<_>>()
                .unwrap();
        assert_eq!(trees.len(), 2);
        let t1 = &trees[0];
        assert_eq!(t1.metadata.name.as_deref(), Some("tree_one"));
        assert_eq!(t1.metadata.rootedness, Rootedness::Rooted);
        assert_eq!(t1.metadata.comments, vec!["Primate gene tree"]);
        assert_eq!(
            t1.leaves()
                .filter_map(|id| t1.label(id))
                .collect::<Vec<_>>(),
            vec!["Human & Chimp", "Pan troglodytes", "Gorilla"]
        );
        let internal = t1.first_child(t1.root()).unwrap();
        assert_eq!(t1.branch_length(internal), Some(0.1));
        assert_eq!(t1.support_value(internal), Some(98.0));
        assert_eq!(t1.branch_annotations(internal).len(), 1);
        assert_eq!(trees[1].metadata.rootedness, Rootedness::Unrooted);
    }

    #[test]
    fn test_phyloxml_skip_and_deep_pectinate_50k() {
        let depth = 50_000usize;
        let mut xml = String::with_capacity(depth * 40);
        xml.push_str(
            "<phyloxml><phylogeny><clade><clade><name>skip</name></clade></clade></phylogeny>",
        );
        xml.push_str("<phylogeny>");
        for _ in 0..depth {
            xml.push_str("<clade>");
        }
        xml.push_str("<clade branch_length=\"0.1\"><name>T0</name></clade>");
        for i in 1..=depth {
            xml.push_str(&format!(
                "<clade branch_length=\"0.1\"><name>T{i}</name></clade></clade>"
            ));
        }
        xml.push_str("</phylogeny></phyloxml>");

        let mut reader = PhyloXmlReader::new(Cursor::new(xml.as_bytes()), ParseOptions::default());
        assert_eq!(reader.skip_trees(1).unwrap(), 1);
        let tree = reader.next().unwrap().unwrap();
        assert_eq!(tree.tip_count(), depth + 1);
        assert!(reader.next().is_none());
    }

    #[test]
    fn test_phyloxml_malformed_errors() {
        for bad in [
            "<phyloxml><phylogeny/></phyloxml>",
            "<phyloxml><phylogeny><name>empty</name></phylogeny></phyloxml>",
            "<phyloxml><phylogeny><clade branch_length=\"bad\"/></phylogeny></phyloxml>",
        ] {
            let mut reader =
                PhyloXmlReader::new(Cursor::new(bad.as_bytes()), ParseOptions::default());
            assert!(
                reader.next().unwrap().is_err(),
                "Expected error for {:?}",
                bad
            );
        }
    }
}
