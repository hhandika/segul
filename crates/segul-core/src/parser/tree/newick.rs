//! Iterative byte-level Newick, NHX, and Extended Newick (eNewick) parser.

use std::borrow::Cow;
use std::io::BufRead;

use anyhow::{Context, Result, bail};
use rayon::prelude::*;

use super::comments::{
    parse_comment_annotation, parse_tree_comment_directives, scan_bracket_comment,
};
use super::model::{AnnotationTarget, HybridInfo, NodeId, Tree};
use super::nexus::TranslateTable;

/// Configuration options for parsing Newick, NHX, and Extended Newick trees.
#[derive(Debug, Clone)]
pub struct ParseOptions {
    /// Convert unquoted underscores (`_`) to spaces. Default: `false`.
    pub underscores_to_spaces: bool,
    /// Parse `[...]` comments into structured annotations. Default: `true`.
    pub parse_annotations: bool,
    /// Allow negative branch lengths (e.g. from distance-based methods). Default: `true`.
    pub allow_negative_branch_lengths: bool,
    /// Allow double-quoted labels in addition to standard single-quoted labels. Default: `true`.
    pub allow_double_quotes: bool,
    /// Parse Extended Newick `#tag` hybrid labels and `:len:support:gamma` fields. Default: `true`.
    pub parse_extended_newick: bool,
}

impl Default for ParseOptions {
    fn default() -> Self {
        Self {
            underscores_to_spaces: false,
            parse_annotations: true,
            allow_negative_branch_lengths: true,
            allow_double_quotes: true,
            parse_extended_newick: true,
        }
    }
}

/// Iterative single-pass byte parser for Newick, NHX, and Extended Newick trees.
pub struct NewickParser;

impl NewickParser {
    /// Parse a single Newick tree from a UTF-8 string.
    pub fn parse_str(input: &str, opts: &ParseOptions) -> Result<Tree> {
        Self::parse_bytes(input.as_bytes(), opts)
    }

    /// Parse a single Newick tree from raw bytes.
    pub fn parse_bytes(input: &[u8], opts: &ParseOptions) -> Result<Tree> {
        Self::parse_bytes_with_translate(input, None, opts)
    }

    /// Parse a single Newick tree, optionally resolving leaf labels via a NEXUS `TRANSLATE` table.
    pub(crate) fn parse_bytes_with_translate(
        input: &[u8],
        translate: Option<&TranslateTable>,
        opts: &ParseOptions,
    ) -> Result<Tree> {
        let mut pos = skip_utf8_bom(input);
        let est_nodes = (input.len() / 6).clamp(8, 1_048_576);
        let est_labels = (input.len() / 2).clamp(32, 4 * 1024 * 1024);
        let mut tree = Tree::with_capacity(est_nodes, est_labels);
        if let Some(table) = translate {
            tree.set_taxa(table.taxa());
        }

        // Consume leading whitespace and tree-level comments such as `[&R]` or `[&W 0.5]`.
        loop {
            pos = skip_ws(input, pos);
            if pos < input.len() && input[pos] == b'[' {
                let (inner, next) = scan_bracket_comment(input, pos)?;
                if !parse_tree_comment_directives(inner, &mut tree.metadata)
                    && opts.parse_annotations
                {
                    let trimmed = inner.trim_ascii();
                    if !trimmed.is_empty() {
                        tree.metadata
                            .comments
                            .push(String::from_utf8_lossy(trimmed).into_owned());
                    }
                }
                pos = next;
            } else {
                break;
            }
        }

        if pos >= input.len() {
            bail!("Empty Newick input");
        }

        #[derive(Clone, Copy)]
        struct Frame {
            node: NodeId,
            last_child: Option<NodeId>,
        }

        let mut stack: Vec<Frame> = Vec::with_capacity(64);
        let mut has_root = false;

        'outer: loop {
            pos = skip_ws(input, pos);
            if pos >= input.len() {
                bail!("Unexpected end of Newick input (missing ';')");
            }

            let node = tree.add_node();
            if let Some(parent) = stack.last_mut() {
                tree.attach_child(parent.node, node, parent.last_child);
                parent.last_child = Some(node);
            } else if !has_root {
                tree.set_root(node);
                has_root = true;
            } else {
                bail!("Unexpected token at byte offset {}", pos);
            }

            // Comments immediately preceding an internal child's `(` belong to that node.
            while pos < input.len() && input[pos] == b'[' {
                let (inner, next) = scan_bracket_comment(input, pos)?;
                if opts.parse_annotations
                    && let Some(annot) = parse_comment_annotation(inner)
                {
                    tree.add_annotation(node, AnnotationTarget::Node, annot);
                }
                pos = skip_ws(input, next);
            }

            if pos < input.len() && input[pos] == b'(' {
                stack.push(Frame {
                    node,
                    last_child: None,
                });
                pos += 1;
                continue 'outer;
            }

            parse_node_metadata(&mut tree, node, input, &mut pos, translate, opts, true)?;

            // Close completed subtrees until we reach `,` (next sibling) or `;` (end of tree).
            loop {
                pos = skip_ws(input, pos);
                if pos >= input.len() {
                    bail!("Unexpected end of Newick input (missing ';')");
                }
                match input[pos] {
                    b',' => {
                        if stack.is_empty() {
                            bail!("Unexpected ',' outside parentheses at byte offset {}", pos);
                        }
                        pos += 1;
                        continue 'outer;
                    }
                    b')' => {
                        let frame = stack.pop().ok_or_else(|| {
                            anyhow::anyhow!("Unmatched ')' at byte offset {}", pos)
                        })?;
                        pos += 1;
                        parse_node_metadata(
                            &mut tree, frame.node, input, &mut pos, translate, opts, false,
                        )?;
                    }
                    b';' => {
                        if !stack.is_empty() {
                            bail!(
                                "Unexpected ';' at byte offset {}: {} unclosed '('",
                                pos,
                                stack.len()
                            );
                        }
                        pos += 1;
                        pos = skip_ws_and_comments(input, pos)?;
                        if pos < input.len() {
                            bail!(
                                "Unexpected trailing token '{}' after ';' at byte offset {}",
                                input[pos] as char,
                                pos
                            );
                        }
                        return Ok(tree);
                    }
                    other => {
                        bail!(
                            "Unexpected character '{}' at byte offset {}",
                            other as char,
                            pos
                        );
                    }
                }
            }
        }
    }
}

/// Streaming reader over a Newick file containing one or more `;`-terminated trees.
pub struct NewickReader<R: BufRead> {
    reader: R,
    opts: ParseOptions,
    buf: Vec<u8>,
    eof: bool,
}

impl<R: BufRead> NewickReader<R> {
    /// Create a streaming reader with the given [`ParseOptions`].
    pub fn new(reader: R, opts: ParseOptions) -> Self {
        Self {
            reader,
            opts,
            buf: Vec::with_capacity(4096),
            eof: false,
        }
    }

    /// Skip up to `count` trees without constructing topology, returning the count skipped.
    pub fn skip_trees(&mut self, count: usize) -> Result<usize> {
        let mut skipped = 0usize;
        while skipped < count {
            if !self.read_next_raw_tree()? {
                break;
            }
            skipped += 1;
        }
        Ok(skipped)
    }

    /// Read up to `chunk_size` raw `;`-terminated tree records without parsing topologies.
    pub(crate) fn next_raw_chunk(&mut self, chunk_size: usize) -> Result<Vec<Vec<u8>>> {
        let mut records = Vec::with_capacity(chunk_size.min(64));
        while records.len() < chunk_size && self.read_next_raw_tree()? {
            records.push(self.buf.clone());
        }
        Ok(records)
    }

    /// Scan the next `;`-terminated Newick record into `self.buf`.
    /// Returns `Ok(false)` at clean EOF.
    fn read_next_raw_tree(&mut self) -> Result<bool> {
        if self.eof {
            return Ok(false);
        }
        self.buf.clear();
        let mut in_single = false;
        let mut in_double = false;
        let mut comment_depth = 0usize;
        let mut has_content = false;

        loop {
            let available = self.reader.fill_buf()?;
            if available.is_empty() {
                self.eof = true;
                if has_content {
                    bail!("Unexpected end of Newick stream before terminating ';'");
                }
                return Ok(false);
            }

            let mut consumed = 0usize;
            let mut found_semi = false;
            for (i, &b) in available.iter().enumerate() {
                if comment_depth > 0 {
                    if b == b'\'' && !in_double {
                        in_single = !in_single;
                    } else if b == b'"' && !in_single {
                        in_double = !in_double;
                    } else if !in_single && !in_double {
                        if b == b'[' {
                            comment_depth += 1;
                        } else if b == b']' {
                            comment_depth -= 1;
                        }
                    }
                    consumed = i + 1;
                    continue;
                }

                if in_single {
                    if b == b'\'' {
                        in_single = false;
                    }
                    has_content = true;
                    consumed = i + 1;
                    continue;
                }

                if in_double {
                    if b == b'"' {
                        in_double = false;
                    }
                    has_content = true;
                    consumed = i + 1;
                    continue;
                }

                match b {
                    b'\'' => {
                        in_single = true;
                        has_content = true;
                    }
                    b'"' if self.opts.allow_double_quotes => {
                        in_double = true;
                        has_content = true;
                    }
                    b'[' => {
                        comment_depth = 1;
                    }
                    b';' => {
                        consumed = i + 1;
                        found_semi = true;
                        break;
                    }
                    b if !b.is_ascii_whitespace() && !is_utf8_bom_byte(b, self.buf.len() + i) => {
                        has_content = true;
                    }
                    _ => {}
                }
                consumed = i + 1;
            }

            self.buf.extend_from_slice(&available[..consumed]);
            self.reader.consume(consumed);

            if found_semi {
                return Ok(true);
            }
        }
    }
}

impl<R: BufRead> Iterator for NewickReader<R> {
    type Item = Result<Tree>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.read_next_raw_tree() {
            Ok(true) => Some(NewickParser::parse_bytes(&self.buf, &self.opts)),
            Ok(false) => None,
            Err(err) => {
                self.eof = true;
                Some(Err(err))
            }
        }
    }
}

/// Parse all Newick trees in an in-memory byte buffer in parallel using `rayon`.
pub fn parse_newick_bytes_par(input: &[u8], opts: &ParseOptions) -> Result<Vec<Tree>> {
    let slices = split_newick_records(input, opts.allow_double_quotes)?;
    slices
        .into_par_iter()
        .map(|slice| NewickParser::parse_bytes(slice, opts))
        .collect()
}

pub(crate) fn split_newick_records(input: &[u8], allow_double_quotes: bool) -> Result<Vec<&[u8]>> {
    let mut records = Vec::new();
    let mut start = skip_utf8_bom(input);
    let mut in_single = false;
    let mut in_double = false;
    let mut comment_depth = 0usize;
    let mut has_content = false;

    for (i, &b) in input.iter().enumerate().skip(start) {
        if comment_depth > 0 {
            if b == b'\'' && !in_double {
                in_single = !in_single;
            } else if b == b'"' && !in_single {
                in_double = !in_double;
            } else if !in_single && !in_double {
                if b == b'[' {
                    comment_depth += 1;
                } else if b == b']' {
                    comment_depth -= 1;
                }
            }
            continue;
        }
        if in_single {
            if b == b'\'' {
                in_single = false;
            }
            has_content = true;
            continue;
        }
        if in_double {
            if b == b'"' {
                in_double = false;
            }
            has_content = true;
            continue;
        }
        match b {
            b'\'' => {
                in_single = true;
                has_content = true;
            }
            b'"' if allow_double_quotes => {
                in_double = true;
                has_content = true;
            }
            b'[' => {
                comment_depth = 1;
            }
            b';' => {
                records.push(&input[start..=i]);
                start = i + 1;
                has_content = false;
            }
            b if !b.is_ascii_whitespace() => {
                has_content = true;
            }
            _ => {}
        }
    }

    if comment_depth > 0 {
        bail!("Unclosed '[' comment in Newick input");
    }
    if in_single || in_double {
        bail!("Unclosed quoted label in Newick input");
    }
    if has_content {
        bail!("Unexpected trailing tokens without terminating ';'");
    }
    Ok(records)
}

fn parse_node_metadata(
    tree: &mut Tree,
    node: NodeId,
    bytes: &[u8],
    pos: &mut usize,
    translate: Option<&TranslateTable>,
    opts: &ParseOptions,
    is_leaf: bool,
) -> Result<()> {
    *pos = skip_ws(bytes, *pos);

    // Pre-label comments attach to the node.
    while *pos < bytes.len() && bytes[*pos] == b'[' {
        let (inner, next) = scan_bracket_comment(bytes, *pos)?;
        if opts.parse_annotations
            && let Some(annot) = parse_comment_annotation(inner)
        {
            tree.add_annotation(node, AnnotationTarget::Node, annot);
        }
        *pos = skip_ws(bytes, next);
    }

    // Optional node label.
    if *pos < bytes.len() && !is_newick_delimiter(bytes[*pos]) {
        let (label, quoted) = parse_label(bytes, pos, opts)?;
        if opts.parse_extended_newick
            && let Some(tag) = extract_hybrid_tag(&label)
        {
            tree.set_hybrid_info(
                node,
                HybridInfo {
                    tag: tag.to_string(),
                    edge_support: None,
                    inheritance_prob: None,
                },
            );
        }
        if is_leaf
            && let Some(table) = translate
            && let Some(taxon_id) = table.lookup(label.as_ref())
        {
            tree.set_interned_label(node, taxon_id, quoted);
        } else {
            tree.set_label(node, label.as_ref(), quoted);
        }
    }

    // Post-label, pre-colon comments also attach to the node.
    loop {
        *pos = skip_ws(bytes, *pos);
        if *pos < bytes.len() && bytes[*pos] == b'[' {
            let (inner, next) = scan_bracket_comment(bytes, *pos)?;
            if opts.parse_annotations
                && let Some(annot) = parse_comment_annotation(inner)
            {
                tree.add_annotation(node, AnnotationTarget::Node, annot);
            }
            *pos = next;
        } else {
            break;
        }
    }

    // Optional branch length (and Extended Newick :len:support:gamma fields).
    if *pos < bytes.len() && bytes[*pos] == b':' {
        *pos += 1;
        let branch_len = parse_optional_branch_float(
            tree,
            node,
            bytes,
            pos,
            opts,
            "branch length",
            opts.allow_negative_branch_lengths,
        )?;
        if let Some(len) = branch_len {
            tree.set_branch_length(node, len);
        }

        if opts.parse_extended_newick && *pos < bytes.len() && bytes[*pos] == b':' {
            *pos += 1;
            let edge_support =
                parse_optional_branch_float(tree, node, bytes, pos, opts, "edge support", true)?;
            let mut inheritance_prob = None;
            if *pos < bytes.len() && bytes[*pos] == b':' {
                *pos += 1;
                inheritance_prob = parse_optional_branch_float(
                    tree,
                    node,
                    bytes,
                    pos,
                    opts,
                    "inheritance probability",
                    false,
                )?;
            }
            if edge_support.is_some() || inheritance_prob.is_some() {
                if let Some(existing) = tree.hybrid_info_mut(node) {
                    if edge_support.is_some() {
                        existing.edge_support = edge_support;
                    }
                    if inheritance_prob.is_some() {
                        existing.inheritance_prob = inheritance_prob;
                    }
                } else {
                    tree.set_hybrid_info(
                        node,
                        HybridInfo {
                            tag: String::new(),
                            edge_support,
                            inheritance_prob,
                        },
                    );
                }
            }
        }
    }

    Ok(())
}

fn parse_optional_branch_float(
    tree: &mut Tree,
    node: NodeId,
    bytes: &[u8],
    pos: &mut usize,
    opts: &ParseOptions,
    field_name: &str,
    allow_negative: bool,
) -> Result<Option<f64>> {
    loop {
        *pos = skip_ws(bytes, *pos);
        if *pos < bytes.len() && bytes[*pos] == b'[' {
            let (inner, next) = scan_bracket_comment(bytes, *pos)?;
            if opts.parse_annotations
                && let Some(annot) = parse_comment_annotation(inner)
            {
                tree.add_annotation(node, AnnotationTarget::Branch, annot);
            }
            *pos = next;
        } else {
            break;
        }
    }

    let start = *pos;
    while *pos < bytes.len()
        && !is_newick_delimiter(bytes[*pos])
        && !bytes[*pos].is_ascii_whitespace()
    {
        *pos += 1;
    }

    let parsed = if *pos > start {
        let raw = std::str::from_utf8(&bytes[start..*pos])
            .with_context(|| format!("Invalid UTF-8 in {} at byte offset {}", field_name, start))?;
        let val: f64 = raw.parse().with_context(|| {
            format!("Invalid {} '{}' at byte offset {}", field_name, raw, start)
        })?;
        if !val.is_finite() {
            bail!(
                "Non-finite {} '{}' at byte offset {}",
                field_name,
                raw,
                start
            );
        }
        if !allow_negative && val < 0.0 {
            bail!("Negative {} '{}' at byte offset {}", field_name, raw, start);
        }
        Some(val)
    } else {
        None
    };

    loop {
        *pos = skip_ws(bytes, *pos);
        if *pos < bytes.len() && bytes[*pos] == b'[' {
            let (inner, next) = scan_bracket_comment(bytes, *pos)?;
            if opts.parse_annotations
                && let Some(annot) = parse_comment_annotation(inner)
            {
                tree.add_annotation(node, AnnotationTarget::Branch, annot);
            }
            *pos = next;
        } else {
            break;
        }
    }

    Ok(parsed)
}

fn parse_label<'a>(
    bytes: &'a [u8],
    pos: &mut usize,
    opts: &ParseOptions,
) -> Result<(Cow<'a, str>, bool)> {
    let start = *pos;
    let quote = bytes[start];
    if quote == b'\'' || (opts.allow_double_quotes && quote == b'"') {
        *pos += 1;
        let content_start = *pos;
        let mut unescaped: Option<String> = None;
        while *pos < bytes.len() {
            let b = bytes[*pos];
            if b == quote {
                if *pos + 1 < bytes.len() && bytes[*pos + 1] == quote {
                    let s = unescaped.get_or_insert_with(|| {
                        String::from_utf8_lossy(&bytes[content_start..*pos]).into_owned()
                    });
                    s.push(quote as char);
                    *pos += 2;
                    continue;
                }
                let text = if let Some(owned) = unescaped {
                    Cow::Owned(owned)
                } else {
                    let raw =
                        std::str::from_utf8(&bytes[content_start..*pos]).with_context(|| {
                            format!("Invalid UTF-8 in quoted label at byte offset {}", start)
                        })?;
                    Cow::Borrowed(raw)
                };
                *pos += 1;
                return Ok((text, true));
            }
            if let Some(s) = &mut unescaped {
                s.push(b as char);
            }
            *pos += 1;
        }
        bail!("Unclosed quoted label starting at byte offset {}", start);
    }

    while *pos < bytes.len()
        && !is_newick_delimiter(bytes[*pos])
        && !bytes[*pos].is_ascii_whitespace()
    {
        *pos += 1;
    }
    let slice = &bytes[start..*pos];
    let text = std::str::from_utf8(slice)
        .with_context(|| format!("Invalid UTF-8 in label at byte offset {}", start))?;
    if opts.underscores_to_spaces && slice.contains(&b'_') {
        Ok((Cow::Owned(text.replace('_', " ")), false))
    } else {
        Ok((Cow::Borrowed(text), false))
    }
}

fn extract_hybrid_tag(label: &str) -> Option<&str> {
    let idx = label.find('#')?;
    let tag = &label[idx..];
    if tag.len() > 1 { Some(tag) } else { None }
}

#[inline]
fn is_newick_delimiter(b: u8) -> bool {
    matches!(b, b'(' | b')' | b',' | b':' | b';' | b'[' | b']')
}

#[inline]
fn skip_ws(bytes: &[u8], mut pos: usize) -> usize {
    while pos < bytes.len() && bytes[pos].is_ascii_whitespace() {
        pos += 1;
    }
    pos
}

fn skip_ws_and_comments(bytes: &[u8], mut pos: usize) -> Result<usize> {
    loop {
        pos = skip_ws(bytes, pos);
        if pos < bytes.len() && bytes[pos] == b'[' {
            let (_, next) = scan_bracket_comment(bytes, pos)?;
            pos = next;
        } else {
            return Ok(pos);
        }
    }
}

#[inline]
fn skip_utf8_bom(bytes: &[u8]) -> usize {
    if bytes.starts_with(&[0xEF, 0xBB, 0xBF]) {
        3
    } else {
        0
    }
}

#[inline]
fn is_utf8_bom_byte(b: u8, total_offset: usize) -> bool {
    matches!((total_offset, b), (0, 0xEF) | (1, 0xBB) | (2, 0xBF))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::tree::model::AnnotationFormat;
    use std::io::Cursor;

    #[test]
    fn test_parse_basic_and_polytomy() {
        let tree = NewickParser::parse_str(
            "(A:0.1,B:0.2,(C:0.3,D:0.4,E:0.5)95:0.6);",
            &ParseOptions::default(),
        )
        .unwrap();
        assert_eq!(tree.tip_count(), 5);
        assert_eq!(tree.len(), 7);
        let leaves: Vec<&str> = tree.leaves().filter_map(|id| tree.label(id)).collect();
        assert_eq!(leaves, vec!["A", "B", "C", "D", "E"]);
        let cl_nodes: Vec<_> = tree.children(tree.root()).collect();
        assert_eq!(cl_nodes.len(), 3);
        assert_eq!(tree.branch_length(cl_nodes[0]), Some(0.1));
        assert_eq!(tree.label(cl_nodes[2]), Some("95"));
        assert_eq!(tree.support_value(cl_nodes[2]), Some(95.0));
        assert_eq!(tree.branch_length(cl_nodes[2]), Some(0.6));
    }

    #[test]
    fn test_parse_quotes_escaping_and_underscores() {
        let input = "('Homo sapiens':0.1,'Canis ''lupus'' (dog):1':1.5e-4,Pan_troglodytes:-0.02);";
        let tree = NewickParser::parse_str(input, &ParseOptions::default()).unwrap();
        let leaves: Vec<(&str, Option<f64>)> = tree
            .leaves()
            .map(|id| (tree.label(id).unwrap(), tree.branch_length(id)))
            .collect();
        assert_eq!(
            leaves,
            vec![
                ("Homo sapiens", Some(0.1)),
                ("Canis 'lupus' (dog):1", Some(1.5e-4)),
                ("Pan_troglodytes", Some(-0.02)),
            ]
        );

        let opts = ParseOptions {
            underscores_to_spaces: true,
            ..Default::default()
        };
        let tree2 = NewickParser::parse_str(input, &opts).unwrap();
        let labels2: Vec<&str> = tree2.leaves().filter_map(|id| tree2.label(id)).collect();
        assert_eq!(
            labels2,
            vec!["Homo sapiens", "Canis 'lupus' (dog):1", "Pan troglodytes"]
        );
    }

    #[test]
    fn test_node_vs_branch_annotations_and_nhx() {
        let input = "((A[&&NHX:S=human:D=N]:0.1,B:0.2)NodeAB[&height=1.2]:0.05[98],C:0.3);";
        let tree = NewickParser::parse_str(input, &ParseOptions::default()).unwrap();
        let root_children: Vec<_> = tree.children(tree.root()).collect();
        let node_ab = root_children[0];
        assert_eq!(tree.label(node_ab), Some("NodeAB"));
        assert_eq!(tree.node_annotations(node_ab).len(), 1);
        assert_eq!(
            tree.node_annotations(node_ab)[0].format,
            AnnotationFormat::Beast
        );
        assert_eq!(tree.branch_annotations(node_ab).len(), 1);
        assert_eq!(tree.branch_annotations(node_ab)[0].raw, "98");
        assert_eq!(tree.support_value(node_ab), Some(98.0));

        let a = tree.first_child(node_ab).unwrap();
        assert_eq!(tree.node_annotations(a)[0].format, AnnotationFormat::Nhx);
    }

    #[test]
    fn test_extended_newick_network() {
        let input = "((A:0.1,(B:0.2,#H1:0.05::0.3):0.1):0.2,(#H1:0.15:95:0.7,C:0.3):0.1);";
        let tree = NewickParser::parse_str(input, &ParseOptions::default()).unwrap();
        assert!(tree.is_network());
        let groups = tree.reticulation_groups();
        assert_eq!(groups.get("#H1").map(Vec::len), Some(2));
        let h_nodes = groups["#H1"].clone();
        assert_eq!(tree.branch_length(h_nodes[0]), Some(0.05));
        assert_eq!(tree.hybrid_info(h_nodes[0]).unwrap().edge_support, None);
        assert_eq!(
            tree.hybrid_info(h_nodes[0]).unwrap().inheritance_prob,
            Some(0.3)
        );
        assert_eq!(tree.branch_length(h_nodes[1]), Some(0.15));
        assert_eq!(
            tree.hybrid_info(h_nodes[1]).unwrap().edge_support,
            Some(95.0)
        );
        assert_eq!(
            tree.hybrid_info(h_nodes[1]).unwrap().inheritance_prob,
            Some(0.7)
        );
    }

    #[test]
    fn test_single_taxon_unary_and_root_branch() {
        let t1 = NewickParser::parse_str("A:0.5;", &ParseOptions::default()).unwrap();
        assert_eq!(t1.tip_count(), 1);
        assert_eq!(t1.label(t1.root()), Some("A"));
        assert_eq!(t1.branch_length(t1.root()), Some(0.5));

        let t2 = NewickParser::parse_str("((A:0.1):0.2):0.3;", &ParseOptions::default()).unwrap();
        assert_eq!(t2.tip_count(), 1);
        assert_eq!(t2.len(), 3);
        assert_eq!(t2.branch_length(t2.root()), Some(0.3));
    }

    #[test]
    fn test_deep_pectinate_100k_no_stack_overflow() {
        let depth = 100_000usize;
        let mut s = String::with_capacity(depth * 6);
        for _ in 0..depth {
            s.push('(');
        }
        s.push_str("T0:0.1");
        for i in 1..=depth {
            s.push_str(&format!(",T{i}:0.1)"));
        }
        s.push(';');

        let tree = NewickParser::parse_str(&s, &ParseOptions::default()).unwrap();
        assert_eq!(tree.tip_count(), depth + 1);
        assert_eq!(tree.preorder().count(), 2 * depth + 1);
        assert_eq!(tree.postorder().count(), 2 * depth + 1);
    }

    #[test]
    fn test_streaming_reader_skip_and_par() {
        let data = "[&R](A:0.1,B:0.2);\r\n(C:0.3,D:0.4);\n(E:0.5,F:0.6);\n";
        let mut reader = NewickReader::new(Cursor::new(data.as_bytes()), ParseOptions::default());
        assert_eq!(reader.skip_trees(1).unwrap(), 1);
        let remaining: Vec<Tree> = reader.collect::<Result<_>>().unwrap();
        assert_eq!(remaining.len(), 2);
        assert_eq!(
            remaining[0]
                .leaves()
                .filter_map(|id| remaining[0].label(id))
                .collect::<Vec<_>>(),
            vec!["C", "D"]
        );

        let par_trees = parse_newick_bytes_par(data.as_bytes(), &ParseOptions::default()).unwrap();
        assert_eq!(par_trees.len(), 3);
        assert_eq!(
            par_trees[0].metadata.rootedness,
            super::super::model::Rootedness::Rooted
        );
    }

    #[test]
    fn test_malformed_newick_errors() {
        for bad in [
            "",
            "   ",
            "(A,B)",
            "(A,B));",
            "((A,B);",
            "(A:abc,B);",
            "('unclosed,B);",
            "(A,B); extra",
        ] {
            assert!(
                NewickParser::parse_str(bad, &ParseOptions::default()).is_err(),
                "Expected error for {:?}",
                bad
            );
        }
    }
}
