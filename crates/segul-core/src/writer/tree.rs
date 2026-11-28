//! Iterative phylogenetic tree writers for Newick, Extended Newick, NHX, and NEXUS.

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

use ahash::AHashMap;
use anyhow::{Context, Result};
use flate2::{Compression, write::GzEncoder};
use indexmap::IndexSet;

use crate::helper::files;
use crate::parser::IO_BUF_SIZE;
use crate::parser::tree::xml_pull::write_xml_escaped;
use crate::parser::tree::{
    Annotation, AnnotationFormat, AnnotationValue, NodeId, Rootedness, Tree, TreeOutputFmt,
};

/// Configuration options for writing phylogenetic trees.
#[derive(Debug, Clone)]
pub struct TreeWriteOptions {
    /// Emit node and branch annotations (`[&...]`, `[&&NHX:...]`, `[raw]`). Default: `true`.
    pub write_annotations: bool,
    /// Emit a NEXUS `TRANSLATE` table when writing NEXUS output. Default: `true`.
    pub use_translate: bool,
    /// Emit an explicit `BEGIN TAXA;` block before `BEGIN TREES;` in NEXUS output. Default: `false`.
    pub include_taxa_block: bool,
    /// Omit branch lengths from output. Default: `false`.
    pub strip_branch_lengths: bool,
    /// Emit Extended Newick `:len:support:gamma` fields when present. Default: `true`.
    pub write_extended_fields: bool,
    /// Compress output with gzip. Default: `false`.
    pub gzip: bool,
}

impl Default for TreeWriteOptions {
    fn default() -> Self {
        Self {
            write_annotations: true,
            use_translate: true,
            include_taxa_block: false,
            strip_branch_lengths: false,
            write_extended_fields: true,
            gzip: false,
        }
    }
}

/// Serialize a single [`Tree`] to a Newick/Extended Newick string using an iterative traversal.
pub fn to_newick_string(tree: &Tree, opts: &TreeWriteOptions) -> String {
    let est = tree.len().saturating_mul(16).max(32);
    let mut buf = Vec::with_capacity(est);
    write_newick_tree(&mut buf, tree, None, opts).expect("Writing to Vec<u8> cannot fail");
    String::from_utf8(buf).expect("Tree Newick output is valid UTF-8")
}

/// Write a single [`Tree`] in Newick/Extended Newick format to `writer` iteratively.
///
/// When `translate` is `Some`, leaf labels present in the map are replaced by
/// their 1-based NEXUS integer keys.
pub fn write_newick_tree<W: Write>(
    writer: &mut W,
    tree: &Tree,
    translate: Option<&AHashMap<&str, usize>>,
    opts: &TreeWriteOptions,
) -> Result<()> {
    if tree.is_empty() {
        return Ok(());
    }

    enum Frame {
        Open(NodeId),
        Comma,
        Close(NodeId),
    }

    let mut stack = Vec::with_capacity(64);
    let mut children_buf = Vec::new();
    stack.push(Frame::Open(tree.root()));

    while let Some(frame) = stack.pop() {
        match frame {
            Frame::Open(id) => {
                if tree.is_leaf(id) {
                    write_node_suffix(writer, tree, id, true, translate, opts)?;
                } else {
                    writer.write_all(b"(")?;
                    stack.push(Frame::Close(id));
                    children_buf.clear();
                    children_buf.extend(tree.children(id));
                    for idx in (0..children_buf.len()).rev() {
                        stack.push(Frame::Open(children_buf[idx]));
                        if idx > 0 {
                            stack.push(Frame::Comma);
                        }
                    }
                }
            }
            Frame::Comma => writer.write_all(b",")?,
            Frame::Close(id) => {
                writer.write_all(b")")?;
                write_node_suffix(writer, tree, id, false, translate, opts)?;
            }
        }
    }

    writer.write_all(b";")?;
    Ok(())
}

/// Write a slice of [`Tree`]s in NEXUS format (`#NEXUS` + `BEGIN TREES; ... END;`).
pub fn write_nexus_trees<W: Write>(
    writer: &mut W,
    trees: &[Tree],
    opts: &TreeWriteOptions,
) -> Result<()> {
    writer.write_all(b"#NEXUS\n")?;

    let mut ordered_taxa: IndexSet<&str> = IndexSet::new();
    if opts.use_translate || opts.include_taxa_block {
        for tree in trees {
            for leaf in tree.leaves() {
                if let Some(lbl) = tree.label(leaf) {
                    ordered_taxa.insert(lbl);
                }
            }
        }
    }

    if opts.include_taxa_block && !ordered_taxa.is_empty() {
        writeln!(writer, "BEGIN TAXA;")?;
        writeln!(writer, "  DIMENSIONS NTAX={};", ordered_taxa.len())?;
        write!(writer, "  TAXLABELS")?;
        for taxon in &ordered_taxa {
            writer.write_all(b" ")?;
            write_escaped_label(writer, taxon, false)?;
        }
        writeln!(writer, ";")?;
        writeln!(writer, "END;")?;
    }

    writeln!(writer, "BEGIN TREES;")?;

    let translate_map: Option<AHashMap<&str, usize>> =
        if opts.use_translate && !ordered_taxa.is_empty() {
            writeln!(writer, "  TRANSLATE")?;
            let mut map = AHashMap::with_capacity(ordered_taxa.len());
            let last_idx = ordered_taxa.len() - 1;
            for (idx, taxon) in ordered_taxa.iter().enumerate() {
                let key = idx + 1;
                map.insert(*taxon, key);
                write!(writer, "    {} ", key)?;
                write_escaped_label(writer, taxon, false)?;
                if idx < last_idx {
                    writer.write_all(b",\n")?;
                } else {
                    writer.write_all(b"\n")?;
                }
            }
            writeln!(writer, "  ;")?;
            Some(map)
        } else {
            None
        };

    for (idx, tree) in trees.iter().enumerate() {
        write_nexus_tree_line(writer, tree, idx, translate_map.as_ref(), opts)?;
    }

    writeln!(writer, "END;")?;
    Ok(())
}

/// File writer for phylogenetic trees supporting Newick, Extended Newick, NEXUS, and gzip.
pub struct TreeWriter<'a> {
    output: &'a Path,
    format: TreeOutputFmt,
    opts: TreeWriteOptions,
}

impl<'a> TreeWriter<'a> {
    /// Create a new writer targeting `output` in `format`.
    pub fn new(output: &'a Path, format: TreeOutputFmt, opts: TreeWriteOptions) -> Self {
        Self {
            output,
            format,
            opts,
        }
    }

    /// Write an in-memory slice of trees to `output`.
    pub fn write_trees(&self, trees: &[Tree]) -> Result<()> {
        let mut writer = self.open_writer()?;
        let effective_opts = self.effective_opts();
        match self.format {
            TreeOutputFmt::Newick | TreeOutputFmt::ExtendedNewick => {
                for tree in trees {
                    write_newick_tree(&mut writer, tree, None, &effective_opts)?;
                    writer.write_all(b"\n")?;
                }
            }
            TreeOutputFmt::Nexus => {
                write_nexus_trees(&mut writer, trees, &effective_opts)?;
            }
            TreeOutputFmt::PhyloXml => {
                write_phyloxml_trees(&mut writer, trees, &effective_opts)?;
            }
            TreeOutputFmt::NeXml => {
                write_nexml_trees(&mut writer, trees, &effective_opts)?;
            }
        }
        writer.flush()?;
        Ok(())
    }

    /// Stream trees from an iterator to `output`, returning the number of trees written.
    ///
    /// For Newick/Extended Newick, PhyloXML, and NEXUS with `use_translate == false`,
    /// writes one tree at a time in constant memory.
    pub fn write_stream<I>(&self, trees: I) -> Result<usize>
    where
        I: Iterator<Item = Result<Tree>>,
    {
        let effective_opts = self.effective_opts();
        match self.format {
            TreeOutputFmt::Newick | TreeOutputFmt::ExtendedNewick => {
                let mut writer = self.open_writer()?;
                let mut count = 0usize;
                for item in trees {
                    let tree = item?;
                    write_newick_tree(&mut writer, &tree, None, &effective_opts)?;
                    writer.write_all(b"\n")?;
                    count += 1;
                }
                writer.flush()?;
                Ok(count)
            }
            TreeOutputFmt::PhyloXml => {
                let mut writer = self.open_writer()?;
                write_phyloxml_header(&mut writer)?;
                let mut count = 0usize;
                for item in trees {
                    let tree = item?;
                    write_phyloxml_phylogeny(&mut writer, &tree, &effective_opts)?;
                    count += 1;
                }
                write_phyloxml_footer(&mut writer)?;
                writer.flush()?;
                Ok(count)
            }
            TreeOutputFmt::Nexus
                if !effective_opts.use_translate && !effective_opts.include_taxa_block =>
            {
                let mut writer = self.open_writer()?;
                writer.write_all(b"#NEXUS\nBEGIN TREES;\n")?;
                let mut count = 0usize;
                for item in trees {
                    let tree = item?;
                    write_nexus_tree_line(&mut writer, &tree, count, None, &effective_opts)?;
                    count += 1;
                }
                writer.write_all(b"END;\n")?;
                writer.flush()?;
                Ok(count)
            }
            TreeOutputFmt::Nexus | TreeOutputFmt::NeXml => {
                let collected: Vec<Tree> = trees.collect::<Result<_>>()?;
                let count = collected.len();
                self.write_trees(&collected)?;
                Ok(count)
            }
        }
    }

    fn effective_opts(&self) -> TreeWriteOptions {
        let mut opts = self.opts.clone();
        if self.format == TreeOutputFmt::Newick {
            opts.write_extended_fields = false;
        }
        opts
    }

    fn open_writer(&self) -> Result<Box<dyn Write>> {
        if let Some(parent) = self.output.parent()
            && !parent.as_os_str().is_empty()
        {
            fs::create_dir_all(parent).with_context(|| {
                format!(
                    "Failed creating output directory for {}",
                    self.output.display()
                )
            })?;
        }
        let file = File::create(self.output)
            .with_context(|| format!("Failed creating {}", self.output.display()))?;
        let use_gz = self.opts.gzip || files::is_gzip(self.output);
        if use_gz {
            let enc = GzEncoder::new(file, Compression::default());
            Ok(Box::new(BufWriter::with_capacity(IO_BUF_SIZE, enc)))
        } else {
            Ok(Box::new(BufWriter::with_capacity(IO_BUF_SIZE, file)))
        }
    }
}

fn write_nexus_tree_line<W: Write>(
    writer: &mut W,
    tree: &Tree,
    idx: usize,
    translate: Option<&AHashMap<&str, usize>>,
    opts: &TreeWriteOptions,
) -> Result<()> {
    let fallback_name;
    let name = match tree.metadata.name.as_deref() {
        Some(n) if !n.is_empty() => n,
        _ => {
            fallback_name = format!("tree_{}", idx + 1);
            &fallback_name
        }
    };
    writer.write_all(b"  TREE ")?;
    if tree.metadata.is_default {
        writer.write_all(b"* ")?;
    }
    write_escaped_label(writer, name, false)?;
    writer.write_all(b" = ")?;
    match tree.metadata.rootedness {
        Rootedness::Rooted => writer.write_all(b"[&R] ")?,
        Rootedness::Unrooted => writer.write_all(b"[&U] ")?,
        Rootedness::Unknown => {}
    }
    if let Some(w) = tree.metadata.weight {
        write!(writer, "[&W {w}] ")?;
    }
    write_newick_tree(writer, tree, translate, opts)?;
    writer.write_all(b"\n")?;
    Ok(())
}

fn write_node_suffix<W: Write>(
    writer: &mut W,
    tree: &Tree,
    id: NodeId,
    is_leaf: bool,
    translate: Option<&AHashMap<&str, usize>>,
    opts: &TreeWriteOptions,
) -> Result<()> {
    if let Some(lbl) = tree.label(id) {
        if is_leaf
            && let Some(map) = translate
            && let Some(&key) = map.get(lbl)
        {
            write!(writer, "{key}")?;
        } else {
            write_escaped_label(writer, lbl, tree.is_quoted_label(id))?;
        }
    }

    if opts.write_annotations {
        for annot in tree.node_annotations(id) {
            write_annotation(writer, annot)?;
        }
    }

    let len = if opts.strip_branch_lengths {
        None
    } else {
        tree.branch_length(id)
    };
    let hybrid = if opts.write_extended_fields {
        tree.hybrid_info(id)
    } else {
        None
    };
    let sup = hybrid.and_then(|h| h.edge_support);
    let gamma = hybrid.and_then(|h| h.inheritance_prob);

    if len.is_some() || sup.is_some() || gamma.is_some() {
        writer.write_all(b":")?;
        if let Some(v) = len {
            write!(writer, "{v}")?;
        }
        if sup.is_some() || gamma.is_some() {
            writer.write_all(b":")?;
            if let Some(s) = sup {
                write!(writer, "{s}")?;
            }
            if let Some(g) = gamma {
                write!(writer, ":{g}")?;
            }
        }
    }

    if opts.write_annotations {
        for annot in tree.branch_annotations(id) {
            write_annotation(writer, annot)?;
        }
    }

    Ok(())
}

fn write_annotation<W: Write>(writer: &mut W, annot: &Annotation) -> Result<()> {
    writer.write_all(b"[")?;
    match annot.format {
        AnnotationFormat::Nhx | AnnotationFormat::Beast | AnnotationFormat::Raw => {
            writer.write_all(annot.raw.as_bytes())?;
        }
    }
    writer.write_all(b"]")?;
    Ok(())
}

fn write_escaped_label<W: Write>(writer: &mut W, label: &str, was_quoted: bool) -> Result<()> {
    let needs_quotes = label.is_empty()
        || label.bytes().any(|b| {
            b.is_ascii_whitespace()
                || b.is_ascii_control()
                || matches!(
                    b,
                    b'(' | b')' | b',' | b':' | b';' | b'[' | b']' | b'\'' | b'"'
                )
        });

    if needs_quotes || (was_quoted && label.is_empty()) {
        writer.write_all(b"'")?;
        if label.contains('\'') {
            for b in label.bytes() {
                if b == b'\'' {
                    writer.write_all(b"''")?;
                } else {
                    writer.write_all(&[b])?;
                }
            }
        } else {
            writer.write_all(label.as_bytes())?;
        }
        writer.write_all(b"'")?;
    } else {
        writer.write_all(label.as_bytes())?;
    }
    Ok(())
}

/// Write a slice of [`Tree`]s in PhyloXML format iteratively.
pub fn write_phyloxml_trees<W: Write>(
    writer: &mut W,
    trees: &[Tree],
    opts: &TreeWriteOptions,
) -> Result<()> {
    write_phyloxml_header(writer)?;
    for tree in trees {
        write_phyloxml_phylogeny(writer, tree, opts)?;
    }
    write_phyloxml_footer(writer)
}

fn write_phyloxml_header<W: Write>(writer: &mut W) -> Result<()> {
    writer.write_all(
        b"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
          <phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns=\"http://www.phyloxml.org\">\n",
    )?;
    Ok(())
}

fn write_phyloxml_footer<W: Write>(writer: &mut W) -> Result<()> {
    writer.write_all(b"</phyloxml>\n")?;
    Ok(())
}

/// Write a single [`Tree`] as a PhyloXML `<phylogeny>` element using an explicit stack.
pub fn write_phyloxml_phylogeny<W: Write>(
    writer: &mut W,
    tree: &Tree,
    opts: &TreeWriteOptions,
) -> Result<()> {
    match tree.metadata.rootedness {
        Rootedness::Rooted => writer.write_all(b"  <phylogeny rooted=\"true\">\n")?,
        Rootedness::Unrooted => writer.write_all(b"  <phylogeny rooted=\"false\">\n")?,
        Rootedness::Unknown => writer.write_all(b"  <phylogeny>\n")?,
    }

    if let Some(name) = tree.metadata.name.as_deref()
        && !name.is_empty()
    {
        writer.write_all(b"    <name>")?;
        write_xml_escaped(writer, name)?;
        writer.write_all(b"</name>\n")?;
    }

    if tree.is_empty() {
        writer.write_all(b"  </phylogeny>\n")?;
        return Ok(());
    }

    enum PhyloFrame {
        Open(NodeId, usize),
        Close(usize),
    }

    let mut stack = Vec::with_capacity(64);
    let mut children_buf = Vec::new();
    stack.push(PhyloFrame::Open(tree.root(), 2));

    while let Some(frame) = stack.pop() {
        match frame {
            PhyloFrame::Open(id, depth) => {
                write_indent(writer, depth)?;
                writer.write_all(b"<clade")?;
                if !opts.strip_branch_lengths
                    && let Some(len) = tree.branch_length(id)
                {
                    write!(writer, " branch_length=\"{len}\"")?;
                }

                let label = tree.label(id);
                let node_annots = if opts.write_annotations {
                    tree.node_annotations(id)
                } else {
                    &[]
                };
                let branch_annots = if opts.write_annotations {
                    tree.branch_annotations(id)
                } else {
                    &[]
                };
                let is_leaf = tree.is_leaf(id);
                let has_body = label.is_some()
                    || !is_leaf
                    || !node_annots.is_empty()
                    || !branch_annots.is_empty();

                if !has_body {
                    writer.write_all(b"/>\n")?;
                    continue;
                }

                writer.write_all(b">\n")?;
                if let Some(lbl) = label {
                    write_indent(writer, depth + 1)?;
                    writer.write_all(b"<name>")?;
                    write_xml_escaped(writer, lbl)?;
                    writer.write_all(b"</name>\n")?;
                }
                for annot in node_annots {
                    write_phyloxml_annotation(writer, annot, "clade", depth + 1)?;
                }
                for annot in branch_annots {
                    write_phyloxml_annotation(writer, annot, "parent_branch", depth + 1)?;
                }

                stack.push(PhyloFrame::Close(depth));
                if !is_leaf {
                    children_buf.clear();
                    children_buf.extend(tree.children(id));
                    for &child in children_buf.iter().rev() {
                        stack.push(PhyloFrame::Open(child, depth + 1));
                    }
                }
            }
            PhyloFrame::Close(depth) => {
                write_indent(writer, depth)?;
                writer.write_all(b"</clade>\n")?;
            }
        }
    }

    writer.write_all(b"  </phylogeny>\n")?;
    Ok(())
}

fn write_phyloxml_annotation<W: Write>(
    writer: &mut W,
    annot: &Annotation,
    applies_to: &str,
    depth: usize,
) -> Result<()> {
    match annot.format {
        AnnotationFormat::Raw => {
            if let Ok(val) = annot.raw.trim().parse::<f64>()
                && val.is_finite()
            {
                write_indent(writer, depth)?;
                writeln!(writer, "<confidence type=\"support\">{val}</confidence>")?;
            } else {
                write_indent(writer, depth)?;
                write!(
                    writer,
                    "<property ref=\"comment\" datatype=\"xsd:string\" applies_to=\"{applies_to}\">"
                )?;
                write_xml_escaped(writer, &annot.raw)?;
                writer.write_all(b"</property>\n")?;
            }
        }
        AnnotationFormat::Beast | AnnotationFormat::Nhx => {
            if let Some((_, conf_val)) = annot
                .entries
                .iter()
                .find(|(k, _)| k.eq_ignore_ascii_case("confidence"))
            {
                let conf_type = annot
                    .entries
                    .iter()
                    .find(|(k, _)| k.eq_ignore_ascii_case("confidence_type"))
                    .map(|(_, v)| v.to_string());
                write_indent(writer, depth)?;
                if let Some(kind) = conf_type {
                    writer.write_all(b"<confidence type=\"")?;
                    write_xml_escaped(writer, &kind)?;
                    writer.write_all(b"\">")?;
                } else {
                    writer.write_all(b"<confidence>")?;
                }
                write_xml_escaped(writer, &conf_val.to_string())?;
                writer.write_all(b"</confidence>\n")?;
            }
            for (k, v) in &annot.entries {
                if k.eq_ignore_ascii_case("confidence") || k.eq_ignore_ascii_case("confidence_type")
                {
                    continue;
                }
                write_indent(writer, depth)?;
                let dtype = xml_xsd_type(v);
                writer.write_all(b"<property ref=\"")?;
                write_xml_escaped(writer, k)?;
                write!(
                    writer,
                    "\" datatype=\"{dtype}\" applies_to=\"{applies_to}\">"
                )?;
                write_xml_escaped(writer, &v.to_string())?;
                writer.write_all(b"</property>\n")?;
            }
        }
    }
    Ok(())
}

/// Write a slice of [`Tree`]s in NeXML format iteratively.
pub fn write_nexml_trees<W: Write>(
    writer: &mut W,
    trees: &[Tree],
    opts: &TreeWriteOptions,
) -> Result<()> {
    writer.write_all(
        b"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
          <nex:nexml version=\"0.9\" xmlns:nex=\"http://www.nexml.org/2009\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n",
    )?;

    let mut ordered_taxa: IndexSet<&str> = IndexSet::new();
    for tree in trees {
        for leaf in tree.leaves() {
            if let Some(lbl) = tree.label(leaf) {
                ordered_taxa.insert(lbl);
            }
        }
    }

    writer.write_all(b"  <otus id=\"otus1\">\n")?;
    let mut otu_map: AHashMap<&str, String> = AHashMap::with_capacity(ordered_taxa.len());
    for (idx, taxon) in ordered_taxa.iter().enumerate() {
        let otu_id = format!("otu{}", idx + 1);
        writer.write_all(b"    <otu id=\"")?;
        writer.write_all(otu_id.as_bytes())?;
        writer.write_all(b"\" label=\"")?;
        write_xml_escaped(writer, taxon)?;
        writer.write_all(b"\"/>\n")?;
        otu_map.insert(*taxon, otu_id);
    }
    writer.write_all(b"  </otus>\n")?;

    writer.write_all(b"  <trees id=\"trees1\" otus=\"otus1\">\n")?;
    for (t_idx, tree) in trees.iter().enumerate() {
        let elem = if tree.is_network() { "network" } else { "tree" };
        let tree_id = format!("tree{}", t_idx + 1);
        write!(writer, "    <{elem} id=\"{tree_id}\"")?;
        if elem == "tree" {
            writer.write_all(b" xsi:type=\"nex:FloatTree\"")?;
        }
        if let Some(name) = tree.metadata.name.as_deref()
            && !name.is_empty()
        {
            writer.write_all(b" label=\"")?;
            write_xml_escaped(writer, name)?;
            writer.write_all(b"\"")?;
        }
        writer.write_all(b">\n")?;

        if !tree.is_empty() {
            let root = tree.root();
            if !opts.strip_branch_lengths
                && let Some(root_len) = tree.branch_length(root)
            {
                writeln!(
                    writer,
                    "      <rootedge id=\"re{}\" target=\"n{}\" length=\"{root_len}\"/>",
                    t_idx + 1,
                    root.index()
                )?;
            }

            for id in tree.preorder() {
                write!(writer, "      <node id=\"n{}\"", id.index())?;
                if tree.is_root(id) {
                    writer.write_all(b" root=\"true\"")?;
                }
                if let Some(lbl) = tree.label(id) {
                    if tree.is_leaf(id)
                        && let Some(otu_id) = otu_map.get(lbl)
                    {
                        write!(writer, " otu=\"{otu_id}\"")?;
                    } else {
                        writer.write_all(b" label=\"")?;
                        write_xml_escaped(writer, lbl)?;
                        writer.write_all(b"\"")?;
                    }
                }
                let annots = if opts.write_annotations {
                    tree.node_annotations(id)
                } else {
                    &[]
                };
                if annots.is_empty() {
                    writer.write_all(b"/>\n")?;
                } else {
                    writer.write_all(b">\n")?;
                    write_nexml_meta_list(writer, annots)?;
                    writer.write_all(b"      </node>\n")?;
                }
            }

            let mut edge_idx = 1usize;
            for id in tree.preorder() {
                let Some(parent) = tree.parent(id) else {
                    continue;
                };
                write!(
                    writer,
                    "      <edge id=\"e{edge_idx}\" source=\"n{}\" target=\"n{}\"",
                    parent.index(),
                    id.index()
                )?;
                edge_idx += 1;
                if !opts.strip_branch_lengths
                    && let Some(len) = tree.branch_length(id)
                {
                    write!(writer, " length=\"{len}\"")?;
                }
                let annots = if opts.write_annotations {
                    tree.branch_annotations(id)
                } else {
                    &[]
                };
                if annots.is_empty() {
                    writer.write_all(b"/>\n")?;
                } else {
                    writer.write_all(b">\n")?;
                    write_nexml_meta_list(writer, annots)?;
                    writer.write_all(b"      </edge>\n")?;
                }
            }
        }

        writeln!(writer, "    </{elem}>")?;
    }
    writer.write_all(b"  </trees>\n</nex:nexml>\n")?;
    Ok(())
}

fn write_nexml_meta_list<W: Write>(writer: &mut W, annots: &[Annotation]) -> Result<()> {
    for annot in annots {
        if annot.entries.is_empty() {
            writer.write_all(
                b"        <meta xsi:type=\"nex:LiteralMeta\" property=\"comment\" content=\"",
            )?;
            write_xml_escaped(writer, &annot.raw)?;
            writer.write_all(b"\"/>\n")?;
        } else {
            for (k, v) in &annot.entries {
                writer.write_all(b"        <meta xsi:type=\"nex:LiteralMeta\" property=\"")?;
                write_xml_escaped(writer, k)?;
                writer.write_all(b"\" content=\"")?;
                write_xml_escaped(writer, &v.to_string())?;
                writer.write_all(b"\"/>\n")?;
            }
        }
    }
    Ok(())
}

fn write_indent<W: Write>(writer: &mut W, depth: usize) -> Result<()> {
    let spaces = depth.min(16) * 2;
    const PAD: &[u8] = b"                                ";
    writer.write_all(&PAD[..spaces])?;
    Ok(())
}

fn xml_xsd_type(val: &AnnotationValue) -> &'static str {
    match val {
        AnnotationValue::Int(_) => "xsd:integer",
        AnnotationValue::Float(_) => "xsd:decimal",
        AnnotationValue::Bool(_) => "xsd:boolean",
        AnnotationValue::Str(_) | AnnotationValue::Array(_) => "xsd:string",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::tree::{NewickParser, ParseOptions, TreeFmt, parse_trees_file};

    #[test]
    fn test_round_trip_newick_and_enewick() {
        let input = "(('Homo sapiens':0.1,'Canis ''lupus''':0.2)95:0.3,((A:0.1,#H1:0.05::0.3):0.1,(#H1:0.15:95:0.7,B:0.2):0.1):0.4);";
        let tree = NewickParser::parse_str(input, &ParseOptions::default()).unwrap();
        let written = to_newick_string(&tree, &TreeWriteOptions::default());
        assert_eq!(written, input);
    }

    #[test]
    fn test_round_trip_deep_pectinate_100k_no_stack_overflow() {
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
        let written = to_newick_string(&tree, &TreeWriteOptions::default());
        assert_eq!(written, s);
    }

    #[test]
    fn test_write_nexus_and_gzip_round_trip() {
        let tree1 = NewickParser::parse_str("((A:0.1,B:0.2):0.3,C:0.4);", &ParseOptions::default())
            .unwrap();
        let tree2 = NewickParser::parse_str("(A:0.2,(B:0.1,C:0.3):0.4);", &ParseOptions::default())
            .unwrap();

        let dir = tempfile::TempDir::with_prefix("segul_tree_write").unwrap();
        let nex_gz = dir.path().join("out.nex.gz");
        let writer = TreeWriter::new(&nex_gz, TreeOutputFmt::Nexus, TreeWriteOptions::default());
        writer.write_trees(&[tree1, tree2]).unwrap();

        let reparsed = parse_trees_file(&nex_gz, TreeFmt::Auto).unwrap();
        assert_eq!(reparsed.len(), 2);
        assert_eq!(
            reparsed[0]
                .leaves()
                .filter_map(|id| reparsed[0].label(id))
                .collect::<Vec<_>>(),
            vec!["A", "B", "C"]
        );
    }

    #[test]
    fn test_round_trip_phyloxml_and_nexml() {
        let input = "(('A & B':0.1,C:0.2)95:0.3,D:0.4);";
        let mut tree = NewickParser::parse_str(input, &ParseOptions::default()).unwrap();
        tree.metadata.name = Some("demo".to_string());
        tree.metadata.rootedness = Rootedness::Rooted;

        let dir = tempfile::TempDir::with_prefix("segul_xml_rt").unwrap();
        let pxml = dir.path().join("demo.phyloxml");
        let nxml = dir.path().join("demo.nexml");

        TreeWriter::new(&pxml, TreeOutputFmt::PhyloXml, TreeWriteOptions::default())
            .write_trees(&[tree.clone()])
            .unwrap();
        TreeWriter::new(&nxml, TreeOutputFmt::NeXml, TreeWriteOptions::default())
            .write_trees(&[tree])
            .unwrap();

        let from_pxml = parse_trees_file(&pxml, TreeFmt::Auto).unwrap();
        let from_nxml = parse_trees_file(&nxml, TreeFmt::Auto).unwrap();
        assert_eq!(from_pxml.len(), 1);
        assert_eq!(from_nxml.len(), 1);
        assert_eq!(
            from_pxml[0]
                .leaves()
                .filter_map(|id| from_pxml[0].label(id))
                .collect::<Vec<_>>(),
            vec!["A & B", "C", "D"]
        );
        assert_eq!(
            from_nxml[0]
                .leaves()
                .filter_map(|id| from_nxml[0].label(id))
                .collect::<Vec<_>>(),
            vec!["A & B", "C", "D"]
        );
    }

    #[test]
    fn test_round_trip_deep_phyloxml_50k_no_stack_overflow() {
        let depth = 50_000usize;
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
        let mut xml_buf = Vec::new();
        write_phyloxml_trees(&mut xml_buf, &[tree], &TreeWriteOptions::default()).unwrap();

        let dir = tempfile::TempDir::with_prefix("segul_pxml_deep").unwrap();
        let pxml = dir.path().join("deep.phyloxml");
        std::fs::write(&pxml, &xml_buf).unwrap();
        let reparsed = parse_trees_file(&pxml, TreeFmt::PhyloXml).unwrap();
        assert_eq!(reparsed.len(), 1);
        assert_eq!(reparsed[0].tip_count(), depth + 1);
    }
}
