//! Streaming NEXUS `TREES` block parser with zero impact on sequence NEXUS parsing.

use std::io::BufRead;
use std::sync::Arc;

use ahash::AHashMap;
use anyhow::{Result, bail};
use rayon::prelude::*;

use super::model::{Rootedness, TaxonId, TaxonTable, Tree};
use super::newick::{NewickParser, ParseOptions};

const MAX_DIRECT_TRANSLATE_KEY: usize = 1_000_000;

/// Mapping from NEXUS `TRANSLATE` keys (numeric or alphanumeric) to interned [`TaxonId`]s.
#[derive(Debug, Clone)]
pub struct TranslateTable {
    taxa: Arc<TaxonTable>,
    numeric: Vec<Option<TaxonId>>,
    named: AHashMap<String, TaxonId>,
}

impl TranslateTable {
    /// Build a translation table from `(key, taxon_label)` pairs.
    pub fn from_pairs(pairs: Vec<(String, String)>) -> Self {
        let mut table = TaxonTable::with_capacity(pairs.len());
        let mut numeric = Vec::new();
        let mut named = AHashMap::new();

        for (key, label) in pairs {
            let taxon_id = table.intern(&label);
            if let Some(num_key) = parse_compact_positive_key(&key) {
                if num_key >= numeric.len() {
                    numeric.resize(num_key + 1, None);
                }
                numeric[num_key] = Some(taxon_id);
            } else {
                named.insert(key, taxon_id);
            }
        }

        Self {
            taxa: Arc::new(table),
            numeric,
            named,
        }
    }

    /// Build a 1-based translation table from `BEGIN TAXA; TAXLABELS ...;`.
    pub fn from_taxlabels(labels: &[String]) -> Self {
        let mut table = TaxonTable::with_capacity(labels.len());
        let mut numeric = Vec::with_capacity(labels.len() + 1);
        numeric.push(None);
        for label in labels {
            numeric.push(Some(table.intern(label)));
        }
        Self {
            taxa: Arc::new(table),
            numeric,
            named: AHashMap::new(),
        }
    }

    /// Return the shared [`TaxonTable`] backing this translation mapping.
    #[inline]
    pub fn taxa(&self) -> Arc<TaxonTable> {
        Arc::clone(&self.taxa)
    }

    /// Look up a leaf token in O(1) for compact integer keys or via hash map otherwise.
    #[inline]
    pub fn lookup(&self, key: &str) -> Option<TaxonId> {
        if let Some(num) = parse_compact_positive_key(key)
            && let Some(Some(id)) = self.numeric.get(num)
        {
            return Some(*id);
        }
        self.named.get(key).copied()
    }
}

#[inline]
fn parse_compact_positive_key(key: &str) -> Option<usize> {
    if key.is_empty() || !key.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let val: usize = key.parse().ok()?;
    if (1..=MAX_DIRECT_TRANSLATE_KEY).contains(&val) {
        Some(val)
    } else {
        None
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TokenKind {
    Word,
    Quoted,
    Punct,
}

/// Low-level streaming byte/token scanner over a `BufRead` NEXUS stream.
struct NexusScanner<R: BufRead> {
    reader: R,
    eof: bool,
}

impl<R: BufRead> NexusScanner<R> {
    fn new(reader: R) -> Self {
        Self { reader, eof: false }
    }

    fn verify_header(&mut self) -> Result<()> {
        let mut tok = Vec::with_capacity(16);
        let Some(kind) = self.next_token(&mut tok, true)? else {
            bail!("Empty file: expected #NEXUS header");
        };
        if kind != TokenKind::Word || !tok.eq_ignore_ascii_case(b"#nexus") {
            let preview = String::from_utf8_lossy(&tok);
            bail!("Invalid NEXUS file: expected '#NEXUS', found '{}'", preview);
        }
        Ok(())
    }

    /// Read the next NEXUS token, skipping whitespace and nested `[...]` comments.
    fn next_token(
        &mut self,
        out: &mut Vec<u8>,
        allow_double_quotes: bool,
    ) -> Result<Option<TokenKind>> {
        out.clear();
        self.skip_ws_and_comments()?;
        let available = self.reader.fill_buf()?;
        if available.is_empty() {
            self.eof = true;
            return Ok(None);
        }

        let first = available[0];
        if matches!(first, b';' | b'=' | b',' | b'*') {
            out.push(first);
            self.reader.consume(1);
            return Ok(Some(TokenKind::Punct));
        }

        if first == b'\'' || (allow_double_quotes && first == b'"') {
            self.read_quoted_token(first, out)?;
            return Ok(Some(TokenKind::Quoted));
        }

        loop {
            let buf = self.reader.fill_buf()?;
            if buf.is_empty() {
                self.eof = true;
                break;
            }
            let mut end = 0usize;
            let mut hit_delim = false;
            while end < buf.len() {
                let b = buf[end];
                if b.is_ascii_whitespace() || matches!(b, b';' | b'=' | b',' | b'*' | b'[') {
                    hit_delim = true;
                    break;
                }
                end += 1;
            }
            out.extend_from_slice(&buf[..end]);
            self.reader.consume(end);
            if hit_delim {
                break;
            }
        }
        Ok(Some(TokenKind::Word))
    }

    fn read_quoted_token(&mut self, quote: u8, out: &mut Vec<u8>) -> Result<()> {
        self.reader.consume(1); // opening quote
        loop {
            let buf = self.reader.fill_buf()?;
            if buf.is_empty() {
                bail!("Unclosed quoted token in NEXUS file");
            }
            if let Some(rel) = memchr::memchr(quote, buf) {
                out.extend_from_slice(&buf[..rel]);
                self.reader.consume(rel + 1);
                let next_buf = self.reader.fill_buf()?;
                if next_buf.first() == Some(&quote) {
                    out.push(quote);
                    self.reader.consume(1);
                    continue;
                }
                return Ok(());
            }
            out.extend_from_slice(buf);
            let len = buf.len();
            self.reader.consume(len);
        }
    }

    fn skip_ws_and_comments(&mut self) -> Result<()> {
        loop {
            let buf = self.reader.fill_buf()?;
            if buf.is_empty() {
                self.eof = true;
                return Ok(());
            }
            let mut consumed = 0usize;
            let mut hit_comment = false;
            while consumed < buf.len() {
                let b = buf[consumed];
                if b == b'[' {
                    hit_comment = true;
                    break;
                }
                if !b.is_ascii_whitespace() && !matches!(b, 0xEF | 0xBB | 0xBF) {
                    self.reader.consume(consumed);
                    return Ok(());
                }
                consumed += 1;
            }
            self.reader.consume(consumed);
            if hit_comment {
                self.skip_bracket_comment()?;
            }
        }
    }

    fn skip_bracket_comment(&mut self) -> Result<()> {
        self.reader.consume(1); // consume '['
        let mut depth = 1usize;
        let mut in_single = false;
        let mut in_double = false;

        while depth > 0 {
            let buf = self.reader.fill_buf()?;
            if buf.is_empty() {
                bail!("Unclosed '[' comment in NEXUS file");
            }
            let mut i = 0usize;
            while i < buf.len() {
                let b = buf[i];
                if b == b'\'' && !in_double {
                    in_single = !in_single;
                } else if b == b'"' && !in_single {
                    in_double = !in_double;
                } else if !in_single && !in_double {
                    if b == b'[' {
                        depth += 1;
                    } else if b == b']' {
                        depth -= 1;
                        if depth == 0 {
                            self.reader.consume(i + 1);
                            return Ok(());
                        }
                    }
                }
                i += 1;
            }
            self.reader.consume(i);
        }
        Ok(())
    }

    /// Skip directly to the next `;` outside quotes and nested `[...]` comments.
    fn skip_until_semicolon(&mut self) -> Result<()> {
        let mut in_single = false;
        let mut in_double = false;
        let mut comment_depth = 0usize;

        loop {
            let buf = self.reader.fill_buf()?;
            if buf.is_empty() {
                self.eof = true;
                return Ok(());
            }
            for (i, &b) in buf.iter().enumerate() {
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
                } else if in_single {
                    if b == b'\'' {
                        in_single = false;
                    }
                } else if in_double {
                    if b == b'"' {
                        in_double = false;
                    }
                } else {
                    match b {
                        b'\'' => in_single = true,
                        b'"' => in_double = true,
                        b'[' => comment_depth = 1,
                        b';' => {
                            self.reader.consume(i + 1);
                            return Ok(());
                        }
                        _ => {}
                    }
                }
            }
            let len = buf.len();
            self.reader.consume(len);
        }
    }

    /// Read a Newick tree body (preserving comments and quotes) up to and including `;`.
    fn read_tree_body(&mut self, out: &mut Vec<u8>, allow_double_quotes: bool) -> Result<()> {
        out.clear();
        let mut in_single = false;
        let mut in_double = false;
        let mut comment_depth = 0usize;

        loop {
            let buf = self.reader.fill_buf()?;
            if buf.is_empty() {
                bail!("Unexpected end of NEXUS TREE statement before ';'");
            }
            let mut consumed = 0usize;
            let mut done = false;
            for (i, &b) in buf.iter().enumerate() {
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
                } else if in_single {
                    if b == b'\'' {
                        in_single = false;
                    }
                } else if in_double {
                    if b == b'"' {
                        in_double = false;
                    }
                } else {
                    match b {
                        b'\'' => in_single = true,
                        b'"' if allow_double_quotes => in_double = true,
                        b'[' => comment_depth = 1,
                        b';' => {
                            consumed = i + 1;
                            done = true;
                            break;
                        }
                        _ => {}
                    }
                }
                consumed = i + 1;
            }
            out.extend_from_slice(&buf[..consumed]);
            self.reader.consume(consumed);
            if done {
                return Ok(());
            }
        }
    }
}

/// Streaming NEXUS `TREES` block reader.
///
/// Skips `DATA`, `CHARACTERS`, `SETS`, `ASSUMPTIONS`, `MRBAYES`, and any other
/// non-`TREES` blocks without buffering alignment matrices, and supports
/// `BEGIN TAXA; TAXLABELS ...;` fallback when `TRANSLATE` is omitted.
pub struct NexusTreeReader<R: BufRead> {
    scanner: NexusScanner<R>,
    opts: ParseOptions,
    header_checked: bool,
    in_trees_block: bool,
    taxa_fallback: Option<Arc<TranslateTable>>,
    block_translate: Option<Arc<TranslateTable>>,
    token_buf: Vec<u8>,
    tree_buf: Vec<u8>,
    eof: bool,
}

impl<R: BufRead> NexusTreeReader<R> {
    /// Create a streaming NEXUS tree reader.
    pub fn new(reader: R, opts: ParseOptions) -> Self {
        Self {
            scanner: NexusScanner::new(reader),
            opts,
            header_checked: false,
            in_trees_block: false,
            taxa_fallback: None,
            block_translate: None,
            token_buf: Vec::with_capacity(256),
            tree_buf: Vec::with_capacity(4096),
            eof: false,
        }
    }

    /// Skip up to `count` trees without parsing Newick topologies (fast MCMC burn-in skip).
    pub fn skip_trees(&mut self, count: usize) -> Result<usize> {
        let mut skipped = 0usize;
        while skipped < count {
            if self.advance_to_next_tree()?.is_none() {
                break;
            }
            self.scanner.skip_until_semicolon()?;
            skipped += 1;
        }
        Ok(skipped)
    }

    fn ensure_header(&mut self) -> Result<()> {
        if !self.header_checked {
            self.scanner.verify_header()?;
            self.header_checked = true;
        }
        Ok(())
    }

    fn advance_to_next_tree(&mut self) -> Result<Option<(String, bool, Rootedness)>> {
        if self.eof {
            return Ok(None);
        }
        self.ensure_header()?;

        let mut tok = Vec::with_capacity(128);
        loop {
            if !self.in_trees_block {
                let Some(kind) = self
                    .scanner
                    .next_token(&mut tok, self.opts.allow_double_quotes)?
                else {
                    self.eof = true;
                    return Ok(None);
                };
                if kind == TokenKind::Punct {
                    continue;
                }
                if tok.eq_ignore_ascii_case(b"begin") {
                    let Some(_) = self
                        .scanner
                        .next_token(&mut tok, self.opts.allow_double_quotes)?
                    else {
                        bail!("Unexpected EOF after 'BEGIN' in NEXUS file");
                    };
                    if tok.eq_ignore_ascii_case(b"trees") {
                        self.scanner.skip_until_semicolon()?;
                        self.in_trees_block = true;
                        self.block_translate = None;
                    } else if tok.eq_ignore_ascii_case(b"taxa") {
                        self.scanner.skip_until_semicolon()?;
                        self.parse_taxa_block()?;
                    } else {
                        self.scanner.skip_until_semicolon()?;
                        self.skip_non_trees_block()?;
                    }
                } else {
                    self.scanner.skip_until_semicolon()?;
                }
                continue;
            }

            // Inside BEGIN TREES;
            let Some(kind) = self
                .scanner
                .next_token(&mut tok, self.opts.allow_double_quotes)?
            else {
                self.eof = true;
                return Ok(None);
            };
            if kind == TokenKind::Punct {
                continue;
            }
            if tok.eq_ignore_ascii_case(b"end") || tok.eq_ignore_ascii_case(b"endblock") {
                self.scanner.skip_until_semicolon()?;
                self.in_trees_block = false;
                continue;
            }
            if tok.eq_ignore_ascii_case(b"translate") {
                self.parse_translate_command()?;
                continue;
            }
            if tok.eq_ignore_ascii_case(b"tree")
                || tok.eq_ignore_ascii_case(b"utree")
                || tok.eq_ignore_ascii_case(b"rtree")
            {
                let default_rooted = if tok.eq_ignore_ascii_case(b"utree") {
                    Rootedness::Unrooted
                } else if tok.eq_ignore_ascii_case(b"rtree") {
                    Rootedness::Rooted
                } else {
                    Rootedness::Unknown
                };
                let (name, is_default) = self.parse_tree_header()?;
                return Ok(Some((name, is_default, default_rooted)));
            }

            // Skip TITLE, LINK, DIMENSIONS, ID, etc.
            self.scanner.skip_until_semicolon()?;
        }
    }

    fn skip_non_trees_block(&mut self) -> Result<()> {
        let mut tok = Vec::with_capacity(64);
        loop {
            let Some(kind) = self.scanner.next_token(&mut tok, true)? else {
                self.eof = true;
                return Ok(());
            };
            if kind == TokenKind::Punct {
                continue;
            }
            if tok.eq_ignore_ascii_case(b"end") || tok.eq_ignore_ascii_case(b"endblock") {
                self.scanner.skip_until_semicolon()?;
                return Ok(());
            }
            self.scanner.skip_until_semicolon()?;
        }
    }

    fn parse_taxa_block(&mut self) -> Result<()> {
        let mut tok = Vec::with_capacity(64);
        loop {
            let Some(kind) = self
                .scanner
                .next_token(&mut tok, self.opts.allow_double_quotes)?
            else {
                self.eof = true;
                return Ok(());
            };
            if kind == TokenKind::Punct {
                continue;
            }
            if tok.eq_ignore_ascii_case(b"end") || tok.eq_ignore_ascii_case(b"endblock") {
                self.scanner.skip_until_semicolon()?;
                return Ok(());
            }
            if tok.eq_ignore_ascii_case(b"taxlabels") {
                let mut labels = Vec::new();
                while let Some(tk) = self
                    .scanner
                    .next_token(&mut tok, self.opts.allow_double_quotes)?
                {
                    if tk == TokenKind::Punct && tok == b";" {
                        break;
                    }
                    let text = decode_nexus_label(
                        &tok,
                        tk == TokenKind::Quoted,
                        self.opts.underscores_to_spaces,
                    );
                    labels.push(text);
                }
                if !labels.is_empty() {
                    self.taxa_fallback = Some(Arc::new(TranslateTable::from_taxlabels(&labels)));
                }
            } else {
                self.scanner.skip_until_semicolon()?;
            }
        }
    }

    fn parse_translate_command(&mut self) -> Result<()> {
        let mut pairs = Vec::new();
        let mut key_buf = Vec::with_capacity(32);
        let mut val_buf = Vec::with_capacity(64);

        loop {
            let Some(key_kind) = self
                .scanner
                .next_token(&mut key_buf, self.opts.allow_double_quotes)?
            else {
                bail!("Unexpected EOF inside NEXUS TRANSLATE command");
            };
            if key_kind == TokenKind::Punct {
                if key_buf == b";" {
                    break;
                }
                if key_buf == b"," {
                    continue;
                }
                bail!(
                    "Unexpected '{}' in NEXUS TRANSLATE command",
                    key_buf[0] as char
                );
            }

            let Some(val_kind) = self
                .scanner
                .next_token(&mut val_buf, self.opts.allow_double_quotes)?
            else {
                bail!("Unexpected EOF after TRANSLATE key");
            };
            if val_kind == TokenKind::Punct {
                bail!("Missing taxon label in NEXUS TRANSLATE command");
            }

            let key = String::from_utf8_lossy(&key_buf).into_owned();
            let label = decode_nexus_label(
                &val_buf,
                val_kind == TokenKind::Quoted,
                self.opts.underscores_to_spaces,
            );
            pairs.push((key, label));
        }

        self.block_translate = Some(Arc::new(TranslateTable::from_pairs(pairs)));
        Ok(())
    }

    fn parse_tree_header(&mut self) -> Result<(String, bool)> {
        let mut is_default = false;
        let Some(mut kind) = self
            .scanner
            .next_token(&mut self.token_buf, self.opts.allow_double_quotes)?
        else {
            bail!("Unexpected EOF in NEXUS TREE statement");
        };

        if kind == TokenKind::Punct && self.token_buf == b"*" {
            is_default = true;
            let Some(next_kind) = self
                .scanner
                .next_token(&mut self.token_buf, self.opts.allow_double_quotes)?
            else {
                bail!("Unexpected EOF after '*' in NEXUS TREE statement");
            };
            kind = next_kind;
        }

        if kind == TokenKind::Punct {
            bail!(
                "Missing tree name in NEXUS TREE statement (found '{}')",
                self.token_buf[0] as char
            );
        }
        let name = String::from_utf8_lossy(&self.token_buf).into_owned();

        let Some(eq_kind) = self
            .scanner
            .next_token(&mut self.token_buf, self.opts.allow_double_quotes)?
        else {
            bail!("Unexpected EOF waiting for '=' after tree '{}'", name);
        };
        if eq_kind != TokenKind::Punct || self.token_buf != b"=" {
            bail!(
                "Expected '=' after tree name '{}', found '{}'",
                name,
                String::from_utf8_lossy(&self.token_buf)
            );
        }

        Ok((name, is_default))
    }

    fn read_next_tree(&mut self) -> Result<Option<Tree>> {
        let Some((name, is_default, default_rooted)) = self.advance_to_next_tree()? else {
            return Ok(None);
        };

        self.scanner
            .read_tree_body(&mut self.tree_buf, self.opts.allow_double_quotes)?;
        let table = self
            .block_translate
            .as_deref()
            .or(self.taxa_fallback.as_deref());
        let mut tree = NewickParser::parse_bytes_with_translate(&self.tree_buf, table, &self.opts)?;
        tree.metadata.name = Some(name);
        tree.metadata.is_default = is_default;
        if tree.metadata.rootedness == Rootedness::Unknown {
            tree.metadata.rootedness = default_rooted;
        }
        Ok(Some(tree))
    }
}

impl<R: BufRead> Iterator for NexusTreeReader<R> {
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

/// Parse all trees in an in-memory NEXUS byte buffer in parallel using `rayon`.
///
/// `TAXA` and `TRANSLATE` tables are parsed once sequentially and shared
/// read-only across worker threads while Newick bodies are parsed in parallel.
pub fn parse_nexus_bytes_par(input: &[u8], opts: &ParseOptions) -> Result<Vec<Tree>> {
    struct RawNexusTree {
        name: String,
        is_default: bool,
        default_rooted: Rootedness,
        translate: Option<Arc<TranslateTable>>,
        body: Vec<u8>,
    }

    let mut reader = NexusTreeReader::new(std::io::Cursor::new(input), opts.clone());
    let mut raw_trees = Vec::new();
    while let Some((name, is_default, default_rooted)) = reader.advance_to_next_tree()? {
        let mut body = Vec::new();
        reader
            .scanner
            .read_tree_body(&mut body, opts.allow_double_quotes)?;
        let translate = reader
            .block_translate
            .clone()
            .or_else(|| reader.taxa_fallback.clone());
        raw_trees.push(RawNexusTree {
            name,
            is_default,
            default_rooted,
            translate,
            body,
        });
    }

    raw_trees
        .into_par_iter()
        .map(|raw| {
            let mut tree = NewickParser::parse_bytes_with_translate(
                &raw.body,
                raw.translate.as_deref(),
                opts,
            )?;
            tree.metadata.name = Some(raw.name);
            tree.metadata.is_default = raw.is_default;
            if tree.metadata.rootedness == Rootedness::Unknown {
                tree.metadata.rootedness = raw.default_rooted;
            }
            Ok(tree)
        })
        .collect()
}

fn decode_nexus_label(bytes: &[u8], quoted: bool, underscores_to_spaces: bool) -> String {
    let text = String::from_utf8_lossy(bytes);
    if !quoted && underscores_to_spaces && bytes.contains(&b'_') {
        text.replace('_', " ")
    } else {
        text.into_owned()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_nexus_translate_and_multiple_trees() {
        let nex = "#NEXUS\nBEGIN TREES;\n\
                   TRANSLATE\n  1 'Homo sapiens',\n  2 Pan_troglodytes,\n  3 Gorilla_gorilla;\n\
                   TREE * t1 = [&R] ((1:0.1,2:0.1):0.2,3:0.3);\n\
                   UTREE t2 = (1:0.2,2:0.2,3:0.4);\n\
                   END;\n";
        let trees: Vec<Tree> =
            NexusTreeReader::new(Cursor::new(nex.as_bytes()), ParseOptions::default())
                .collect::<Result<_>>()
                .unwrap();
        assert_eq!(trees.len(), 2);
        assert_eq!(trees[0].metadata.name.as_deref(), Some("t1"));
        assert!(trees[0].metadata.is_default);
        assert_eq!(trees[0].metadata.rootedness, Rootedness::Rooted);
        assert_eq!(trees[1].metadata.name.as_deref(), Some("t2"));
        assert_eq!(trees[1].metadata.rootedness, Rootedness::Unrooted);

        let leaves0: Vec<&str> = trees[0]
            .leaves()
            .filter_map(|id| trees[0].label(id))
            .collect();
        assert_eq!(
            leaves0,
            vec!["Homo sapiens", "Pan_troglodytes", "Gorilla_gorilla"]
        );
        assert!(Arc::ptr_eq(
            trees[0].taxa().unwrap(),
            trees[1].taxa().unwrap()
        ));
    }

    #[test]
    fn test_nexus_taxa_fallback_and_block_skipping() {
        let nex = "#NEXUS\n\
                   [Top comment [nested]]\n\
                   BEGIN TAXA;\n  DIMENSIONS NTAX=3;\n  TAXLABELS Alpha Beta Gamma;\nEND;\n\
                   BEGIN DATA;\n  DIMENSIONS NTAX=3 NCHAR=4;\n  FORMAT DATATYPE=DNA;\n  MATRIX\n    Alpha ACGT\n    Beta ACGT\n    Gamma ACGT;\nEND;\n\
                   BEGIN TREES;\n  TREE t1 = ((1:0.1,2:0.2):0.3,3:0.4);\nEND;\n";
        let trees: Vec<Tree> =
            NexusTreeReader::new(Cursor::new(nex.as_bytes()), ParseOptions::default())
                .collect::<Result<_>>()
                .unwrap();
        assert_eq!(trees.len(), 1);
        let leaves: Vec<&str> = trees[0]
            .leaves()
            .filter_map(|id| trees[0].label(id))
            .collect();
        assert_eq!(leaves, vec!["Alpha", "Beta", "Gamma"]);
    }

    #[test]
    fn test_nexus_skip_burnin_and_parallel() {
        let nex = "#NEXUS\nBEGIN TREES;\n\
                   TRANSLATE 1 A, 2 B, 3 C;\n\
                   TREE gen_0 = (1,2,3);\n\
                   TREE gen_100 = (1,(2,3));\n\
                   TREE gen_200 = ((1,2),3);\n\
                   END;\n";
        let mut reader = NexusTreeReader::new(Cursor::new(nex.as_bytes()), ParseOptions::default());
        assert_eq!(reader.skip_trees(2).unwrap(), 2);
        let rest: Vec<Tree> = reader.collect::<Result<_>>().unwrap();
        assert_eq!(rest.len(), 1);
        assert_eq!(rest[0].metadata.name.as_deref(), Some("gen_200"));

        let par = parse_nexus_bytes_par(nex.as_bytes(), &ParseOptions::default()).unwrap();
        assert_eq!(par.len(), 3);
        assert_eq!(par[1].metadata.name.as_deref(), Some("gen_100"));
    }

    #[test]
    fn test_invalid_nexus_header_errors() {
        let mut reader = NexusTreeReader::new(Cursor::new(b"(A,B);"), ParseOptions::default());
        assert!(reader.next().unwrap().is_err());
    }
}
