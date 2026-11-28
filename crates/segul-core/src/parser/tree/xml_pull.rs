//! Zero-dependency streaming XML pull tokenizer and writer helpers for PhyloXML and NeXML.

use std::io::{BufRead, Write};

use anyhow::{Context, Result, bail};

/// Streaming XML token event with namespace-stripped local element names.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum XmlEvent {
    Start {
        name: String,
        attrs: Vec<(String, String)>,
        empty: bool,
    },
    End(String),
    Text(String),
    Eof,
}

/// Streaming pull-based XML tokenizer backed by a reusable byte scratch buffer.
pub(crate) struct XmlPullReader<R: BufRead> {
    reader: R,
    scratch: Vec<u8>,
    bom_checked: bool,
}

impl<R: BufRead> XmlPullReader<R> {
    pub(crate) fn new(reader: R) -> Self {
        Self {
            reader,
            scratch: Vec::with_capacity(512),
            bom_checked: false,
        }
    }

    /// Read the next XML event, skipping processing instructions, comments, and DOCTYPE.
    pub(crate) fn next_event(&mut self) -> Result<XmlEvent> {
        if !self.bom_checked {
            self.skip_utf8_bom()?;
            self.bom_checked = true;
        }

        loop {
            let Some(first) = self.peek_byte()? else {
                return Ok(XmlEvent::Eof);
            };

            if first != b'<' {
                let text = self.read_text_until_lt()?;
                return Ok(XmlEvent::Text(text));
            }

            self.consume_one(); // '<'
            let Some(second) = self.next_byte()? else {
                bail!("Unexpected EOF after '<' in XML stream");
            };

            match second {
                b'?' => self.skip_processing_instruction()?,
                b'!' => {
                    if let Some(text) = self.handle_bang_construct()? {
                        return Ok(XmlEvent::Text(text));
                    }
                }
                b'/' => return self.read_end_tag(),
                _ => return self.read_start_tag(second),
            }
        }
    }

    /// Skip everything inside the current non-empty start element until its matching end tag.
    pub(crate) fn skip_element(&mut self) -> Result<()> {
        let mut depth = 1usize;
        while depth > 0 {
            match self.next_event()? {
                XmlEvent::Start { empty: false, .. } => depth += 1,
                XmlEvent::End(_) => depth -= 1,
                XmlEvent::Eof => bail!("Unexpected EOF while skipping XML element"),
                _ => {}
            }
        }
        Ok(())
    }

    /// Read concatenated text content of the current element until its matching end tag,
    /// skipping any unexpected nested child elements.
    pub(crate) fn read_element_text(&mut self) -> Result<String> {
        let mut text = String::new();
        loop {
            match self.next_event()? {
                XmlEvent::Text(t) => text.push_str(&t),
                XmlEvent::Start { empty: false, .. } => self.skip_element()?,
                XmlEvent::Start { empty: true, .. } => {}
                XmlEvent::End(_) => return Ok(text),
                XmlEvent::Eof => bail!("Unexpected EOF inside XML element text"),
            }
        }
    }

    fn skip_utf8_bom(&mut self) -> Result<()> {
        if self.peek_byte()? == Some(0xEF) {
            self.consume_one();
            if self.next_byte()? != Some(0xBB) || self.next_byte()? != Some(0xBF) {
                bail!("Invalid UTF-8 BOM at start of XML stream");
            }
        }
        Ok(())
    }

    fn peek_byte(&mut self) -> Result<Option<u8>> {
        let buf = self.reader.fill_buf()?;
        Ok(buf.first().copied())
    }

    fn next_byte(&mut self) -> Result<Option<u8>> {
        let buf = self.reader.fill_buf()?;
        let Some(&b) = buf.first() else {
            return Ok(None);
        };
        self.reader.consume(1);
        Ok(Some(b))
    }

    #[inline]
    fn consume_one(&mut self) {
        self.reader.consume(1);
    }

    fn read_text_until_lt(&mut self) -> Result<String> {
        self.scratch.clear();
        loop {
            let buf = self.reader.fill_buf()?;
            if buf.is_empty() {
                break;
            }
            if let Some(pos) = memchr::memchr(b'<', buf) {
                self.scratch.extend_from_slice(&buf[..pos]);
                self.reader.consume(pos);
                break;
            }
            self.scratch.extend_from_slice(buf);
            let len = buf.len();
            self.reader.consume(len);
        }
        unescape_xml(&self.scratch)
    }

    fn skip_processing_instruction(&mut self) -> Result<()> {
        let mut prev_q = false;
        while let Some(b) = self.next_byte()? {
            if prev_q && b == b'>' {
                return Ok(());
            }
            prev_q = b == b'?';
        }
        bail!("Unclosed XML processing instruction '<?'")
    }

    /// Returns `Some(text)` for CDATA sections, `None` for comments and DOCTYPE.
    fn handle_bang_construct(&mut self) -> Result<Option<String>> {
        let Some(b) = self.next_byte()? else {
            bail!("Unexpected EOF after '<!' in XML stream");
        };
        if b == b'-' {
            if self.next_byte()? != Some(b'-') {
                bail!("Malformed XML comment declaration");
            }
            self.skip_comment()?;
            return Ok(None);
        }
        if b == b'[' {
            for expected in b"CDATA[" {
                if self.next_byte()? != Some(*expected) {
                    bail!("Unsupported XML '<![' construct (expected CDATA)");
                }
            }
            let cdata = self.read_cdata()?;
            return Ok(Some(cdata));
        }
        self.skip_doctype()?;
        Ok(None)
    }

    fn skip_comment(&mut self) -> Result<()> {
        let mut dashes = 0usize;
        while let Some(b) = self.next_byte()? {
            if b == b'-' {
                dashes += 1;
            } else if b == b'>' && dashes >= 2 {
                return Ok(());
            } else {
                dashes = 0;
            }
        }
        bail!("Unclosed XML comment '<!--'")
    }

    fn read_cdata(&mut self) -> Result<String> {
        self.scratch.clear();
        let mut brackets = 0usize;
        while let Some(b) = self.next_byte()? {
            if b == b']' {
                brackets += 1;
                self.scratch.push(b);
            } else if b == b'>' && brackets >= 2 {
                let new_len = self.scratch.len() - 2;
                self.scratch.truncate(new_len);
                let text = std::str::from_utf8(&self.scratch)
                    .context("Invalid UTF-8 inside XML CDATA")?
                    .to_string();
                return Ok(text);
            } else {
                brackets = 0;
                self.scratch.push(b);
            }
        }
        bail!("Unclosed XML CDATA section")
    }

    fn skip_doctype(&mut self) -> Result<()> {
        let mut in_single = false;
        let mut in_double = false;
        let mut bracket_depth = 0usize;
        while let Some(b) = self.next_byte()? {
            if b == b'\'' && !in_double {
                in_single = !in_single;
            } else if b == b'"' && !in_single {
                in_double = !in_double;
            } else if !in_single && !in_double {
                if b == b'[' {
                    bracket_depth += 1;
                } else if b == b']' && bracket_depth > 0 {
                    bracket_depth -= 1;
                } else if b == b'>' && bracket_depth == 0 {
                    return Ok(());
                }
            }
        }
        bail!("Unclosed XML DOCTYPE declaration")
    }

    fn read_end_tag(&mut self) -> Result<XmlEvent> {
        self.skip_ws()?;
        let raw_name = self.read_name()?;
        if raw_name.is_empty() {
            bail!("Empty XML end tag name");
        }
        self.skip_ws()?;
        if self.next_byte()? != Some(b'>') {
            bail!("Expected '>' closing XML end tag '</{}'", raw_name);
        }
        Ok(XmlEvent::End(local_name(&raw_name).to_string()))
    }

    fn read_start_tag(&mut self, first_byte: u8) -> Result<XmlEvent> {
        self.scratch.clear();
        if first_byte.is_ascii_whitespace() || matches!(first_byte, b'>' | b'/' | b'=') {
            bail!("Invalid start of XML tag name '{}'", first_byte as char);
        }
        self.scratch.push(first_byte);
        self.read_name_rest()?;
        let raw_name = std::str::from_utf8(&self.scratch)
            .context("Invalid UTF-8 in XML tag name")?
            .to_string();
        let name = local_name(&raw_name).to_string();

        let mut attrs = Vec::new();
        loop {
            self.skip_ws()?;
            let Some(b) = self.peek_byte()? else {
                bail!("Unexpected EOF inside XML start tag '<{}'", raw_name);
            };
            if b == b'>' {
                self.consume_one();
                return Ok(XmlEvent::Start {
                    name,
                    attrs,
                    empty: false,
                });
            }
            if b == b'/' {
                self.consume_one();
                if self.next_byte()? != Some(b'>') {
                    bail!("Expected '>' after '/' in XML empty tag '<{}'", raw_name);
                }
                return Ok(XmlEvent::Start {
                    name,
                    attrs,
                    empty: true,
                });
            }

            let attr_key = self.read_name()?;
            if attr_key.is_empty() {
                bail!("Invalid XML attribute in tag '<{}'", raw_name);
            }
            self.skip_ws()?;
            if self.next_byte()? != Some(b'=') {
                bail!(
                    "Expected '=' after attribute '{}' in '<{}'",
                    attr_key,
                    raw_name
                );
            }
            self.skip_ws()?;
            let Some(quote) = self.next_byte()? else {
                bail!("Unexpected EOF before attribute value for '{}'", attr_key);
            };
            if quote != b'"' && quote != b'\'' {
                bail!("Expected quoted attribute value for '{}'", attr_key);
            }
            let attr_val = self.read_attr_value(quote)?;
            attrs.push((attr_key, attr_val));
        }
    }

    fn skip_ws(&mut self) -> Result<()> {
        loop {
            let buf = self.reader.fill_buf()?;
            if buf.is_empty() {
                return Ok(());
            }
            let count = buf.iter().take_while(|b| b.is_ascii_whitespace()).count();
            if count == 0 {
                return Ok(());
            }
            let reached_end = count == buf.len();
            self.reader.consume(count);
            if !reached_end {
                return Ok(());
            }
        }
    }

    fn read_name(&mut self) -> Result<String> {
        self.scratch.clear();
        self.read_name_rest()?;
        let s = std::str::from_utf8(&self.scratch).context("Invalid UTF-8 in XML name")?;
        Ok(s.to_string())
    }

    fn read_name_rest(&mut self) -> Result<()> {
        loop {
            let buf = self.reader.fill_buf()?;
            if buf.is_empty() {
                return Ok(());
            }
            let count = buf
                .iter()
                .take_while(|&&b| !b.is_ascii_whitespace() && !matches!(b, b'>' | b'/' | b'='))
                .count();
            self.scratch.extend_from_slice(&buf[..count]);
            let hit_end = count < buf.len();
            self.reader.consume(count);
            if hit_end {
                return Ok(());
            }
        }
    }

    fn read_attr_value(&mut self, quote: u8) -> Result<String> {
        self.scratch.clear();
        loop {
            let buf = self.reader.fill_buf()?;
            if buf.is_empty() {
                bail!("Unclosed XML attribute value");
            }
            if let Some(pos) = memchr::memchr(quote, buf) {
                self.scratch.extend_from_slice(&buf[..pos]);
                self.reader.consume(pos + 1);
                return unescape_xml(&self.scratch);
            }
            self.scratch.extend_from_slice(buf);
            let len = buf.len();
            self.reader.consume(len);
        }
    }
}

/// Strip an optional `prefix:` XML namespace qualifier.
#[inline]
pub(crate) fn local_name(raw: &str) -> &str {
    match raw.split_once(':') {
        Some((_, local)) => local,
        None => raw,
    }
}

/// Look up an attribute value by exact key or namespace-stripped local key.
pub(crate) fn attr_value<'a>(attrs: &'a [(String, String)], key: &str) -> Option<&'a str> {
    attrs
        .iter()
        .find(|(k, _)| k == key || local_name(k) == key)
        .map(|(_, v)| v.as_str())
}

/// Unescape standard XML 1.0 entities and numeric character references.
pub(crate) fn unescape_xml(raw: &[u8]) -> Result<String> {
    let text = std::str::from_utf8(raw).context("Invalid UTF-8 in XML text")?;
    if !text.contains('&') {
        return Ok(text.to_string());
    }
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(amp) = rest.find('&') {
        out.push_str(&rest[..amp]);
        let after_amp = &rest[amp + 1..];
        let Some(semi) = after_amp.find(';') else {
            bail!(
                "Unclosed XML entity '&{}'",
                &after_amp[..after_amp.len().min(16)]
            );
        };
        let entity = &after_amp[..semi];
        match entity {
            "amp" => out.push('&'),
            "lt" => out.push('<'),
            "gt" => out.push('>'),
            "quot" => out.push('"'),
            "apos" => out.push('\''),
            _ if entity.starts_with("#x") || entity.starts_with("#X") => {
                let code = u32::from_str_radix(&entity[2..], 16)
                    .with_context(|| format!("Invalid XML hex entity '&{};'", entity))?;
                let ch = char::from_u32(code)
                    .with_context(|| format!("Invalid Unicode code point in '&{};'", entity))?;
                out.push(ch);
            }
            _ if entity.starts_with('#') => {
                let code: u32 = entity[1..]
                    .parse()
                    .with_context(|| format!("Invalid XML numeric entity '&{};'", entity))?;
                let ch = char::from_u32(code)
                    .with_context(|| format!("Invalid Unicode code point in '&{};'", entity))?;
                out.push(ch);
            }
            other => bail!("Unsupported XML entity '&{};'", other),
        }
        rest = &after_amp[semi + 1..];
    }
    out.push_str(rest);
    Ok(out)
}

/// Write XML-escaped text (`&`, `<`, `>`, `"`, `'`) to `writer`.
pub(crate) fn write_xml_escaped<W: Write>(writer: &mut W, text: &str) -> Result<()> {
    let bytes = text.as_bytes();
    let mut start = 0usize;
    for (i, &b) in bytes.iter().enumerate() {
        let esc = match b {
            b'&' => "&amp;",
            b'<' => "&lt;",
            b'>' => "&gt;",
            b'"' => "&quot;",
            b'\'' => "&apos;",
            _ => continue,
        };
        if i > start {
            writer.write_all(&bytes[start..i])?;
        }
        writer.write_all(esc.as_bytes())?;
        start = i + 1;
    }
    if start < bytes.len() {
        writer.write_all(&bytes[start..])?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_xml_pull_entities_comments_cdata_and_skip() {
        let xml = b"\xEF\xBB\xBF<?xml version=\"1.0\"?>\n\
                   <!DOCTYPE nexml [ <!ENTITY foo \"bar\"> ]>\n\
                   <!-- leading comment -->\n\
                   <nex:root attr='a &amp; b' empty=\"&#x41;\">\n\
                     <skip><nested/></skip>\n\
                     <keep><![CDATA[1 < 2 & 3]]></keep>\n\
                   </nex:root>";
        let mut reader = XmlPullReader::new(Cursor::new(xml.as_slice()));
        let next_non_ws = |r: &mut XmlPullReader<_>| loop {
            match r.next_event().unwrap() {
                XmlEvent::Text(t) if t.trim().is_empty() => continue,
                ev => return ev,
            }
        };
        match next_non_ws(&mut reader) {
            XmlEvent::Start { name, attrs, empty } => {
                assert_eq!(name, "root");
                assert!(!empty);
                assert_eq!(attr_value(&attrs, "attr"), Some("a & b"));
                assert_eq!(attr_value(&attrs, "empty"), Some("A"));
            }
            other => panic!("Unexpected {:?}", other),
        }
        match next_non_ws(&mut reader) {
            XmlEvent::Start { name, .. } => {
                assert_eq!(name, "skip");
                reader.skip_element().unwrap();
            }
            other => panic!("Unexpected {:?}", other),
        }
        match next_non_ws(&mut reader) {
            XmlEvent::Start { name, .. } => {
                assert_eq!(name, "keep");
                assert_eq!(reader.read_element_text().unwrap(), "1 < 2 & 3");
            }
            other => panic!("Unexpected {:?}", other),
        }
    }

    #[test]
    fn test_xml_malformed_errors() {
        for bad in [
            "<root",
            "<root attr=unquoted/>",
            "<root><!-- unclosed",
            "<r>&bogus;</r>",
        ] {
            let mut reader = XmlPullReader::new(Cursor::new(bad.as_bytes()));
            let mut errored = false;
            loop {
                match reader.next_event() {
                    Ok(XmlEvent::Eof) => break,
                    Ok(_) => {}
                    Err(_) => {
                        errored = true;
                        break;
                    }
                }
            }
            assert!(errored, "Expected error for {:?}", bad);
        }
    }
}
