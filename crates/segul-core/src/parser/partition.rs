//! Parse RAxML and NEXUS partition files, including merged, strided, and charpartition schemes.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, bail};
use indexmap::IndexMap;

use crate::helper::types::{Partition, PartitionFmt, SiteRange};

const DISALLOWED_NAME_CHARS: &str = r#"()/\,"';:?!<>|*"#;

/// Parse a partition file and return a vector of partitions.
pub struct PartitionParser<'a> {
    path: &'a Path,
    partition_fmt: &'a PartitionFmt,
    is_uncheck: bool,
}

impl<'a> PartitionParser<'a> {
    pub fn new(path: &'a Path, partition_fmt: &'a PartitionFmt, is_uncheck: bool) -> Self {
        Self {
            path,
            partition_fmt,
            is_uncheck,
        }
    }

    /// Parse the partition file, panicking with a descriptive message on failure.
    pub fn parse(&self) -> Vec<Partition> {
        self.try_parse().unwrap_or_else(|err| panic!("{err}"))
    }

    /// Parse and validate the partition file, returning an error on invalid input.
    pub fn try_parse(&self) -> Result<Vec<Partition>> {
        let file = File::open(self.path)
            .with_context(|| format!("Unable to open partition file {}", self.path.display()))?;
        let mut reader = BufReader::new(file);
        match self.partition_fmt {
            PartitionFmt::Nexus
            | PartitionFmt::NexusCodon
            | PartitionFmt::Charset
            | PartitionFmt::CharsetCodon => self.parse_nexus(&mut reader),
            PartitionFmt::Raxml | PartitionFmt::RaxmlCodon => self.parse_raxml(&mut reader),
        }
    }

    /// Parse NEXUS partition text from `reader` instead of the file at `path`.
    /// Used to read only the tail of a large NEXUS alignment.
    pub(crate) fn try_parse_nexus_from<R: BufRead>(
        &self,
        reader: &mut R,
    ) -> Result<Vec<Partition>> {
        self.parse_nexus(reader)
    }

    /// Parse a RAxML / RAxML-NG / PartitionFinder RAxML-style partition file.
    fn parse_raxml<R: BufRead>(&self, reader: &mut R) -> Result<Vec<Partition>> {
        let mut partitions = Vec::new();
        for (idx, line_res) in reader.lines().enumerate() {
            let line = line_res?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with("//") {
                continue;
            }
            let line_num = idx + 1;
            let (lhs, rhs) = trimmed.split_once('=').with_context(|| {
                format!("Invalid RAxML partition at line {line_num}: missing '=' in '{trimmed}'")
            })?;
            let (gene, model) = if let Some((prefix, name)) = lhs.split_once(',') {
                let raw_model = prefix.trim();
                let model = if !raw_model.is_empty() && !is_generic_raxml_dna(raw_model) {
                    Some(raw_model.to_string())
                } else {
                    None
                };
                (sanitize_gene_name(name)?, model)
            } else {
                (sanitize_gene_name(lhs)?, None)
            };
            let ranges = parse_site_ranges(rhs).with_context(|| {
                format!("Invalid site range at line {line_num} for partition '{gene}'")
            })?;
            partitions.push(Partition::from_ranges(gene, ranges, model));
        }
        self.validate_partitions(&partitions)?;
        Ok(partitions)
    }

    /// Parse a NEXUS `sets` / `assumptions` / `mrbayes` partition block.
    fn parse_nexus<R: BufRead>(&self, reader: &mut R) -> Result<Vec<Partition>> {
        let mut raw = String::new();
        reader.read_to_string(&mut raw)?;
        let cleaned = strip_nexus_comments(&raw);
        let statements = split_nexus_statements(&cleaned);

        let mut charsets: IndexMap<String, Vec<SiteRange>> = IndexMap::new();
        let mut charpartition_raw: Option<String> = None;
        let mut current_block: Option<String> = None;
        let mut seen_sets_block = false;
        let mut sets_block_has_data = false;

        for stmt in statements {
            let trimmed = stmt.trim();
            if trimmed.is_empty() {
                continue;
            }
            if let Some(block) = parse_begin_block(trimmed) {
                current_block = Some(block);
                sets_block_has_data = false;
                continue;
            }
            if is_end_block(trimmed) {
                if matches!(current_block.as_deref(), Some("sets") | Some("assumptions"))
                    && sets_block_has_data
                {
                    seen_sets_block = true;
                }
                current_block = None;
                continue;
            }
            let is_active = match current_block.as_deref() {
                None | Some("sets") | Some("assumptions") => true,
                Some("mrbayes") => !seen_sets_block,
                Some(_) => false,
            };
            if !is_active {
                continue;
            }

            if let Some(rest) = strip_keyword_ci(trimmed, "charset") {
                let (lhs, rhs) = rest
                    .split_once('=')
                    .with_context(|| format!("Invalid NEXUS charset statement: '{trimmed}'"))?;
                let gene = sanitize_gene_name(lhs)?;
                let ranges = parse_nexus_charset_rhs(rhs, &charsets)
                    .with_context(|| format!("Invalid site range for charset '{gene}'"))?;
                charsets.insert(gene, ranges);
                sets_block_has_data = true;
            } else if let Some(rest) = strip_keyword_ci(trimmed, "charpartition")
                && charpartition_raw.is_none()
            {
                let (_, rhs) = rest.split_once('=').with_context(|| {
                    format!("Invalid NEXUS charpartition statement: '{trimmed}'")
                })?;
                charpartition_raw = Some(rhs.to_string());
                sets_block_has_data = true;
            }
        }

        let partitions = if let Some(cp_rhs) = charpartition_raw {
            let resolved = resolve_charpartition(&cp_rhs, &charsets)?;
            if resolved.is_empty() {
                charsets_to_partitions(charsets)
            } else {
                resolved
            }
        } else {
            charsets_to_partitions(charsets)
        };

        self.validate_partitions(&partitions)?;
        Ok(partitions)
    }

    fn validate_partitions(&self, partitions: &[Partition]) -> Result<()> {
        if partitions.is_empty() {
            bail!("Failed parsing partition. No partition found.");
        }
        for part in partitions {
            if part.gene.is_empty() || part.gene.chars().any(char::is_whitespace) {
                bail!("Invalid partition name '{}'.", part.gene);
            }
            let ranges = part.effective_ranges();
            if ranges.is_empty() {
                bail!("Partition '{}' has no site ranges.", part.gene);
            }
            for range in ranges.iter() {
                if !range.is_valid() {
                    bail!(
                        "Invalid site range '{}' in partition '{}'.",
                        range,
                        part.gene
                    );
                }
            }
        }
        if self.is_uncheck {
            return Ok(());
        }

        let mut names = HashSet::with_capacity(partitions.len());
        for part in partitions {
            if !names.insert(part.gene.as_str()) {
                bail!("Duplicate partition name '{}'.", part.gene);
            }
        }

        if partitions.iter().any(Partition::has_stride) {
            validate_strided_coverage(partitions)
        } else {
            validate_contiguous_intervals(partitions)
        }
    }
}

fn charsets_to_partitions(charsets: IndexMap<String, Vec<SiteRange>>) -> Vec<Partition> {
    charsets
        .into_iter()
        .map(|(gene, ranges)| Partition::from_ranges(gene, ranges, None))
        .collect()
}

fn validate_contiguous_intervals(partitions: &[Partition]) -> Result<()> {
    let mut intervals: Vec<(usize, usize, &str)> = Vec::new();
    for part in partitions {
        for range in part.effective_ranges().iter() {
            intervals.push((range.start, range.end, part.gene.as_str()));
        }
    }
    intervals.sort_unstable_by_key(|&(start, end, _)| (start, end));

    let first_start = intervals[0].0;
    if first_start != 1 {
        bail!(
            "Invalid partition input. First partition start position is {} not 1.",
            first_start
        );
    }

    let mut prev_end = intervals[0].1;
    let mut prev_gene = intervals[0].2;
    for &(start, end, gene) in intervals.iter().skip(1) {
        if start <= prev_end {
            bail!(
                "Overlapping partition sites at position {} between '{}' and '{}'.",
                start,
                prev_gene,
                gene
            );
        }
        if start > prev_end + 1 {
            bail!(
                "Invalid partition format. Gap in coverage between position {} and {}.",
                prev_end,
                start
            );
        }
        prev_end = end;
        prev_gene = gene;
    }
    Ok(())
}

fn validate_strided_coverage(partitions: &[Partition]) -> Result<()> {
    let min_start = partitions
        .iter()
        .filter_map(|p| p.effective_ranges().iter().map(|r| r.start).min())
        .min()
        .unwrap_or(0);
    if min_start != 1 {
        bail!(
            "Invalid partition input. First partition start position is {} not 1.",
            min_start
        );
    }
    let max_end = partitions
        .iter()
        .map(Partition::max_site)
        .max()
        .unwrap_or(0);
    let words = max_end.div_ceil(64);
    let mut seen = vec![0u64; words];

    for part in partitions {
        for range in part.effective_ranges().iter() {
            for site in range.iter_sites() {
                let idx = site - 1;
                let mask = 1u64 << (idx & 63);
                if seen[idx >> 6] & mask != 0 {
                    bail!(
                        "Overlapping partition site at position {} in partition '{}'.",
                        site,
                        part.gene
                    );
                }
                seen[idx >> 6] |= mask;
            }
        }
    }

    for (word_idx, &word) in seen.iter().enumerate() {
        let valid_bits = if word_idx + 1 == words && max_end % 64 != 0 {
            (1u64 << (max_end % 64)) - 1
        } else {
            u64::MAX
        };
        if word != valid_bits {
            let missing_bit = (word ^ valid_bits).trailing_zeros() as usize;
            let missing_site = word_idx * 64 + missing_bit + 1;
            bail!(
                "Invalid partition format. Missing site at position {} (gap in coverage).",
                missing_site
            );
        }
    }
    Ok(())
}

fn is_generic_raxml_dna(label: &str) -> bool {
    label.eq_ignore_ascii_case("dna") || label.eq_ignore_ascii_case("dnax")
}

fn sanitize_gene_name(raw: &str) -> Result<String> {
    let mut name = strip_outer_quotes(raw.trim()).trim().to_string();
    if name.chars().any(char::is_whitespace) {
        bail!(
            "Failed parsing partition file. Gene name cannot contain spaces: '{}'",
            name
        );
    }
    name.retain(|c| !DISALLOWED_NAME_CHARS.contains(c));
    if name.contains('.') {
        name = name.replace('.', "_");
    }
    if name.is_empty() {
        bail!("Failed parsing partition file. Partition name is empty.");
    }
    Ok(name)
}

fn strip_outer_quotes(s: &str) -> &str {
    if s.len() >= 2
        && ((s.starts_with('\'') && s.ends_with('\'')) || (s.starts_with('"') && s.ends_with('"')))
    {
        &s[1..s.len() - 1]
    } else {
        s
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum RangeToken {
    Number(usize),
    Dash,
    Backslash,
}

/// Parse one or more site ranges (`1-100`, `1-100, 201-300`, `1-300\3 2-300\3`, `15`).
pub(crate) fn parse_site_ranges(spec: &str) -> Result<Vec<SiteRange>> {
    let tokens = tokenize_site_ranges(spec)?;
    if tokens.is_empty() {
        bail!("Empty site range specification");
    }
    let mut ranges = Vec::new();
    let mut i = 0;
    while i < tokens.len() {
        let start = match tokens[i] {
            RangeToken::Number(n) => n,
            RangeToken::Dash => bail!("Unexpected '-' in site range '{spec}'"),
            RangeToken::Backslash => bail!("Unexpected '\\' in site range '{spec}'"),
        };
        i += 1;

        let end = if i < tokens.len() && tokens[i] == RangeToken::Dash {
            i += 1;
            match tokens.get(i) {
                Some(RangeToken::Number(n)) => {
                    i += 1;
                    *n
                }
                _ => bail!("Expected end position after '-' in site range '{spec}'"),
            }
        } else {
            start
        };

        let step = if i < tokens.len() && tokens[i] == RangeToken::Backslash {
            i += 1;
            match tokens.get(i) {
                Some(RangeToken::Number(n)) => {
                    i += 1;
                    *n
                }
                _ => bail!("Expected stride number after '\\' in site range '{spec}'"),
            }
        } else {
            1
        };

        let range = SiteRange::with_step(start, end, step);
        if !range.is_valid() {
            bail!(
                "Invalid site range '{range}' (start={start}, end={end}, step={step}) in '{spec}'"
            );
        }
        ranges.push(range);
    }
    Ok(ranges)
}

fn tokenize_site_ranges(spec: &str) -> Result<Vec<RangeToken>> {
    let mut tokens = Vec::new();
    let bytes = spec.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let ch = bytes[i] as char;
        if ch.is_ascii_whitespace() || ch == ',' {
            i += 1;
        } else if ch == '-' {
            tokens.push(RangeToken::Dash);
            i += 1;
        } else if ch == '\\' {
            tokens.push(RangeToken::Backslash);
            i += 1;
        } else if ch.is_ascii_digit() {
            let start = i;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                i += 1;
            }
            let num: usize = spec[start..i].parse().with_context(|| {
                format!("Invalid site number '{}' in '{spec}'", &spec[start..i])
            })?;
            tokens.push(RangeToken::Number(num));
        } else {
            bail!("Unexpected character '{ch}' in site range '{spec}'");
        }
    }
    Ok(tokens)
}

fn strip_nexus_comments(input: &str) -> String {
    let mut out = String::with_capacity(input.len());
    let mut in_single = false;
    let mut in_double = false;
    let mut comment_depth: usize = 0;

    for ch in input.chars() {
        if comment_depth == 0 {
            if ch == '\'' && !in_double {
                in_single = !in_single;
                out.push(ch);
            } else if ch == '"' && !in_single {
                in_double = !in_double;
                out.push(ch);
            } else if ch == '[' && !in_single && !in_double {
                comment_depth = 1;
                out.push(' ');
            } else {
                out.push(ch);
            }
        } else if ch == '[' {
            comment_depth += 1;
        } else if ch == ']' {
            comment_depth -= 1;
        }
    }
    out
}

fn split_nexus_statements(input: &str) -> Vec<String> {
    let mut stmts = Vec::new();
    let mut current = String::new();
    let mut in_single = false;
    let mut in_double = false;

    for ch in input.chars() {
        if ch == '\'' && !in_double {
            in_single = !in_single;
            current.push(ch);
        } else if ch == '"' && !in_single {
            in_double = !in_double;
            current.push(ch);
        } else if ch == ';' && !in_single && !in_double {
            let trimmed = current.trim();
            if !trimmed.is_empty() {
                stmts.push(trimmed.to_string());
            }
            current.clear();
        } else {
            current.push(ch);
        }
    }
    let trimmed = current.trim();
    if !trimmed.is_empty() {
        stmts.push(trimmed.to_string());
    }
    stmts
}

fn strip_keyword_ci<'a>(stmt: &'a str, keyword: &str) -> Option<&'a str> {
    if stmt.len() < keyword.len() {
        return None;
    }
    let (head, tail) = stmt.split_at(keyword.len());
    if head.eq_ignore_ascii_case(keyword) && tail.starts_with(char::is_whitespace) {
        Some(tail.trim_start())
    } else {
        None
    }
}

fn parse_begin_block(stmt: &str) -> Option<String> {
    strip_keyword_ci(stmt, "begin").map(|rest| rest.trim().to_ascii_lowercase())
}

fn is_end_block(stmt: &str) -> bool {
    stmt.eq_ignore_ascii_case("end") || stmt.eq_ignore_ascii_case("endblock")
}

fn lookup_charset<'a>(
    charsets: &'a IndexMap<String, Vec<SiteRange>>,
    token: &str,
) -> Option<(&'a String, &'a Vec<SiteRange>)> {
    let unquoted = strip_outer_quotes(token.trim());
    if let Some(found) = charsets.get_key_value(unquoted) {
        return Some(found);
    }
    if let Ok(sanitized) = sanitize_gene_name(unquoted)
        && let Some(found) = charsets.get_key_value(&sanitized)
    {
        return Some(found);
    }
    charsets
        .iter()
        .find(|(key, _)| key.eq_ignore_ascii_case(unquoted))
}

fn split_quote_aware_tokens(input: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    let mut current = String::new();
    let mut in_single = false;
    let mut in_double = false;

    for ch in input.chars() {
        if ch == '\'' && !in_double {
            in_single = !in_single;
            current.push(ch);
        } else if ch == '"' && !in_single {
            in_double = !in_double;
            current.push(ch);
        } else if !in_single && !in_double && (ch.is_whitespace() || ch == ',') {
            if !current.is_empty() {
                tokens.push(std::mem::take(&mut current));
            }
        } else {
            current.push(ch);
        }
    }
    if !current.is_empty() {
        tokens.push(current);
    }
    tokens
}

fn parse_nexus_charset_rhs(
    rhs: &str,
    charsets: &IndexMap<String, Vec<SiteRange>>,
) -> Result<Vec<SiteRange>> {
    let mut ranges = Vec::new();
    let mut raw_buf = String::new();
    for token in split_quote_aware_tokens(rhs) {
        if let Some((_, ref_ranges)) = lookup_charset(charsets, &token) {
            if !raw_buf.trim().is_empty() {
                ranges.extend(parse_site_ranges(&raw_buf)?);
                raw_buf.clear();
            }
            ranges.extend(ref_ranges.iter().copied());
        } else {
            raw_buf.push_str(&token);
            raw_buf.push(' ');
        }
    }
    if !raw_buf.trim().is_empty() {
        ranges.extend(parse_site_ranges(&raw_buf)?);
    }
    if ranges.is_empty() {
        bail!("Empty charset definition '{rhs}'");
    }
    Ok(ranges)
}

fn split_charpartition_entries(rhs: &str) -> Vec<String> {
    let mut entries = Vec::new();
    let mut current = String::new();
    let mut in_single = false;
    let mut in_double = false;
    let mut brace_depth: usize = 0;

    for ch in rhs.chars() {
        if ch == '\'' && !in_double {
            in_single = !in_single;
            current.push(ch);
        } else if ch == '"' && !in_single {
            in_double = !in_double;
            current.push(ch);
        } else if ch == '{' && !in_single && !in_double {
            brace_depth += 1;
            current.push(ch);
        } else if ch == '}' && !in_single && !in_double && brace_depth > 0 {
            brace_depth -= 1;
            current.push(ch);
        } else if ch == ',' && !in_single && !in_double && brace_depth == 0 {
            let trimmed = current.trim();
            if !trimmed.is_empty() {
                entries.push(trimmed.to_string());
            }
            current.clear();
        } else {
            current.push(ch);
        }
    }
    let trimmed = current.trim();
    if !trimmed.is_empty() {
        entries.push(trimmed.to_string());
    }
    entries
}

fn resolve_charpartition(
    rhs: &str,
    charsets: &IndexMap<String, Vec<SiteRange>>,
) -> Result<Vec<Partition>> {
    let mut partitions = Vec::new();
    let mut used_names = HashSet::new();
    let mut subset_counter = 1usize;

    for entry in split_charpartition_entries(rhs) {
        let (raw_label, spec) = entry.split_once(':').with_context(|| {
            format!("Invalid charpartition entry '{entry}': missing ':' separator")
        })?;
        let clean_label = strip_outer_quotes(raw_label.trim()).trim().to_string();

        let mut charset_refs: Vec<String> = Vec::new();
        let mut ranges: Vec<SiteRange> = Vec::new();
        let mut raw_buf = String::new();
        let mut has_raw_ranges = false;

        for token in split_quote_aware_tokens(spec) {
            if let Some((cname, ref_ranges)) = lookup_charset(charsets, &token) {
                if !raw_buf.trim().is_empty() {
                    ranges.extend(parse_site_ranges(&raw_buf)?);
                    raw_buf.clear();
                    has_raw_ranges = true;
                }
                charset_refs.push(cname.clone());
                ranges.extend(ref_ranges.iter().copied());
            } else {
                raw_buf.push_str(&token);
                raw_buf.push(' ');
            }
        }
        if !raw_buf.trim().is_empty() {
            ranges.extend(parse_site_ranges(&raw_buf)?);
            has_raw_ranges = true;
        }
        if ranges.is_empty() {
            bail!("Empty site specification in charpartition entry '{entry}'");
        }

        let (gene, model) = if charset_refs.is_empty() {
            (sanitize_gene_name(&clean_label)?, None)
        } else if charset_refs.len() == 1 && !has_raw_ranges {
            let cname = charset_refs[0].clone();
            let model = if clean_label.is_empty() || clean_label.eq_ignore_ascii_case(&cname) {
                None
            } else {
                Some(clean_label)
            };
            (cname, model)
        } else {
            let model = if clean_label.is_empty() {
                None
            } else {
                Some(clean_label)
            };
            let gene = pick_merged_subset_name(
                &charset_refs,
                has_raw_ranges,
                &used_names,
                &mut subset_counter,
            );
            (gene, model)
        };

        used_names.insert(gene.clone());
        partitions.push(Partition::from_ranges(gene, ranges, model));
    }

    Ok(partitions)
}

fn pick_merged_subset_name(
    charset_refs: &[String],
    has_raw_ranges: bool,
    used_names: &HashSet<String>,
    subset_counter: &mut usize,
) -> String {
    if !has_raw_ranges && charset_refs.len() <= 3 {
        let joined = charset_refs.join("_");
        if joined.len() <= 64 && !used_names.contains(&joined) {
            return joined;
        }
    }
    loop {
        let candidate = format!("Subset{}", *subset_counter);
        *subset_counter += 1;
        if !used_names.contains(&candidate) {
            return candidate;
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::io::Cursor;

    macro_rules! test_partition_parser {
        ($input:expr_2021, $format:ident) => {
            let path = Path::new($input);
            let parser = PartitionParser::new(path, &PartitionFmt::$format, true);
            let partitions = parser.parse();
            assert_eq!(partitions.len(), 3);
            assert_eq!(partitions[0].gene, "locus1");
            assert_eq!(partitions[0].start, 1);
            assert_eq!(partitions[0].end, 5);
            assert_eq!(partitions[1].gene, "locus2");
            assert_eq!(partitions[1].start, 6);
            assert_eq!(partitions[1].end, 10);
            assert_eq!(partitions[2].gene, "locus3");
            assert_eq!(partitions[2].start, 11);
            assert_eq!(partitions[2].end, 15);
        };
    }

    #[test]
    fn test_parse_partition_raxml() {
        test_partition_parser!("tests/files/partition/partition.txt", Raxml);
    }

    #[test]
    fn test_parse_partition_nexus() {
        test_partition_parser!("tests/files/partition/partition.nex", Nexus);
    }

    #[test]
    fn test_parse_partition_raxml_with_whitespaces() {
        let path = Path::new("tests/files/partition/partition_whitespaces.txt");
        test_partition_parser!(path, Raxml);
    }

    #[test]
    fn test_partition_nexus_with_whitespaces() {
        let path = Path::new("tests/files/partition/partition_whitespaces.nex");
        test_partition_parser!(path, Nexus);
    }

    #[test]
    fn test_partition_raxml_no_datatype() {
        let path = Path::new("tests/files/partition/partition_no_datatype.txt");
        test_partition_parser!(path, Raxml);
    }

    #[test]
    #[should_panic]
    fn test_parse_partition_raxml_with_invalid_format() {
        let path = Path::new("tests/files/partition/partition_invalid.txt");
        let parser = PartitionParser::new(path, &PartitionFmt::Raxml, false);
        let _ = parser.parse();
    }

    #[test]
    fn test_parse_partition_raxml_codon() {
        let path = Path::new("tests/files/partition/partition_codon.txt");
        let parser = PartitionParser::new(path, &PartitionFmt::Raxml, false);
        let partitions = parser.try_parse().unwrap();
        assert_eq!(partitions.len(), 9);
        assert_eq!(partitions[0].gene, "locus1_Subset1");
        assert_eq!(partitions[0].ranges, vec![SiteRange::with_step(1, 5, 3)]);
        assert_eq!(partitions[0].model, None);
        assert_eq!(partitions[1].gene, "locus1_Subset2");
        assert_eq!(partitions[1].ranges, vec![SiteRange::with_step(2, 5, 3)]);
        assert_eq!(partitions[8].gene, "locus3_Subset3");
        assert_eq!(partitions[8].ranges, vec![SiteRange::with_step(13, 15, 3)]);
    }

    #[test]
    #[should_panic]
    fn test_parse_partition_raxml_inv() {
        let path = Path::new("tests/files/partition/partition_inv_pos.txt");
        let parser = PartitionParser::new(path, &PartitionFmt::Raxml, false);
        let _ = parser.parse();
    }

    #[test]
    #[should_panic]
    fn test_parse_partition_raxml_inv_start() {
        let path = Path::new("tests/files/partition/partition_inv_start.txt");
        let parser = PartitionParser::new(path, &PartitionFmt::Raxml, false);
        let _ = parser.parse();
    }

    #[test]
    fn test_parse_merged_raxml_with_models_and_comments() {
        let input = "\
# PartitionFinder scheme
GTR+I+G, Subset1 = 1-5, 11-15
// Comment line

HKY+G, Subset2 = 6-10
";
        let parser = PartitionParser::new(Path::new("dummy.txt"), &PartitionFmt::Raxml, false);
        let mut cursor = Cursor::new(input);
        let parts = parser.parse_raxml(&mut cursor).unwrap();
        assert_eq!(parts.len(), 2);
        assert_eq!(parts[0].gene, "Subset1");
        assert_eq!(parts[0].model.as_deref(), Some("GTR+I+G"));
        assert_eq!(
            parts[0].ranges,
            vec![SiteRange::new(1, 5), SiteRange::new(11, 15)]
        );
        assert_eq!(parts[0].start, 1);
        assert_eq!(parts[0].end, 15);
        assert_eq!(parts[1].gene, "Subset2");
        assert_eq!(parts[1].model.as_deref(), Some("HKY+G"));
    }

    #[test]
    fn test_parse_merged_strided_raxml() {
        let input = "\
GTR+G, Subset1 = 1-6\\3, 7-12\\3
GTR+G, Subset2 = 2-6\\3, 3-6\\3, 8-12\\3, 9-12\\3
";
        let parser = PartitionParser::new(Path::new("dummy.txt"), &PartitionFmt::Raxml, false);
        let mut cursor = Cursor::new(input);
        let parts = parser.parse_raxml(&mut cursor).unwrap();
        assert_eq!(parts.len(), 2);
        assert_eq!(
            parts[0].ranges,
            vec![
                SiteRange::with_step(1, 6, 3),
                SiteRange::with_step(7, 12, 3)
            ]
        );
        assert_eq!(parts[1].ranges.len(), 4);
    }

    #[test]
    fn test_parse_merged_nexus_with_comments_and_multiline() {
        let input = "\
#nexus
[ Partition file generated by [nested] tool ]
begin sets;
    charset 'Subset_1' = 1-4
        9-12;
    charset Subset_2 = 5-8, 13;
end;
";
        let parser = PartitionParser::new(Path::new("dummy.nex"), &PartitionFmt::Nexus, false);
        let mut cursor = Cursor::new(input);
        let parts = parser.parse_nexus(&mut cursor).unwrap();
        assert_eq!(parts.len(), 2);
        assert_eq!(parts[0].gene, "Subset_1");
        assert_eq!(
            parts[0].ranges,
            vec![SiteRange::new(1, 4), SiteRange::new(9, 12)]
        );
        assert_eq!(parts[1].gene, "Subset_2");
        assert_eq!(
            parts[1].ranges,
            vec![SiteRange::new(5, 8), SiteRange::single(13)]
        );
    }

    #[test]
    fn test_parse_iqtree_charpartition_merges_and_models() {
        let input = "\
#nexus
begin sets;
  charset locus1 = 1-5;
  charset locus2 = 6-10;
  charset locus3 = 11-15;
  charset locus4 = 16-20;
  charset locus5 = 21-25;
  charset locus6 = 26-30;
  charpartition mymodels =
    GTR+F+I+G4: locus1 locus3,
    HKY+F+G4: locus2,
    GTR+F+I+G4: locus4 locus5 locus6;
end;
";
        let parser = PartitionParser::new(Path::new("best.nex"), &PartitionFmt::Nexus, false);
        let mut cursor = Cursor::new(input);
        let parts = parser.parse_nexus(&mut cursor).unwrap();
        assert_eq!(parts.len(), 3);
        assert_eq!(parts[0].gene, "locus1_locus3");
        assert_eq!(parts[0].model.as_deref(), Some("GTR+F+I+G4"));
        assert_eq!(
            parts[0].ranges,
            vec![SiteRange::new(1, 5), SiteRange::new(11, 15)]
        );
        assert_eq!(parts[1].gene, "locus2");
        assert_eq!(parts[1].model.as_deref(), Some("HKY+F+G4"));
        assert_eq!(parts[2].gene, "locus4_locus5_locus6");
        assert_eq!(parts[2].model.as_deref(), Some("GTR+F+I+G4"));
    }

    #[test]
    fn test_parse_iqtree_large_merge_uses_subset_name() {
        let input = "\
#nexus
begin sets;
  charset g1 = 1-2;
  charset g2 = 3-4;
  charset g3 = 5-6;
  charset g4 = 7-8;
  charset g5 = 9-10;
  charpartition p = GTR+G: g1 g2 g3 g4, HKY: g5;
end;
";
        let parser = PartitionParser::new(Path::new("best.nex"), &PartitionFmt::Nexus, false);
        let mut cursor = Cursor::new(input);
        let parts = parser.parse_nexus(&mut cursor).unwrap();
        assert_eq!(parts.len(), 2);
        assert_eq!(parts[0].gene, "Subset1");
        assert_eq!(parts[0].ranges.len(), 4);
        assert_eq!(parts[1].gene, "g5");
    }

    #[test]
    fn test_parse_partitionfinder_skips_duplicate_mrbayes_block() {
        let input = "\
#nexus
begin sets;
    charset Subset1 = 1-5 11-15;
    charset Subset2 = 6-10;
    charpartition PartitionFinder = GTR+I+G:Subset1, HKY+G:Subset2;
end;
begin mrbayes;
    charset Subset1 = 1-5 11-15;
    charset Subset2 = 6-10;
end;
";
        let parser = PartitionParser::new(Path::new("pf.nex"), &PartitionFmt::Nexus, false);
        let mut cursor = Cursor::new(input);
        let parts = parser.parse_nexus(&mut cursor).unwrap();
        assert_eq!(parts.len(), 2);
        assert_eq!(parts[0].gene, "Subset1");
        assert_eq!(parts[0].model.as_deref(), Some("GTR+I+G"));
    }

    #[test]
    fn test_validate_rejects_overlap_gap_and_duplicates() {
        let parser = PartitionParser::new(Path::new("dummy.txt"), &PartitionFmt::Raxml, false);

        let overlap = "p1 = 1-5, 10-12\np2 = 4-9\n";
        assert!(parser.parse_raxml(&mut Cursor::new(overlap)).is_err());

        let gap = "p1 = 1-4\np2 = 6-10\n";
        assert!(parser.parse_raxml(&mut Cursor::new(gap)).is_err());

        let strided_gap = "p1 = 1-6\\3\np2 = 2-6\\3\n";
        assert!(parser.parse_raxml(&mut Cursor::new(strided_gap)).is_err());

        let dup_name = "p1 = 1-5\np1 = 6-10\n";
        assert!(parser.parse_raxml(&mut Cursor::new(dup_name)).is_err());
    }
}
