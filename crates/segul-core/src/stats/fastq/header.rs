//! Lightweight FASTQ header inspection for platform, filter, and run metadata.
//!
//! Headers are parsed from the borrowed record slice without allocating on the
//! hot path. Distinct run/lane identifiers are bounded so memory stays small
//! even on malformed inputs.
use std::collections::BTreeSet;

/// Maximum number of distinct run/flowcell/lane groups retained per file.
const MAX_RUN_GROUPS: usize = 64;
const MAX_GROUP_LEN: usize = 128;

const PLATFORM_ILLUMINA: usize = 0;
const PLATFORM_NANOPORE: usize = 1;
const PLATFORM_PACBIO: usize = 2;

/// Streaming, mergeable summary of FASTQ header metadata.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct HeaderAccumulator {
    platform_counts: [u64; 3],
    pub(crate) chastity_failed: u64,
    run_groups: BTreeSet<String>,
}

impl HeaderAccumulator {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Inspect one raw FASTQ header line (including the leading `@`).
    pub(crate) fn add(&mut self, header: &[u8]) {
        let line = header.strip_prefix(b"@").unwrap_or(header);
        if line.is_empty() {
            return;
        }
        let (id, comment) = split_first_whitespace(line);
        if let Some(comment) = comment
            && is_nanopore_comment(comment)
        {
            self.platform_counts[PLATFORM_NANOPORE] += 1;
            if let Some(run_id) = extract_kv(comment, b"runid=") {
                self.record_group(&[b"runid=", run_id]);
            }
            return;
        }
        if self.try_illumina(id, comment) {
            return;
        }
        self.try_pacbio(id);
    }

    /// Combine header counts from another block of the same file.
    pub(crate) fn merge(&mut self, other: Self) {
        for (a, b) in self.platform_counts.iter_mut().zip(other.platform_counts) {
            *a += b;
        }
        self.chastity_failed += other.chastity_failed;
        self.run_groups.extend(other.run_groups);
        while self.run_groups.len() > MAX_RUN_GROUPS {
            self.run_groups.pop_last();
        }
    }

    /// Dominant sequencing platform inferred from headers, if recognized.
    pub(crate) fn platform(&self) -> Option<&'static str> {
        let names = ["Illumina", "Nanopore", "PacBio"];
        self.platform_counts
            .iter()
            .enumerate()
            .filter(|(_, count)| **count > 0)
            .max_by_key(|(_, count)| *count)
            .map(|(i, _)| names[i])
    }

    pub(crate) fn run_groups(&self) -> Vec<String> {
        self.run_groups.iter().cloned().collect()
    }

    pub(crate) fn run_group_count(&self) -> u64 {
        self.run_groups.len() as u64
    }

    fn try_illumina(&mut self, id: &[u8], comment: Option<&[u8]>) -> bool {
        let mut parts = id.split(|&b| b == b':');
        let (Some(p0), Some(p1), Some(p2), Some(p3), Some(p4)) = (
            parts.next(),
            parts.next(),
            parts.next(),
            parts.next(),
            parts.next(),
        ) else {
            return false;
        };
        match (parts.next(), parts.next(), parts.next()) {
            // Casava 1.8+ / DRAGEN: instrument:run:flowcell:lane:tile:x:y
            (Some(p5), Some(p6), None)
                if !p0.is_empty()
                    && !p2.is_empty()
                    && is_digits(p3)
                    && is_digits(p4)
                    && is_digits(p5)
                    && is_digits(strip_coord_suffix(p6)) =>
            {
                self.platform_counts[PLATFORM_ILLUMINA] += 1;
                self.record_group(&[p0, b":", p2, b":", p3]);
                if let Some(comment) = comment
                    && is_chastity_failed(comment)
                {
                    self.chastity_failed += 1;
                }
                true
            }
            // Legacy Casava 1.4: instrument:lane:tile:x:y
            (None, None, None)
                if !p0.is_empty()
                    && is_digits(p1)
                    && is_digits(p2)
                    && is_digits(p3)
                    && is_digits(strip_coord_suffix(p4)) =>
            {
                self.platform_counts[PLATFORM_ILLUMINA] += 1;
                self.record_group(&[p0, b":", p1]);
                true
            }
            _ => false,
        }
    }

    fn try_pacbio(&mut self, id: &[u8]) -> bool {
        let mut parts = id.split(|&b| b == b'/');
        let (Some(movie), Some(zmw), Some(kind), None) =
            (parts.next(), parts.next(), parts.next(), parts.next())
        else {
            return false;
        };
        if movie.is_empty() || !is_digits(zmw) || !is_pacbio_subread_or_ccs(kind) {
            return false;
        }
        self.platform_counts[PLATFORM_PACBIO] += 1;
        self.record_group(&[movie]);
        true
    }

    fn record_group(&mut self, pieces: &[&[u8]]) {
        let mut buf = [0u8; MAX_GROUP_LEN];
        let mut len = 0;
        for piece in pieces {
            let take = piece.len().min(MAX_GROUP_LEN - len);
            buf[len..len + take].copy_from_slice(&piece[..take]);
            len += take;
        }
        let Ok(group) = std::str::from_utf8(&buf[..len]) else {
            return;
        };
        if self.run_groups.contains(group) {
            return;
        }
        if self.run_groups.len() < MAX_RUN_GROUPS {
            self.run_groups.insert(group.to_string());
        } else if let Some(last) = self.run_groups.last()
            && group < last.as_str()
        {
            self.run_groups.pop_last();
            self.run_groups.insert(group.to_string());
        }
    }
}

fn split_first_whitespace(line: &[u8]) -> (&[u8], Option<&[u8]>) {
    match line.iter().position(|&b| b == b' ' || b == b'\t') {
        Some(pos) => {
            let rest = &line[pos + 1..];
            let trimmed = rest
                .iter()
                .position(|&b| b != b' ' && b != b'\t')
                .map(|start| &rest[start..]);
            (&line[..pos], trimmed)
        }
        None => (line, None),
    }
}

fn is_digits(bytes: &[u8]) -> bool {
    !bytes.is_empty() && bytes.iter().all(u8::is_ascii_digit)
}

fn strip_coord_suffix(bytes: &[u8]) -> &[u8] {
    let end = bytes
        .iter()
        .position(|&b| b == b'/' || b == b'#')
        .unwrap_or(bytes.len());
    &bytes[..end]
}

fn is_chastity_failed(comment: &[u8]) -> bool {
    let (token, _) = split_first_whitespace(comment);
    let mut fields = token.split(|&b| b == b':');
    let (Some(read_num), Some(filtered), Some(_control)) =
        (fields.next(), fields.next(), fields.next())
    else {
        return false;
    };
    is_digits(read_num) && filtered == b"Y"
}

fn is_nanopore_comment(comment: &[u8]) -> bool {
    comment.split(|&b| b == b' ' || b == b'\t').any(|tok| {
        tok.starts_with(b"runid=") || tok.starts_with(b"ch=") || tok.starts_with(b"start_time=")
    })
}

fn extract_kv<'a>(comment: &'a [u8], prefix: &[u8]) -> Option<&'a [u8]> {
    comment
        .split(|&b| b == b' ' || b == b'\t')
        .find_map(|tok| tok.strip_prefix(prefix))
        .filter(|val| !val.is_empty())
}

fn is_pacbio_subread_or_ccs(kind: &[u8]) -> bool {
    if kind == b"ccs" {
        return true;
    }
    let mut bounds = kind.split(|&b| b == b'_');
    matches!(
        (bounds.next(), bounds.next(), bounds.next()),
        (Some(start), Some(end), None) if is_digits(start) && is_digits(end)
    )
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_illumina_headers_and_chastity() {
        let mut acc = HeaderAccumulator::new();
        acc.add(b"@A00123:45:HXYZ1DSXX:1:1101:1234:5678 1:N:0:ATCACG");
        acc.add(b"@A00123:45:HXYZ1DSXX:1:1101:1235:5679 1:Y:0:ATCACG");
        acc.add(b"@A00123:45:HXYZ1DSXX:2:1101:1236:5680 1:N:0:ATCACG");
        assert_eq!(acc.platform(), Some("Illumina"));
        assert_eq!(acc.chastity_failed, 1);
        assert_eq!(acc.run_group_count(), 2);
        assert_eq!(
            acc.run_groups(),
            vec!["A00123:HXYZ1DSXX:1", "A00123:HXYZ1DSXX:2"]
        );
    }

    #[test]
    fn test_nanopore_headers() {
        let mut acc = HeaderAccumulator::new();
        acc.add(b"@f7a8b9c0-1234-5678-9abc-def012345678 runid=abc123 read=12 ch=42 start_time=2025-01-01T00:00:00Z");
        assert_eq!(acc.platform(), Some("Nanopore"));
        assert_eq!(acc.chastity_failed, 0);
        assert_eq!(acc.run_groups(), vec!["runid=abc123"]);
    }

    #[test]
    fn test_pacbio_headers() {
        let mut acc = HeaderAccumulator::new();
        acc.add(b"@m64011_190830_220126/12345/ccs");
        acc.add(b"@m64011_190830_220126/12346/0_15000");
        assert_eq!(acc.platform(), Some("PacBio"));
        assert_eq!(acc.run_group_count(), 1);
        assert_eq!(acc.run_groups(), vec!["m64011_190830_220126"]);
    }

    #[test]
    fn test_generic_and_merge() {
        let mut first = HeaderAccumulator::new();
        first.add(b"@SRR12345.1 length=150");
        assert_eq!(first.platform(), None);
        assert_eq!(first.run_group_count(), 0);

        let mut second = HeaderAccumulator::new();
        second.add(b"@M001:1:FC1:1:1101:10:20 1:Y:0:AAAA");
        first.merge(second);
        assert_eq!(first.platform(), Some("Illumina"));
        assert_eq!(first.chastity_failed, 1);
        assert_eq!(first.run_group_count(), 1);
    }
}
