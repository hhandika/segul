//! All the parsers for the supported formats.
#[cfg(feature = "maf")]
pub mod bed;
pub mod delimited;
pub mod fasta;
#[cfg(feature = "maf")]
pub mod maf;
pub mod nexus;
pub mod partition;
pub mod phylip;
#[cfg(feature = "read")]
pub mod qscores;
#[cfg(feature = "plot")]
pub mod summary_csv;
pub mod txt;

#[cfg(feature = "maf")]
const END_OF_LINE: u8 = b'\n';
#[cfg(feature = "maf")]
const EOF: usize = 0;

/// Default buffered I/O size (64 KB) for sequence parsers.
pub(crate) const IO_BUF_SIZE: usize = 64 * 1024;

/// Upper bound for pre-allocating interleaved sequences,
/// guarding against invalid NCHAR values in the header.
pub(crate) const MAX_RESERVE: usize = 1 << 30;

/// Split a matrix line into its ID and sequence bytes.
/// Returns None unless the line has exactly two whitespace-separated tokens.
pub(crate) fn split_id_seq(line: &[u8]) -> Option<(&[u8], &[u8])> {
    let line = line.trim_ascii();
    let split = line.iter().position(u8::is_ascii_whitespace)?;
    let (id, seq) = (&line[..split], line[split..].trim_ascii_start());
    if seq.is_empty() || has_whitespace(seq) {
        return None;
    }
    Some((id, seq))
}

fn has_whitespace(bytes: &[u8]) -> bool {
    // ASCII whitespace and control characters are all <= b' '.
    // Branch-free per chunk so the compiler can vectorize the scan.
    // Matches are then confirmed with the exact check.
    bytes
        .chunks(64)
        .any(|c| c.iter().fold(false, |acc, &b| acc | (b <= b' ')))
        && bytes.iter().any(|b| b.is_ascii_whitespace())
}

/// Pre-allocated builder for interleaved NEXUS and PHYLIP alignments.
/// Writes sequence chunks directly into a single flat `ntax * nchar` buffer
/// when header dimensions are valid, falling back to per-row buffers if a
/// sequence or taxon count exceeds the declared header.
pub(crate) struct InterleavedMatrixBuilder {
    ids: indexmap::IndexSet<String>,
    ntax: usize,
    nchar: usize,
    flat: Option<FlatInterleaved>,
    fallback: Vec<Vec<u8>>,
}

struct FlatInterleaved {
    data: Vec<u8>,
    cursors: Vec<usize>,
}

impl InterleavedMatrixBuilder {
    pub(crate) fn new(ntax: usize, nchar: usize) -> Self {
        let flat = if ntax > 0 && nchar > 0 {
            ntax.checked_mul(nchar)
                .filter(|&tot| tot <= MAX_RESERVE)
                .map(|tot| FlatInterleaved {
                    data: vec![0u8; tot],
                    cursors: vec![0; ntax],
                })
        } else {
            None
        };
        Self {
            ids: indexmap::IndexSet::with_capacity(ntax),
            ntax,
            nchar,
            flat,
            fallback: Vec::new(),
        }
    }

    #[inline]
    pub(crate) fn id_at(&self, idx: usize) -> Option<&str> {
        self.ids.get_index(idx).map(String::as_str)
    }

    #[inline]
    pub(crate) fn get_by_id(&self, id: &str) -> Option<&[u8]> {
        let idx = self.ids.get_index_of(id)?;
        Some(self.row_slice(idx))
    }

    /// Append `chunk` to the row for `id`, inserting `id` if not seen yet.
    pub(crate) fn append_by_id(&mut self, id: &str, chunk: &[u8]) {
        if let Some(idx) = self.ids.get_index_of(id) {
            self.append_by_index(idx, chunk);
        } else {
            let idx = self.ids.len();
            self.ids.insert(id.to_string());
            if let Some(flat) = &mut self.flat
                && idx < self.ntax
                && chunk.len() <= self.nchar
            {
                let start = idx * self.nchar;
                flat.data[start..start + chunk.len()].copy_from_slice(chunk);
                flat.cursors[idx] = chunk.len();
                return;
            }
            if self.flat.is_some() {
                self.spill_to_fallback(idx);
            }
            self.fallback.push(chunk.to_vec());
        }
    }

    /// Append `chunk` to an already-inserted row `idx`.
    #[inline]
    pub(crate) fn append_by_index(&mut self, idx: usize, chunk: &[u8]) {
        if let Some(flat) = &mut self.flat {
            let cursor = flat.cursors[idx];
            if cursor + chunk.len() <= self.nchar {
                let start = idx * self.nchar + cursor;
                flat.data[start..start + chunk.len()].copy_from_slice(chunk);
                flat.cursors[idx] += chunk.len();
                return;
            }
            self.spill_to_fallback(self.ids.len());
        }
        self.fallback[idx].extend_from_slice(chunk);
    }

    pub(crate) fn into_matrix(self) -> crate::helper::types::SeqMatrix {
        if let Some(flat) = self.flat {
            let uniform =
                self.ids.len() == self.ntax && flat.cursors.iter().all(|&len| len == self.nchar);
            if uniform {
                return crate::helper::types::SeqMatrix::from_uniform(
                    self.ids, flat.data, self.nchar,
                );
            }
            let mut offsets = Vec::with_capacity(self.ids.len() + 1);
            offsets.push(0);
            let mut data = Vec::with_capacity(flat.cursors.iter().sum());
            for (idx, &len) in flat.cursors.iter().take(self.ids.len()).enumerate() {
                let start = idx * self.nchar;
                data.extend_from_slice(&flat.data[start..start + len]);
                offsets.push(data.len());
            }
            return crate::helper::types::SeqMatrix::from_parts(self.ids, data, offsets);
        }

        let total: usize = self.fallback.iter().map(Vec::len).sum();
        let mut data = Vec::with_capacity(total);
        let mut offsets = Vec::with_capacity(self.ids.len() + 1);
        offsets.push(0);
        for row in self.fallback {
            data.extend_from_slice(&row);
            offsets.push(data.len());
        }
        crate::helper::types::SeqMatrix::from_parts(self.ids, data, offsets)
    }

    fn row_slice(&self, idx: usize) -> &[u8] {
        match &self.flat {
            Some(flat) => {
                let start = idx * self.nchar;
                &flat.data[start..start + flat.cursors[idx]]
            }
            None => &self.fallback[idx],
        }
    }

    #[cold]
    fn spill_to_fallback(&mut self, existing_count: usize) {
        if let Some(flat) = self.flat.take() {
            self.fallback = Vec::with_capacity(existing_count.max(self.ntax));
            for (idx, &len) in flat.cursors.iter().take(existing_count).enumerate() {
                let start = idx * self.nchar;
                self.fallback.push(flat.data[start..start + len].to_vec());
            }
        }
    }
}

macro_rules! insert_matrix {
    ($self: ident, $id: ident, $seq: ident ) => {
        match $self.matrix.get($id) {
            Some(original_seq) => panic!(
                "Found duplicate IDs for file {}. \
                First duplicate ID found: {}. \
                Both sequences are the same: {}.",
                $self.input.display(),
                $id,
                original_seq == $seq
            ),
            None => {
                $self.matrix.insert($id.to_string(), $seq);
            }
        }
    };
}

macro_rules! warn_duplicate_ids {
    ($self: ident, $ids: ident) => {
        if $ids.len() != $self.header.ntax {
            log::warn!(
                "{} Found problematic alignment: {}.\n\
                Number of taxa in the matrix ({} taxa) \
                does not match the number of taxa in the header ({} taxa).\n\
                If this is `id` subcommand, try to use the same file as an input \
                for other subcommands.\n\
                If the cause is duplicate IDs, \
                the other subcommands will show the first found duplicate.\n\n",
                "WARNING!".red(),
                $self.input.display(),
                $ids.len(),
                $self.header.ntax
            );
        }
    };
}

pub(crate) use insert_matrix;
pub(crate) use warn_duplicate_ids;

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_split_id_seq() {
        assert_eq!(
            Some((b"ABCD".as_slice(), b"ACGT".as_slice())),
            split_id_seq(b"ABCD ACGT")
        );
        assert_eq!(
            Some((b"ABCD".as_slice(), b"ACGT".as_slice())),
            split_id_seq(b" ABCD \t  ACGT\r\n")
        );
        assert_eq!(None, split_id_seq(b"AB CD ACGT"));
        assert_eq!(None, split_id_seq(b"ABCD ACG T"));
        assert_eq!(None, split_id_seq(b"ABCD"));
        assert_eq!(None, split_id_seq(b""));
    }
}
