//! Normalize MAF blocks: sample naming and reference orientation.
use std::collections::HashSet;
use std::path::Path;

use anyhow::{Context, Result, bail};
use indexmap::IndexMap;

use crate::{
    helper::{alphabet, files, types::DnaStrand},
    parser::maf::{MafAlignment, MafParagraph, MafReader},
};

/// How MAF source names (e.g. `hg18.chr7`) become sample names.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum SampleNaming {
    /// Use the text before the first dot (`hg18.chr7` becomes `hg18`).
    /// Sources without a dot are used as they are.
    #[default]
    AssemblyPrefix,
    /// Use the full source name. Each chromosome becomes its own sample.
    FullSource,
}

/// Map MAF source names to sample names and chromosomes.
pub(crate) struct SpeciesResolver {
    naming: SampleNaming,
    // Listed species, longest first, so the most specific name matches.
    species: Vec<String>,
    keep: Option<HashSet<String>>,
}

impl SpeciesResolver {
    /// Listed species name sources that contain dots, such as
    /// `GCF_000001405.40.NC_000001.11`, and limit output to those species.
    pub(crate) fn new(naming: SampleNaming, species: Option<&[String]>) -> Self {
        let mut listed: Vec<String> = species.map(|s| s.to_vec()).unwrap_or_default();
        listed.sort_by(|a, b| b.len().cmp(&a.len()).then(a.cmp(b)));
        Self {
            naming,
            keep: species.map(|s| s.iter().cloned().collect()),
            species: listed,
        }
    }

    /// Sample name of a MAF source.
    pub(crate) fn sample(&self, source: &str) -> String {
        if let Some(name) = self.listed_prefix(source) {
            return name.to_string();
        }
        match self.naming {
            SampleNaming::FullSource => source.to_string(),
            SampleNaming::AssemblyPrefix => source.split('.').next().unwrap_or(source).to_string(),
        }
    }

    /// Chromosome part of a MAF source: the text after the sample prefix and a dot,
    /// or the whole source when there is no dot.
    pub(crate) fn chrom<'a>(&self, source: &'a str) -> &'a str {
        let prefix_len = match self.listed_prefix(source) {
            Some(name) => name.len(),
            None => source.find('.').unwrap_or(source.len()),
        };
        match source.get(prefix_len + 1..) {
            Some(chrom) if !chrom.is_empty() => chrom,
            _ => source,
        }
    }

    /// Whether a sample is written to the output.
    pub(crate) fn is_kept(&self, sample: &str) -> bool {
        self.keep.as_ref().is_none_or(|keep| keep.contains(sample))
    }

    fn listed_prefix(&self, source: &str) -> Option<&str> {
        self.species
            .iter()
            .find(|name| {
                source == name.as_str()
                    || source
                        .strip_prefix(name.as_str())
                        .is_some_and(|rest| rest.starts_with('.'))
            })
            .map(|name| name.as_str())
    }
}

/// A MAF block oriented to the forward strand of its reference sequence.
pub(crate) struct RefBlock {
    /// Full MAF source of the reference row, e.g. `hg18.chr7`.
    pub source: String,
    /// Reference chromosome, e.g. `chr7`.
    pub chrom: String,
    /// Zero-based, half-open reference interval on the forward strand.
    pub start: usize,
    pub end: usize,
    pub src_size: usize,
    pub score: Option<f64>,
    /// The reference row, used for coordinates even if not written.
    pub ref_text: Vec<u8>,
    /// Rows written to the output, reference first when kept.
    pub rows: IndexMap<String, Vec<u8>>,
    /// Rows dropped because a sample appeared more than once.
    pub duplicates: usize,
}

impl RefBlock {
    /// Normalize an alignment block. Returns `None` when the block has no reference row.
    pub(crate) fn from_alignment(
        alignment: MafAlignment,
        resolver: &SpeciesResolver,
        ref_sample: &str,
    ) -> Result<Option<Self>> {
        let Some(ref_idx) = alignment
            .sequences
            .iter()
            .position(|seq| resolver.sample(&seq.source) == ref_sample)
        else {
            return Ok(None);
        };
        let ncols = alignment.sequences[ref_idx].text.len();
        if let Some(seq) = alignment
            .sequences
            .iter()
            .find(|seq| seq.text.len() != ncols)
        {
            bail!(
                "rows in a block have different lengths: {} has {} columns, expected {}",
                seq.source,
                seq.text.len(),
                ncols
            );
        }
        let reference = &alignment.sequences[ref_idx];
        let ref_bases = reference.text.iter().filter(|&&b| b != b'-').count();
        if ref_bases != reference.size {
            bail!(
                "reference {} at {} has {} bases but its size is {}",
                reference.source,
                reference.start,
                ref_bases,
                reference.size
            );
        }
        let is_reverse = reference.strand == DnaStrand::Reverse;
        let start = if is_reverse {
            reference
                .src_size
                .checked_sub(reference.start + reference.size)
                .with_context(|| format!("invalid coordinates for {}", reference.source))?
        } else {
            reference.start
        };
        let orient = |text: &[u8]| {
            if is_reverse {
                alphabet::reverse_complement_dna(text)
            } else {
                text.to_vec()
            }
        };
        let mut rows = IndexMap::new();
        let mut duplicates = 0;
        let order = std::iter::once(ref_idx)
            .chain((0..alignment.sequences.len()).filter(|&i| i != ref_idx));
        for idx in order {
            let seq = &alignment.sequences[idx];
            let sample = resolver.sample(&seq.source);
            if rows.contains_key(&sample) || (idx != ref_idx && sample == ref_sample) {
                duplicates += 1;
                continue;
            }
            if resolver.is_kept(&sample) {
                rows.insert(sample, orient(&seq.text));
            }
        }
        Ok(Some(Self {
            chrom: resolver.chrom(&reference.source).to_string(),
            source: reference.source.clone(),
            start,
            end: start + reference.size,
            src_size: reference.src_size,
            score: alignment.score,
            ref_text: orient(&reference.text),
            rows,
            duplicates,
        }))
    }

    /// Number of alignment columns.
    pub(crate) fn ncols(&self) -> usize {
        self.ref_text.len()
    }
}

/// Sample name of the first sequence in the first alignment block of the files.
pub(crate) fn first_sample(
    files: &[std::path::PathBuf],
    resolver: &SpeciesResolver,
) -> Result<String> {
    for file in files {
        if let Some(source) = first_source(file)? {
            return Ok(resolver.sample(&source));
        }
    }
    bail!("No MAF alignment blocks found in the input files")
}

fn first_source(file: &Path) -> Result<Option<String>> {
    let reader = MafReader::new(files::open_reader(file)?);
    for paragraph in reader {
        let paragraph = paragraph.with_context(|| format!("Failed reading {}", file.display()))?;
        if let MafParagraph::Alignment(alignment) = paragraph {
            return Ok(alignment.sequences.first().map(|seq| seq.source.clone()));
        }
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::maf::MafSequence;

    pub(crate) fn sequence(
        source: &str,
        start: usize,
        strand: DnaStrand,
        text: &str,
    ) -> MafSequence {
        MafSequence {
            source: source.to_string(),
            start,
            size: text.bytes().filter(|&b| b != b'-').count(),
            strand,
            src_size: 100,
            text: text.as_bytes().to_vec(),
        }
    }

    fn alignment(sequences: Vec<MafSequence>) -> MafAlignment {
        let mut aln = MafAlignment::new();
        aln.sequences = sequences;
        aln
    }

    #[test]
    fn test_default_naming() {
        let resolver = SpeciesResolver::new(SampleNaming::AssemblyPrefix, None);
        assert_eq!(resolver.sample("hg18.chr7"), "hg18");
        assert_eq!(resolver.chrom("hg18.chr7"), "chr7");
        assert_eq!(resolver.sample("baboon"), "baboon");
        assert_eq!(resolver.chrom("baboon"), "baboon");
        assert_eq!(resolver.chrom("hg18.chrUn.random"), "chrUn.random");
    }

    #[test]
    fn test_full_source_naming() {
        let resolver = SpeciesResolver::new(SampleNaming::FullSource, None);
        assert_eq!(resolver.sample("hg18.chr7"), "hg18.chr7");
        assert_eq!(resolver.chrom("hg18.chr7"), "chr7");
    }

    #[test]
    fn test_listed_species_with_dots() {
        let species = vec!["GCF_000001405".to_string(), "GCF_000001405.40".to_string()];
        let resolver = SpeciesResolver::new(SampleNaming::AssemblyPrefix, Some(&species));
        let source = "GCF_000001405.40.NC_000001.11";
        assert_eq!(resolver.sample(source), "GCF_000001405.40");
        assert_eq!(resolver.chrom(source), "NC_000001.11");
        assert!(resolver.is_kept("GCF_000001405.40"));
        assert!(!resolver.is_kept("mm10"));
    }

    #[test]
    fn test_reference_first_and_duplicates() {
        let resolver = SpeciesResolver::new(SampleNaming::AssemblyPrefix, None);
        let aln = alignment(vec![
            sequence("mm4.chr6", 0, DnaStrand::Forward, "AC-T"),
            sequence("hg18.chr7", 10, DnaStrand::Forward, "ACGT"),
            sequence("mm4.chr1", 0, DnaStrand::Forward, "ACGA"),
        ]);
        let block = RefBlock::from_alignment(aln, &resolver, "hg18")
            .unwrap()
            .unwrap();
        assert_eq!(block.rows.keys().collect::<Vec<_>>(), ["hg18", "mm4"]);
        assert_eq!(block.rows["mm4"], b"AC-T".to_vec());
        assert_eq!(block.duplicates, 1);
        assert_eq!(
            (block.chrom.as_str(), block.start, block.end),
            ("chr7", 10, 14)
        );
    }

    #[test]
    fn test_minus_strand_reference_is_normalized() {
        let resolver = SpeciesResolver::new(SampleNaming::AssemblyPrefix, None);
        let aln = alignment(vec![
            sequence("ref.chr1", 90, DnaStrand::Reverse, "ACGTT"),
            sequence("other.chr1", 0, DnaStrand::Forward, "AC-TT"),
        ]);
        let block = RefBlock::from_alignment(aln, &resolver, "ref")
            .unwrap()
            .unwrap();
        assert_eq!((block.start, block.end), (5, 10));
        assert_eq!(block.ref_text, b"AACGT".to_vec());
        assert_eq!(block.rows["other"], b"AA-GT".to_vec());
    }

    #[test]
    fn test_missing_reference_and_unequal_rows() {
        let resolver = SpeciesResolver::new(SampleNaming::AssemblyPrefix, None);
        let aln = alignment(vec![sequence("mm4.chr6", 0, DnaStrand::Forward, "ACGT")]);
        assert!(
            RefBlock::from_alignment(aln, &resolver, "hg18")
                .unwrap()
                .is_none()
        );
        let aln = alignment(vec![
            sequence("hg18.chr7", 0, DnaStrand::Forward, "ACGT"),
            sequence("mm4.chr6", 0, DnaStrand::Forward, "ACG"),
        ]);
        assert!(RefBlock::from_alignment(aln, &resolver, "hg18").is_err());
    }

    #[test]
    fn test_first_sample() {
        let resolver = SpeciesResolver::new(SampleNaming::AssemblyPrefix, None);
        let files = vec![std::path::PathBuf::from("tests/files/maf/simple.maf")];
        assert_eq!(first_sample(&files, &resolver).unwrap(), "hg18");
    }
}
