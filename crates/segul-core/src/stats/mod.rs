//! Sequence summary statistics
//!
//! ## Supported statistics:
//!
//! 1. Read sequences in FASTQ format, including gzipped FASTQ
//! 2. Contig sequences in FASTA format, including gzipped FASTA
//! 3. Sequence Alignment in FASTA, NEXUS, and relaxed-PHYLIP formats
//! 4. Pairwise genetic distances of sequence alignments
pub mod common;
#[cfg(feature = "contig")]
pub mod contigs;
#[cfg(feature = "align")]
pub mod distance;
#[cfg(feature = "read")]
pub mod fastq;
#[cfg(feature = "read")]
pub mod qscores;
#[cfg(feature = "read")]
pub mod read;
#[cfg(feature = "align")]
pub mod sequence;
