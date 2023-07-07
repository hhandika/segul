//! Sequence summary statistics
//!
//! ## Supported statistics:
//!
//! 1. Read sequences in FASTQ format, including gzipped FASTQ
//! 2. Contig sequences in FASTA format, including gzipped FASTA
//! 3. Sequence Alignment in FASTA, NEXUS, and relaxed-PHYLIP formats
pub mod common;
pub mod contigs;
pub mod fastq;
pub mod qscores;
pub mod read;
pub mod sequence;
