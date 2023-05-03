//! # Utilities for alignment manipulation and summary statistics
//! SEGUL provides fast and memory efficient utilities for alignment manipulation and summary statistics.
//!
//! ## Features:
//! 1. Computing alignment summary statistics.
//! 2. Concatenating alignments with partition settings.
//! 3. Converting alignments to different formats.
//! 4. Converting partition formats.
//! 5. Extracting sequences from a collection of alignments based on user-defined IDs (include regular expression support).
//! 6. Filtering alignments based on minimal taxon completeness, alignment length, or numbers of parsimony informative sites.
//! 7. Getting sample IDs from a collection of alignments.
//! 8. Mapping sample distribution in a collection of alignments.
//! 9. Batch removing sequence based user-defined IDs.
//! 10. Batch renaming sequence IDs (include regular expression support).
//! 11. Splitting alignments by partitions.
//! 12. Translating DNA sequences to amino acid sequences
//!
//! ## Example
//!
//! ### Convert fasta to philip alignment
//! ```rust, ignore
//! use std::path::Path;
//! use segul::helper::types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix};
//! use segul::helper::sequence::SeqParser;
//! use segul::writer::sequences::SeqWriter;
//!
//! let input = Path::new("tests/files/simple.fas");
//! let output = Path::new("tests/data/alignment.phy");
//! let input_fmt = InputFmt::Fasta;
//! let output_fmt = OutputFmt::Phylip;
//! let datatype = DataType::Dna;
//!
//! let (sequence_matrix, header) = SeqParser::new(input, &datatype).parse(&input_fmt);
//! let writer = SeqWriter::new(output, &sequence_matrix, &header);
//! writer.write_sequence(&output_fmt).unwrap()
//! ```
pub mod cli;
pub mod handler;
pub mod helper;
pub mod parser;
pub mod stats;
pub mod writer;
