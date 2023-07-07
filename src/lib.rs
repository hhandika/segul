//! # Utilities for alignment manipulation and summary statistics
//! SEGUL provides fast and memory efficient utilities for alignment manipulation and summary statistics.
//!
//! ## Features:
//! 1. Alignment concatenation
//! 2. Alignment conversion
//! 3. Alignment filtering
//! 4. Alignment partition conversion
//! 5. Alignment splitting
//! 6. Alignment summary statistics
//! 7. Contiguous sequence summary statistics
//! 8. Raw read summary statistics
//! 9. Sample distribution mapping
//! 10. Sequence extraction
//! 11. Sequence ID renaming
//! 12. Sequence removal
//! 13. Sequence translation
//! 14. Sequence unique ID parsing
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
