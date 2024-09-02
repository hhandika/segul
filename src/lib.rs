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
//! 7. Unalign alignments
//! 7. Contiguous sequence summary statistics
//! 8. Raw read summary statistics
//! 9. Sample distribution mapping
//! 10. Sequence extraction
//! 11. Sequence filtering
//! 12. Sequence ID renaming
//! 13. Sequence removal
//! 14. Sequence translation
//! 15. Sequence unique ID parsing
//!
//! ## Example
//!
//! ### Convert fasta to philip alignment
//! ```rust
//! use std::path::Path;
//!
//! use tempdir::TempDir;
//!
//! use segul::helper::types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix};
//! use segul::helper::sequence::SeqParser;
//! use segul::writer::sequences::SeqWriter;
//!
//! let input = Path::new("tests/files/simple.fas");
//! let input_fmt = InputFmt::Fasta;
//! let datatype = DataType::Dna;
//! // Replace binding with your output directory, for example:
//! // let output_path = Path::new("output").join("my_output_alignment");
//! let binding = TempDir::new("temp").unwrap();
//! let output_path = binding.path().join("my_output_alignment");
//! let output_fmt = OutputFmt::Phylip;
//! let (sequence_matrix, header) = SeqParser::new(input, &datatype).parse(&input_fmt);
//! let mut writer = SeqWriter::new(&output_path, &sequence_matrix, &header);
//! writer.write_sequence(&output_fmt).unwrap()
//! ```
pub mod cli;
pub mod core;
pub mod helper;
pub mod parser;
pub mod stats;
pub mod writer;
