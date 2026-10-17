//! # Utilities for alignment manipulation and summary statistics
//! SEGUL provides fast and memory efficient utilities for alignment manipulation and summary statistics.
//!
//! ## Features:
//! 1. Alignment concatenation
//! 2. Alignment conversion
//! 3. Alignment filtering
//! 4. Alignment pairwise genetic distances
//! 5. Alignment partition conversion
//! 6. Alignment splitting
//! 7. Alignment summary statistics
//! 8. Unalign alignments
//! 9. Contiguous sequence summary statistics
//! 10. Raw read summary statistics
//! 11. Sample distribution mapping
//! 12. Sequence extraction
//! 13. Sequence filtering
//! 14. Sequence ID renaming
//! 15. Sequence removal
//! 16. Sequence translation
//! 17. Sequence transposition
//! 18. Sequence unique ID parsing
//!
//! ## Cargo features
//!
//! All features are enabled by default. Disable default features and pick
//! the ones you need to reduce compile time and dependencies:
//!
//! | Feature    | Provides                                                  |
//! |------------|-----------------------------------------------------------|
//! | `align`    | Alignment concat, convert, distance, filter, split, summary, trim |
//! | `sequence` | Sequence extract, filter, rename, remove, translate, transpose, ID |
//! | `read`     | Raw read (FASTQ) summary statistics                       |
//! | `contig`   | Contiguous sequence summary statistics                    |
//! | `maf`      | Multiple alignment format (MAF) conversion                |
//!
//! ```toml
//! segul-core = { version = "0.24", default-features = false, features = ["read"] }
//! ```
//!
//! ## Example
//!
//! ### Convert fasta to philip alignment
//! ```rust
//! use std::path::Path;
//!
//! use tempdir::TempDir;
//!
//! use segul_core::helper::types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix};
//! use segul_core::helper::sequence::SeqParser;
//! use segul_core::writer::sequences::SeqWriter;
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
#![cfg_attr(docsrs, feature(doc_cfg))]

pub mod core;
pub mod helper;
pub mod parser;
pub mod stats;
pub mod writer;
