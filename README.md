# SEGUL

![Segul-Tests](https://github.com/hhandika/segul/workflows/Segul-Tests/badge.svg)
![Crate-IO](https://img.shields.io/crates/v/segul)
![Crates-Download](https://img.shields.io/crates/d/segul?color=orange&label=crates.io-downloads)
![GH-Release](https://img.shields.io/github/v/tag/hhandika/segul?label=gh-releases)
![GH-Downloads](https://img.shields.io/github/downloads/hhandika/segul/total?color=blue&label=gh-release-downloads)
[![LoC](https://tokei.rs/b1/github/hhandika/segul?category=code)](https://github.com/XAMPPRocky/tokei)
![last-commit](https://img.shields.io/github/last-commit/hhandika/segul)
![License](https://img.shields.io/github/license/hhandika/segul)

SEGUL is an ultra-fast, memory-efficient application for working with sequence alignments. It is a cross-platform, single executable app, with zero runtime dependency. In simple words, you need not worry about dependencies. It is designed to handle operations on large genomic datasets, while using minimal computational resources. However, it also provides convenient features for working on smaller datasets (e.g., Sanger datasets). In our tests, it consistently offers a faster and more efficient (low memory footprint) alternative to existing applications for a variety of sequence alignment manipulations ([see benchmark](https://github.com/hhandika/segul-bench)).

The app is designed with data security and repeatability in mind. By default, it prevents overwriting files and checks whether input sequences contain only valid IUPAC characters. When concatenating alignment, the app checks if all the sequences in each input file are aligned. It will abort processing if the input sequences are not aligned. It also carries the memory and multi-threading safety features provided by [the Rust programming language](https://www.rust-lang.org/). The program also logs its terminal output for record keeping.

Features:

1. Computing alignment summary statistics.
2. Concatenating alignments with partition settings.
3. Converting alignments to different formats.
4. Converting partition formats.
5. Extracting sequences from a collection of alignments based on user-defined IDs (include regular expression support).
6. Filtering alignments based on minimal taxon completeness, alignment length, or numbers of parsimony informative sites.
7. Getting sample IDs from a collection of alignments.
8. Mapping sample distribution in a collection of alignments.
9. Batch removing sequence based user-defined IDs.
10. Batch renaming sequence IDs (include regular expression support).
11. Splitting alignments by partitions.
12. Translating DNA sequences to amino acid sequences

Upcoming features and bug fixes ([try](https://github.com/hhandika/segul/wiki/2.-Installation#try-development-features)):

1. GUI version is coming soon (see [SEGUL-GUI](https://github.com/hhandika/segui)).
2. New command line structure.
3. New feature to generate summary statistics for raw-read genomic sequences.
4. Bug fixes and update on deprecated dependencies.

Supported sequence formats:

1. NEXUS
2. Relaxed PHYLIP
3. FASTA

All of the formats are supported in interleave and sequential versions. The app supports both DNA and amino acid sequences.

Supported partition formats:

1. RaXML
2. NEXUS

The NEXUS partition can be written as a charset block embedded in NEXUS formatted sequences or a separate file.

Documentation: [Link](https://docs.page/hhandika/segul-docs/)

Citation: [pre-print](https://doi.org/10.22541/au.165167823.30911834/v1)

## Supported Platforms

The app may work in any Rust supported [platform](https://doc.rust-lang.org/nightly/rustc/platform-support.html). Below is a list of operating system that we tested and is guaranteed to work:

- Linux
- MacOS
- Windows
- Windows Subsystem for Linux (WSL)

> :warning: **SEGUL modern terminal output comes with a cost of requiring a terminal application that supports [UTF-8 encoding](https://en.wikipedia.org/wiki/UTF-8)**. For MacOS and native Linux, your default terminal should have supported UTF-8 encoding by default. For Windows (including WSL) users, we recommend using [Windows Terminal](https://www.microsoft.com/en-us/p/windows-terminal/9n0dx20hk701#activetab=pivot:overviewtab) to ensure consistent terminal output. Windows Terminal requires separate installation for Windows 10. It should come pre-installed on Windows 11.

## Contribution

We welcome any kind of contribution, from issue reporting, ideas to improve the app, to code contribution. For ideas and issue reporting please post in [the Github issues page](https://github.com/hhandika/segul/issues). For code contribution, please fork the repository and send pull requests to this repo

<!-- ## Acknowledgment

We thank Giovani for testing the earlier version of the app and provide some ideas to further develop it. Some `SEGUL` features are inspired by Phyluce pipelines.  -->
