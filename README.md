# SEGUL <img src="https://raw.githubusercontent.com/hhandika/segui/main/assets/images/logo.png" alt="segul logo" align="right" width="150"/>

![Segul-Tests](https://github.com/hhandika/segul/workflows/Segul-Tests/badge.svg)
![Crate-IO](https://img.shields.io/crates/v/segul)
![Crates-Download](https://img.shields.io/crates/d/segul?color=orange&label=crates.io-downloads)
![GH-Release](https://img.shields.io/github/v/tag/hhandika/segul?label=gh-releases)
![GH-Downloads](https://img.shields.io/github/downloads/hhandika/segul/total?color=blue&label=gh-release-downloads)
[![LoC](https://tokei.rs/b1/github/hhandika/segul?category=code)](https://github.com/XAMPPRocky/tokei)
![last-commit](https://img.shields.io/github/last-commit/hhandika/segul)
![License](https://img.shields.io/github/license/hhandika/segul)

SEGUL is an ultra-fast, memory-efficient application for working with phylogenomic datasets. It is available as standalone, zero dependency command line, GUI applications (called SEGUI), and library/packages for Rust and other programming languages. It runs from your smartphone to High Performance Computers (see platform support below). It is safe, multi threaded, and easy to use.

It is designed to handle operations on large genomic datasets, while using minimal computational resources. However, it also provides convenient features for working on smaller datasets (e.g., Sanger datasets). In our tests, it consistently offers a faster and more efficient (low memory footprint) alternative to existing applications for a variety of sequence alignment manipulations ([see benchmark](https://github.com/hhandika/segul-bench)).

## Links

- App Documentation: [[EN]](https://docs.page/hhandika/segul-docs/)
- Source Code Documentation: [[EN]](https://docs.rs/segul/0.18.1/segul/)
- GUI source code: [[LINK]](https://github.com/hhandika/segui)
- Citation: [pre-print](https://doi.org/10.22541/au.165167823.30911834/v1)

## Features

| Feature                          | Quick Link                                               |
| -------------------------------- | -------------------------------------------------------- |
| Alignment concatenation          | [[CLI]](https://docs.page/hhandika/segul-docs/usage_concat) [[GUI]]()|
| Alignment conversion             | [[CLI]]() [[GUI]]()                                                         |
| Alignment filtering              | [[CLI]]() [[GUI]]()                                                         |
| Alignment splitting              | [[CLI]]() [[GUI]]()                                                         |
| Alignment summary statistics     | [[CLI]]() [[GUI]]()                                                         |
| Contiguous sequence statistics   | [[CLI]]() [[GUI]]()                                                         |
| Partition format conversion      | [[CLI]]() [[GUI]]()                                                         |
| Sequence ID extraction           | [[CLI]]() [[GUI]]()                                                         |
| Sequence ID renaming             | [[CLI]]() [[GUI]]()                                                         |
| Sequence read summary statistics | [[CLI]]() [[GUI]]()                                                         |
| Sequence removal                 | [[CLI]]() [[GUI]]()                                                         |
| Sequence translation             | [[CLI]]() [[GUI]]()                                                         |

Supported common sequence formats:

1. NEXUS
2. Relaxed PHYLIP
3. FASTA
4. FASTQ (gzipped and uncompressed)

All of the formats are supported in interleave and sequential versions. Except for FASTQ (DNA only), the app supports both DNA and amino acid sequences.

Supported partition formats:

1. RaXML
2. NEXUS

The NEXUS partition can be written as a charset block embedded in NEXUS formatted sequences or a separate file.

## Supported Platforms

The CLI app may work in any Rust supported [platform](https://doc.rust-lang.org/nightly/rustc/platform-support.html). For both CLI and GUI, below is a list of operating system that we tested and is guaranteed to work:

- Linux (CLI)
- MacOS (GUI and CLI)
- Windows (GUI and CLI)
- Windows Subsystem for Linux (WSL) (CLI)
- iOS (GUI)
- iPadOS (GUI)
- Android (GUI and CLI using [Termux](https://termux.com/))

> For Windows (including WSL) users, we recommend using [Windows Terminal](https://www.microsoft.com/en-us/p/windows-terminal/9n0dx20hk701#activetab=pivot:overviewtab) to ensure consistent terminal output. Windows Terminal requires separate installation for Windows 10. It should come pre-installed on Windows 11.

## Contribution

We welcome any kind of contribution, from issue reporting, ideas to improve the app, to code contribution. For ideas and issue reporting please post in [the Github issues page](https://github.com/hhandika/segul/issues). For code contribution, please fork the repository and send pull requests to this repository.
