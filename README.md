# SEGUL <img src="https://raw.githubusercontent.com/hhandika/segui/main/assets/launcher/iconDesktop.png" alt="segul logo" align="right" width="150"/>

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

- App Documentation: [[EN]](https://segul.app/)
- API Documentation: [[Rust]](https://docs.rs/segul/0.18.1/segul/)
- GUI source code: [[Repository]](https://github.com/hhandika/segui)
- Citation: [[Pre-print]](https://www.authorea.com/doi/full/10.22541/au.165167823.30911834/v1)

## Features

### Big changes in the version 0.19.0+ 💪🏼

- New command structure. Check [quick start](https://docs.page/hhandika/segul-docs/quick_start#cli-command-list) instruction to see the most up to date commands. We are working on updating the documentation throughout the website.
- New command to calculate summary statistics for raw reads and contigs.
- Input dir `-d` or `--dir` now can infer input format based on the file extension. No need to specify the input format anymore.
- Beta GUI version coming soon. We are cleaning up some bugs. Test alpha version [here](https://docs.page/hhandika/segul-docs/gui_install).

Bug fixes:

- Fix translation table errors.
- Update deprecated dependencies.
- Fix extract output issues.

### Full feature list

| Feature                        | Quick Link                                                                                                                                |
| ------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------- |
| Alignment concatenation        | [CLI](https://docs.page/hhandika/segul-docs/usage_concat) / [GUI](https://docs.page/hhandika/segul-docs/gui_usage_concat)                 |
| Alignment conversion           | [CLI](https://docs.page/hhandika/segul-docs/usage_convert) / [GUI](https://docs.page/hhandika/segul-docs/gui_usage_convert)               |
| Alignment filtering            | [CLI](https://docs.page/hhandika/segul-docs/usage_filter) / [GUI](https://docs.page/hhandika/segul-docs/gui_usage_filter)                 |
| Alignment splitting            | [CLI](https://docs.page/hhandika/segul-docs/usage_split) / [GUI](https://docs.page/hhandika/segul-docs/gui_usage_split)                   |
| Alignment partition conversion | [CLI](https://docs.page/hhandika/segul-docs/usage_part) / [GUI](https://docs.page/hhandika/segul-docs/gui_usage_part)                     |
| Alignment summary statistics   | [CLI](https://docs.page/hhandika/segul-docs/usage_summary) / [GUI](https://docs.page/hhandika/segul-docs/gui_usage_summary)               |
| Contig summary statistics      | [CLI](https://docs.page/hhandika/segul-docs/usage_contig_summary) / [GUI](https://docs.page/hhandika/segul-docs/gui_usage_contig_summary) |
| Raw read summary statistics    | [CLI](https://docs.page/hhandika/segul-docs/usage_raw_summary) / [GUI](https://docs.page/hhandika/segul-docs/gui_usage_raw_summary)       |
| Sample distribution mapping    | [CLI](https://docs.page/hhandika/segul-docs/usage_map) / [GUI](https://docs.page/hhandika/segul-docs/gui_usage_map)                       |
| Sequence extraction            | [CLI](https://docs.page/hhandika/segul-docs/usage_extract) / [GUI](https://docs.page/hhandika/segul-docs/gui_usage_extract)               |
| Sequence ID parsing            | [CLI](https://docs.page/hhandika/segul-docs/usage_id) / [GUI](https://docs.page/hhandika/segul-docs/gui_usage_id)                         |
| Sequence ID renaming           | [CLI](https://docs.page/hhandika/segul-docs/usage_rename) / [GUI](https://docs.page/hhandika/segul-docs/gui_usage_rename)                 |
| Sequence removal               | [CLI](https://docs.page/hhandika/segul-docs/usage_remove) / [GUI](https://docs.page/hhandika/segul-docs/gui_usage_remove)                 |
| Sequence translation           | [CLI](https://docs.page/hhandika/segul-docs/usage_translate) / [GUI](https://docs.page/hhandika/segul-docs/gui_usage_translate)           |

Supported sequence formats:

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
