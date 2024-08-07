# SEGUL <img src="https://raw.githubusercontent.com/hhandika/segui/main/assets/launcher/iconDesktop.png" alt="segul logo" align="right" width="150"/>

![Segul-Tests](https://github.com/hhandika/segul/workflows/Segul-Tests/badge.svg)
![Crate-IO](https://img.shields.io/crates/v/segul)
![Crates-Download](https://img.shields.io/crates/d/segul?color=orange&label=crates.io-downloads)
![GH-Release](https://img.shields.io/github/v/tag/hhandika/segul?label=gh-releases)
![GH-Downloads](https://img.shields.io/github/downloads/hhandika/segul/total?color=blue&label=gh-release-downloads)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/segul/README.html)
![last-commit](https://img.shields.io/github/last-commit/hhandika/segul)
![License](https://img.shields.io/github/license/hhandika/segul)
[![LoC](https://tokei.rs/b1/github/hhandika/segul?category=code)](https://github.com/XAMPPRocky/tokei)

SEGUL is an ultra-fast, memory-efficient application for working with phylogenomic datasets. It is available as standalone, zero dependency command line, GUI applications (called SEGUI), and library/packages for Rust and other programming languages. It runs from your smartphone to High Performance Computers (see platform support below). It is safe, multi threaded, and easy to use.

It is designed to handle operations on large genomic datasets, while using minimal computational resources. However, it also provides convenient features for working on smaller datasets (e.g., Sanger datasets). In our tests, it consistently offers a faster and more efficient (low memory footprint) alternative to existing applications for a variety of sequence alignment manipulations ([see benchmark](https://www.segul.app/docs/cli_gui#performance)).

Learn more about SEGUL in the [documentation](https://www.segul.app/).

## Citation

> Handika, H., and J. A. Esselstyn. 2024. SEGUL: Ultrafast, memory-efficient and mobile-friendly software for manipulating and summarizing phylogenomic datasets. _Molecular Ecology Resources_. [https://doi.org/10.1111/1755-0998.13964](https://doi.org/10.1111/1755-0998.13964).

## Links

- App Documentation: [[EN]](https://segul.app/)
- API Documentation: [[Rust]](https://docs.rs/segul/0.18.1/segul/)
- GUI: [[Source code]](https://github.com/hhandika/segui)
- Python API: [[Source code]](https://github.com/hhandika/pysegul) [[PyPI]](https://pypi.org/project/pysegul/)

## Installation

### GUI Version

### Desktop

[<img style="padding-left: 15" alt="Microsoft Store download" src="https://get.microsoft.com/images/en-us%20dark.svg" width="200" />](https://apps.microsoft.com/detail/SEGUI/9NP1BQ6FW9PW?mode=direct)

[<img
    style="padding: 15"
    src="https://tools.applemediaservices.com/api/badges/download-on-the-mac-app-store/black/en-us?size=250x83&amp;releaseDate=1716076800"
    alt="Download on the Mac App Store"
    width="220"
  />](https://apps.apple.com/us/app/segui/id6447999874?mt=12&amp;itsct=apps_box_badge&amp;itscg=30200)

[![Get it from the Snap Store](https://snapcraft.io/static/images/badges/en/snap-store-black.svg)](https://snapcraft.io/segui)

### Mobile

[<img style="padding-left: 15" src="https://tools.applemediaservices.com/api/badges/download-on-the-app-store/black/en-us?size=250x83&amp;releaseDate=1716076800" alt="Download on the App Store" width="180">](https://apps.apple.com/us/app/segui/id6447999874?itsct=apps_box_badge&amp;itscg=30200)

[<img
    alt="Get it on Google Play"
    src="https://play.google.com/intl/en_us/badges/static/images/badges/en_badge_web_generic.png"
    height="80"
  />](https://play.google.com/store/apps/details?id=com.hhandika.segui&pcampaignid=pcampaignidMKT-Other-global-all-co-prtnr-py-PartBadge-Mar2515-1)

Learn more about device requirements and GUI app installation in the [documentation](https://www.segul.app/docs/installation/install_gui).

### CLI Version

The CLI app may work in any Rust supported [platform](https://doc.rust-lang.org/nightly/rustc/platform-support.html). However, we only tested and officially support the following platforms:

- Linux
- MacOS
- Windows
- Windows Subsystem for Linux (WSL)

#### CLI Installation Methods

- Pre-compiled binaries: [[Releases]](https://github.com/hhandika/segul/releases) [[Docs]](https://www.segul.app/docs/installation/install_binary)
- Conda: [[Package]](https://anaconda.org/bioconda/segul) [[Docs]](https://www.segul.app/docs/installation/install_conda)
- Package manager: [[Docs]](https://www.segul.app/docs/installation/install_cargo)
- From source: [[Docs]](https://www.segul.app/docs/installation/install_source)
- Beta version: [[Docs]](https://www.segul.app/docs/installation/install_dev)

### API version

The API version is available for Rust and other programming languages. For Rust users, you can install it via Cargo:

```bash
cargo add segul
```

#### Python

We provide binding for Python (called [pysegul](https://pypi.org/project/pysegul/)). Use SEGUL just like any other Python package:

```python
pip install pysegul
```

Learn more about using SEGUL API in the [documentation](https://www.segul.app/docs/api-usage/python/intro).

## Features

| Feature                        | Quick Link                                                                                                             |
| ------------------------------ | ---------------------------------------------------------------------------------------------------------------------- |
| Alignment concatenation        | [CLI](https://www.segul.app/docs/cli-usage/concat) / [GUI](https://www.segul.app/docs/gui-usage/align-concat)          |
| Alignment conversion           | [CLI](https://www.segul.app/docs/cli-usage/convert) / [GUI](https://www.segul.app/docs/gui-usage/align-convert)        |
| Alignment filtering            | [CLI](https://www.segul.app/docs/cli-usage/filter) / [GUI](https://www.segul.app/docs/gui-usage/align-filter)          |
| Alignment splitting            | [CLI](https://www.segul.app/docs/cli-usage/split) / [GUI](https://www.segul.app/docs/gui-usage/align-split)            |
| Alignment partition conversion | [CLI](https://www.segul.app/docs/cli-usage/part) / [GUI](https://www.segul.app/docs/gui-usage/align-partition)         |
| Alignment summary statistics   | [CLI](https://www.segul.app/docs/cli-usage/summary) / [GUI](https://www.segul.app/docs/gui-usage/align-summary)        |
| Genomic summary statistics     | [CLI](https://www.segul.app/docs/cli-usage/genomic) / [GUI](https://www.segul.app/docs/gui-usage/genomic)              |
| Sequence extraction            | [CLI](https://www.segul.app/docs/cli-usage/extract) / [GUI](https://www.segul.app/docs/gui-usage/sequence-extract)     |
| Sequence filtering             | [CLI](https://www.segul.app/docs/cli-usage/seq_filter) / GUI in development                                            |
| Sequence ID extraction         | [CLI](https://www.segul.app/docs/cli-usage/id) / [GUI](https://www.segul.app/docs/gui-usage/sequence-id)               |
| Sequence ID mapping            | [CLI](https://www.segul.app/docs/cli-usage/map) / [GUI](https://www.segul.app/docs/gui-usage/sequence-id-map)          |
| Sequence ID renaming           | [CLI](https://www.segul.app/docs/cli-usage/rename) / [GUI](https://www.segul.app/docs/gui-usage/sequence-rename)       |
| Sequence removal               | [CLI](https://www.segul.app/docs/cli-usage/remove) / [GUI](https://www.segul.app/docs/gui-usage/sequence-remove)       |
| Sequence translation           | [CLI](https://www.segul.app/docs/cli-usage/translate) / [GUI](https://www.segul.app/docs/gui-usage/sequence-translate) |

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

## Contribution

We welcome any kind of contribution, from issue reporting, ideas to improve the app, to code contribution. For ideas and issue reporting please post in [the Github issues page](https://github.com/hhandika/segul/issues). For code contribution, please fork the repository and send pull requests to this repository.
