# SEGUL

![Segul-Tests](https://github.com/hhandika/segul/workflows/Segul-Tests/badge.svg)
![Crate-IO](https://img.shields.io/crates/v/segul)
![GH-Release](https://img.shields.io/github/v/release/hhandika/segul)
![License](https://img.shields.io/github/license/hhandika/segul)

SEGUL is an ultrafast and memory efficient command-line application for alignment manipulation tasks that typically done using interpreted programming languages, such as Python, R, or Perl. It is designed to handle genomic datasets, but just as capable for Sanger datasets. In our test using a dataset with 4060 UCE loci, for instance, compare to a program written using biopython library, SEGUL is >40x faster for alignment concatenation while using 4x less RAM space.

Available features:

1. Converting alignments to different formats.
2. Concatenating alignments.
3. Filtering alignments based on minimal taxon completeness, alignment length, or numbers of parsimony informative sites.
4. Computing alignment summary statistics.
5. Getting sample IDs from a collection of alignments.

Planned features:

1. Filtering sequences from a collection of alignments based on user-defined IDs.
2. Converting dna sequences to amino acid and vice versa.

Supported sequence formats:

1. Nexus
2. Relaxed phylip
3. Fasta

All of the formats are supported in interleave and sequential.

It is now in active development. Our goal is to provide as many functionalities possible for alignment manipulation tasks.

## Supported Platforms

The program may work in any rust supported [platform](https://doc.rust-lang.org/nightly/rustc/platform-support.html). Below is a list of operating system that we tested and is guaranteed to work:

- Linux
- MacOS
- Windows
- Windows Sub-System for Linux

## Quick Start

If you already using a rust program and familiar with its [toolchain](https://www.rust-lang.org/learn/get-started), the best option is to install the app using cargo:

```{Bash}
cargo install segul
```

You can also use the pre-compiled binary available in [the release page](https://github.com/hhandika/segul/releases/). The installation is similar to any other single executable command line app, such as IQ-Tree or RaXML. All you need to do is to make sure the path to the app is registered in your environment variable, so that the app can be called from anywhere in your system. If you are not familiar with it, please see the detailed instruction below.

To program command structure is similar to git or gh-cli, and other program that use subcommands:

```{Bash}
segul <SUBCOMMAND> [OPTIONS] <VALUES> <A FLAG IF APPLICABLE>
```

For example to concat all the alignments in a directory named `nexus-alignments`:

```{Bash}
segul concat --dir nexus-alignments --format nexus
```

It is also available in short options:

```{Bash}
segul concat -d nexus-alignments -f nexus
```

To check for available subcommand:

```{Bash}
segul --help
```

To check for available options and flags for each sub-command:

```{Bash}
segul <SUBCOMMAND> --help
```

## Installation

We want the installation as flexible as possible. We offer three ways to install the program:

1. Using a pre-compiled library.
2. Through the rust package manager, cargo. It is similar to pip for python, germ for ruby, npm for JavaScript, etc.
3. Compiling it from the source available through Github repository.

### Using a pre-compiled library

This is the quickest and the most straigtforward option. The pre-compiled library is available at [the release page](https://github.com/hhandika/segul/releases/). You only need to download the zip file and extract it in your computer. Then, register it to your environmental variable.

For Linux/MacOS, this is an example of installing the app without leaving your terminal. We use version 0.3.1 as an example.

Donwload the library.

```{Bash}

wget https://github.com/hhandika/segul/releases/download/v0.3.1/segul-MacOS-x86_64.zip
```

Unzip the file.

```{Bash}
unzip segul-MacOS-x86_64.zip
```

Copy to a folder registered in your path variable.

### Installing through Cargo

This is the recommended option. Cargo will compile the program and fine-tuned it for your specific hardware. It also allows to easily updating the app. You will need the [rust compiler tool-chain](https://www.rust-lang.org/learn/get-started). It requires rust version 1.5 or higher. Then, check if the tool-chain installation successful:

```{Bash}
cargo --version
```

It should show the cargo version number. Then, install the app:

```{Bash}
cargo install segul
```

If you encounter a compiling  issue, you may need to install a C-development tool. For Debian-based Linux distribution, such as Debian, Ubuntu, PopOS, etc., the easiest way is to install build-essential:

```{Bash}
sudo apt install build-essential
```

For Windows, you only need to install the gnu toolchain for rust. The installation should be straighforward using rustup. Rustup comes as a part of the rust-compiler toolchain. It should be available in your system at the same time as you install cargo.

```{Bash}
rustup toolchain install stable-x86_64-pc-windows-gnu

rustup default stable-x86_64-pc-windows-gnu
```

To install the development version for any supported platform:

```{Bash}
cargo install --git https://github.com/hhandika/segul.git
```

## Usages

### Available command options

```{Bash}
USAGE:
    segul <SUBCOMMAND>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    concat     Concatenates alignments
    convert    Converts sequence formats
    filter     Filter alignments with specified min taxon completeness, alignment length, or parsimony informative sites
    help       Prints this message or the help of the given subcommand(s)
    id         Gets sample ids from multiple alignments
    summary    Gets alignment summary stats
```

More detailed instructions coming soon...
