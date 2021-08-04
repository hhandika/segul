# SEGUL

![Segul-Tests](https://github.com/hhandika/segul/workflows/Segul-Tests/badge.svg)
![Crate-IO](https://img.shields.io/crates/v/segul)
![GH-Release](https://img.shields.io/github/v/release/hhandika/segul)
![Download](https://img.shields.io/crates/d/segul?color=yellow)
![License](https://img.shields.io/github/license/hhandika/segul)

SEGUL is an ultrafast and memory efficient command-line (cli) application for working with sequence alignments that typically done using interpreted programming languages, such as Python, R, or Perl. It is designed to handle genomic datasets, but just as capable for Sanger datasets. In our test using a dataset with 4060 UCE loci, for instance, compare to a app written using biopython library, SEGUL is >40x faster for alignment concatenation while using 3x less RAM space.

Available features:

1. Converting alignments to different formats.
2. Concatenating alignments with partition settings.
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

All of the formats are supported in interleave and sequential. The app supports both DNA and amino acid sequences.

Supported partition formats:

1. RaXML
2. Embedded nexus (labeled `charset` in the app)
3. Separate nexus

It is now in active development. Our goal is to provide as many functionalities possible for alignment manipulation tasks.

Documentation: [GitHub Wiki](https://github.com/hhandika/segul/wiki)

## Supported Platforms

The app may work in any Rust supported [platform](https://doc.rust-lang.org/nightly/rustc/platform-support.html). Below is a list of operating system that we tested and is guaranteed to work:

- Linux
- MacOS
- Windows
- Windows Subsystem for Linux (WSL)

> Notes for Windows users: We recommend using [Windows Terminal](https://www.microsoft.com/en-us/p/windows-terminal/9n0dx20hk701#activetab=pivot:overviewtab) to take full advantages of the app features.

## Quick Start

If you are familiar with [the Rust programming language](https://en.wikipedia.org/wiki/Rust_(programming_language)) and its development [toolchain](https://www.rust-lang.org/learn/get-started), the best option is to install the app using the Rust package manager: [cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html). In addition to cargo, for Linux system (including WSL), it requires the C-development toolkit, `build-essential` for Debian-based distributions or its equivalent in other Linux distributions.

If you are new to using a command line application, installing through cargo is also the easiest route ([see details in the installation instruction](https://github.com/hhandika/segul/wiki/2.-Installation)). After you have cargo installed in your computer, installing SEGUL is one command away:

```Bash
cargo install segul
```

You can also use the pre-compiled binary available in [the release page](https://github.com/hhandika/segul/releases/). The installation is similar to any other single executable command line app, such as the phylogenetic programs IQ-Tree and RaXML. You only need to make sure the path to the app is registered in your environment variable, so that the app can be called from anywhere in your system ([see instructions](https://github.com/hhandika/segul/wiki/2.-Installation#using-a-pre-compiled-binary)).

The app command structure is similar to git, gh-cli, or any other app that use subcommands. The app file name will be `segul` for Linux/MacOS/WSL and `segul.exe` for Windows.

```Bash
[THE-PROGRAM-FILENAME] <SUBCOMMAND> [OPTIONS] <VALUES> <A-FLAG-IF-APPLICABLE>
```

To check for available subcommand:

```Bash
segul --help
```

To check for available options and flags for each sub-command:

```Bash
segul <SUBCOMMAND> --help
```

For example, to concat all the alignments in a directory named `nexus-alignments`:

```Bash
segul concat --dir nexus-alignments --input-format nexus
```

It is also available in short options:

```Bash
segul concat -d nexus-alignments -f nexus
```

The app outputs are the resulting files from each task and a log file. Most information that is printed to the terminal is written to the log file. Unlike the terminal output that we try to keep it clean and only show the most important information, the log file will also contain the dates, times, and the log level status. Each time you run the app, the app will append the log output to the same log file (named `segul.log`) if the file exists in the same directory. Rename this file or move it to a different folder if you would like to keep a different log file for each task.

<!-- ## Acknowledgment -->
