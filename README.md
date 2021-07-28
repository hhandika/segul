# SEGUL

![Segul-Tests](https://github.com/hhandika/segul/workflows/Segul-Tests/badge.svg)

SEGUL is an ultrafast and efficient solution for alignment manipulation tasks that typically done using interpreted programming languages, such as Python, R, and Perl. It is designed to handle genomic datasets, but just as capable for Sanger datasets. In our test using a dataset with 4060 UCE loci, for instance, compare to a program written using biopython library, SEGUL is >40x faster for alignment concatenation while using 4x less RAM space.

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

All formats are supported in interleaved and sequential forms.

The program is now in active development. Our goal is to provide as many functionalities possible for alignment manipulation tasks.

## Installation

Supported operating system:

- Linux
- MacOS
- Windows
- Windows Sub-System for Linux

This program is still under development. For now, it is available for testing only. You will need the [rust compiler tool-chain](https://www.rust-lang.org/learn/get-started) to install SEGUL. It requires rust version 1.5 or higher. Then, check if the tool-chain installation successful:

```{Bash}
cargo --version
```

It should show the cargo version number.

For Linux, install build-essential if you don't have it install yet. For Debian based distros, such as Debian, Ubuntu, PopOS, Linux Mint, etc:

```{Bash}
sudo apt install build-essential
```

For windows users:

```{Bash}
rustup toolchain install stable-x86_64-pc-windows-gnu

# then
rustup default stable-x86_64-pc-windows-gnu
```

Then, install the program:

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
