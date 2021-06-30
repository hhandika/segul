# SEGUL

![Segul-Tests](https://github.com/hhandika/segul/workflows/Segul-Tests/badge.svg)

SEGUL is an ultrafast and efficient solution for alignment manipulation tasks that typically done using interpreted programming languages, such as Python, R, and Perl. It is designed to handle genomic datasets, but just as capable for Sanger datasets. In our test using a dataset with 4060 UCE loci, for instance, compare to a program written using biopython library, SEGUL is >30x faster for alignment concatenation while using 4x less RAM space. 

SEGUL brings the benefit of the Rust programming languages that guarantees no invalid RAM access and multi-core performance without data races. There is no extra effort from the users to specify the numbers of cores. The program will assess the available system resources and use multi-cores when the job requires them. It only uses as many cores that the task needs. This avoids the common problem that users over estimates the number of cores required for a task, which often lead to the program being slower or even worse having data race issues. The program is rigorously tested, manually or using an automated testing framework.

Several tasks that currently SEGUL can do:

1. Converting alignments to different formats.
2. Concatenating alignments.
3. Filtering alignments based on the taxon completeness.
4. Computing alignment summary statistics.
5. Getting sample IDs from a collection of alignments.

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

This program is still under development. For now, it is available for testing only. You will need the [rust compiler tool-chain](https://www.rust-lang.org/learn/get-started) to install SEGUL. It requires rust version 1.5 or higher. Then, check if the tool-chain installation succesful:

```{Bash}
cargo --version
```
It should show the cargo version number. Then, install the program:

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
    help       Prints this message or the help of the given subcommand(s)
    id         Gets unique ids from multiple alignments
    pick       Picks alignments with specified min taxa
    summary    Gets alignment summary stats
```
