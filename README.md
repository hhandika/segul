# segul

![Segul-Tests](https://github.com/hhandika/segul/workflows/Segul-Tests/badge.svg)

SEGUL is an ultrafast solution for alignment manipulation tasks that typically done using interpreted programming languages. For instance, compare to a program written using biopython library, segul is >30x faster for alignment concatenation while using 4x less ram space. SEGUL takes advantage of multi-core computers without extra effort from the users. It is guaranteed without issues of data races. The program is rigorously tested, manually or using an automated testing framework.

## Installation

Supported operating system:

- Linux
- MacOS
- Windows
- Windows Sub-System for Linux

This program is still under development. For now, it is available for testing only. You will need the [rust compiler tool-chain](https://www.rust-lang.org/learn/get-started) to install it. Then, use cargo to install the program:

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
