# SEGUL

![Segul-Tests](https://github.com/hhandika/segul/workflows/Segul-Tests/badge.svg)
![Crate-IO](https://img.shields.io/crates/v/segul)
![Crates-Download](https://img.shields.io/crates/d/segul?color=orange&label=crates.io-downloads)
![GH-Release](https://img.shields.io/github/v/release/hhandika/segul)
![GH-Downloads](https://img.shields.io/github/downloads/hhandika/segul/total?color=blue&label=binary-downloads)
[![LoC](https://tokei.rs/b1/github/hhandika/segul?category=code)](https://github.com/XAMPPRocky/tokei)
![last-commit](https://img.shields.io/github/last-commit/hhandika/segul)
![License](https://img.shields.io/github/license/hhandika/segul)

SEGUL is an ultrafast and memory efficient command-line (cli) application for working with sequence alignments that is typically done using interpreted programming languages, such as Python, R, or Perl. It is a cross-platform single executable app and has zero runtime dependency on MacOS and Windows and only relies on a library provided by the OS on Linux. It is designed to handle genomic datasets, but just as capable for Sanger datasets. In our test using a dataset with 4060 UCE loci, for instance, compare to an app written using [the biopython library](https://biopython.org/), SEGUL is >40x faster for alignment concatenation while using 3x less RAM space.

Features:

1. Converting alignments to different formats.
2. Concatenating alignments with partition settings.
3. Splitting alignments by partitions.
4. Filtering alignments based on minimal taxon completeness, alignment length, or numbers of parsimony informative sites.
5. Computing alignment summary statistics.
6. Getting sample IDs from a collection of alignments.
7. Map sample distribution in a collection of alignments.
8. Extracting sequences from a collection of alignments based on user-defined IDs. The feature comes with regular expression support.
9. Batch renaming sequence IDs.
10. Converting dna sequences to amino acid.

Supported sequence formats:

1. Nexus
2. Relaxed phylip
3. Fasta

All of the formats are supported in interleave and sequential. The app supports both DNA and amino acid sequences.

Supported partition formats:

1. RaXML
2. Nexus

The Nexus partition can be written as a charset block embedded in Nexus formatted sequences or be written in a separate file.

Documentation: [GitHub Wiki](https://github.com/hhandika/segul/wiki)

Citation:

> Handika, H. and Esselstyn, J. A. In prep. SEGUL: An ultrafast, memory efficient, and cross-platform alignment manipulation tool for phylogenomics.

## Quick Links

- [SEGUL](#segul)
  - [Quick Links](#quick-links)
  - [Supported Platforms](#supported-platforms)
  - [Quick Start](#quick-start)
    - [Installation](#installation)
      - [Using pre-compiled binary](#using-pre-compiled-binary)
      - [Install from the Rust Package Manager](#install-from-the-rust-package-manager)
    - [Command Structure](#command-structure)
    - [Input Options](#input-options)
    - [Datatype](#datatype)
    - [Output](#output)
    - [Converting alignments](#converting-alignments)
    - [Concatenating alignments](#concatenating-alignments)
    - [Splitting alignments by partitions](#splitting-alignments-by-partitions)
    - [Computing sequence summary statistics](#computing-sequence-summary-statistics)
    - [Getting sample IDs from a collection of alignments](#getting-sample-ids-from-a-collection-of-alignments)
    - [Map sample distribution in a collection of alignments](#map-sample-distribution-in-a-collection-of-alignments)
    - [Filtering alignments](#filtering-alignments)
    - [Extracting sequences in alignments](#extracting-sequences-in-alignments)
    - [Batch renaming sequence IDs](#batch-renaming-sequence-ids)
    - [Translating DNA sequences](#translating-dna-sequences)
  - [Logging](#logging)
  - [Contribution](#contribution)

## Supported Platforms

The app may work in any Rust supported [platform](https://doc.rust-lang.org/nightly/rustc/platform-support.html). Below is a list of operating system that we tested and is guaranteed to work:

- Linux
- MacOS
- Windows
- Windows Subsystem for Linux (WSL)

> :warning: **SEGUL modern terminal output comes with a cost of requiring a terminal application that supports [UTF-8 encoding](https://en.wikipedia.org/wiki/UTF-8)**. For MacOS and native Linux, your default terminal should have supported UTF-8 encoding by default. For Windows (including WSL) users, we recommend using [Windows Terminal](https://www.microsoft.com/en-us/p/windows-terminal/9n0dx20hk701#activetab=pivot:overviewtab) to ensure consistent terminal output. Windows Terminal requires separate installation for Windows 10. It should come pre-installed on Windows 11.

## Quick Start

The instruction below assumes familiarity with command line application and only highlight some common features that users may need for alignment manipulation and generating sequence statistics tasks. We provide more detailed instruction in the [documentation](https://github.com/hhandika/segul/wiki).

### Installation

#### Using pre-compiled binary

For a quick installation, we provide pre-compiled binaries in [the release page](https://github.com/hhandika/segul/releases/). For WSL, either the ManyLinux or Linux binary should work. In our test system, the ManyLinux binary is a little faster. For native Linux OS, first check your GLIBC version:

```Bash
ldd --version
```

If your system GLIBC is >=2.18, use the Linux binary. If lower, use the ManyLinux binary. The installation is similar to any other single executable command-line app, such as the phylogenetic programs IQ-Tree or RaXML. You only need to make sure the path to the app is registered in your environment variable, so that the app can be called from anywhere in your system ([see instructions](https://github.com/hhandika/segul/wiki/2.-Installation#using-a-pre-compiled-binary)). If you are still having issues running the program, try to install it using the package manager. This installation method will optimize the compiled binary for your system (see below).

> ATTENTION!: For MacOS users, when you run `segul` for the first time, MacOS gatekeeper will prevent the program to run because `segul` is not signed by Apple. Go to Security Setting to allow `segul` to run. More details are here in [Apple Website](https://support.apple.com/en-us/HT202491).

#### Install from the Rust Package Manager

The Rust package manager is called [cargo](https://crates.io/). Cargo is easy to install (also easy to uninstall) and will help you to manage the app ([see details in the installation instruction](https://github.com/hhandika/segul/wiki/2.-Installation)). Installing SEGUL through Cargo is similar to installing it from source code, except that it only use the stable version of the code. The source code is managed on [crates.io](https://crates.io/). The badge at top of this Readme has information on the latest version of the app available on [crates.io](https://crates.io/).

After you have Cargo installed in your computer, in Linux system (including WSL), first install the C-development toolkit, `build-essential` for Debian-based distributions (Debian, Ubuntu, PopOS, Linux Mint, etc.) or its equivalent in other Linux distributions:

```Bash
sudo apt install build-essential
```

On Windows, you only need to install the GNU compiler toolchain available using Rustup. Rustup is installed automatically when you install Cargo. To install the toolchain:

```Bash
rustup toolchain install stable-x86_64-pc-windows-gnu

rustup default stable-x86_64-pc-windows-gnu
```

Then, install SEGUL:

```Bash
cargo install segul
```

You could also install SEGUL from the GitHub repository. Learn more about SEGUL installation [here](https://github.com/hhandika/segul/wiki/2.-Installation).

### Command Structure

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

Across the app functions, most generic arguments are also available in short format to save time typing them. For example, below we use short arguments to concat alignments in a directory named `nexus-alignments`:

```Bash
segul concat -d nexus-alignments -f nexus
```

Learn more about SEGUL command structure and expected behaviors for each argument [here](https://github.com/hhandika/segul/wiki/4.-Command-Structure).

### Input Options

The app has two input options. The standard input `--input` or `-i` in short format and `--dir` or `-d` in short format. If your input files are all in a single directory, you should use the `--dir` or `-d` option and specify the file format:

```Bash
segul <SUBCOMMAND> -d alignment_dir -f nexus
```

When dealing with a single file, more complex folder structure, or unusual file extensions, use the `--input` or `-i` option.

For a single file:

```Bash
segul <SUBCOMMAND> -i alignment-dir/alignment_file.fasta
```

Multiple file in a directory using wildcard:

```Bash
segul <SUBCOMMAND> -i alignment-dir/*.fasta
```

Multiple files in multiple directories:

```Bash
segul <SUBCOMMAND> -i alignment-dir1/*.fasta alignment-dir2/*.fasta
```

For unusual file extensions or if the app failed to detect the file format, specify the input format:

```Bash
segul <SUBCOMMAND> -i alignment-dir/*.aln -f fasta
```

Both of the input options are available in all subcommands. To keep it simple, the command examples below use `--dir` as an input.

### Datatype

The app support both DNA and amino acid sequences. It will check whether the sequences contain only valid IUPAC characters of the datatype. By default, it sets to DNA sequences. Use the option `--datatype aa` if your input is amino acid sequences. For example:

```Bash
segul convert -d /alignments -f nexus --datatype aa
```

### Output

Most functions will save into their default directory. For example, the concat function will default to create `SEGUL-concat` directory and will save its output files into the directory. To specify the output directory, use the `--output` or `-o` option. For example:

```Bash
segul convert -d /alignments -f nexus -o alignments_concat
```

The app avoids over-writting files with similar names. The app will check if a such file or directory exists and will ask if you like to remove it. The app will exit if you decide to not remove it.

### Converting alignments

Segul can convert a single sequence file or multiple sequence files in a directory:

```Bash
segul convert --dir [path-to-your-repository] --input-format [sequence-format-keyword] --output-format [sequence-format-keyword]
```

### Concatenating alignments

To concat all alignments in a directory:

```Bash
segul concat --dir [a-path-to-a-directory] --input-format [sequence-format-keyword]
```

### Splitting alignments by partitions

To split alignment by partions, you need the alignment file and the alignment partion in a separate file:

```Bash
segul split -i [a-path-to-an-alignment] --input-partition [a-path-to-partition-file]
```

If it is not provided, `segul` will use the alignment name as an output directory. To provide the output directory name, use the `--output` or `-o` option.

### Computing sequence summary statistics

To generate sequence summary statistics of alignments in a directory:

```Bash
segul summary --dir [a-path-to-a-directory] --input-format [sequence-format-keyword]
```

### Getting sample IDs from a collection of alignments

You have multiple alignments and want to know what are samples you have in all of those alignment. You can easily do it using `segul`. The app can find all the unique IDs across thousands of alignments within seconds.

```Bash
segul id --dir [a-path-to-a-directory] --input-format [sequence-format-keyword]
```

It will generate a text file that contains all the unique IDs across your alignments.

### Map sample distribution in a collection of alignments

If you would like to know how the samples distributed across your alignments, you only need to add the `--map` flag when searching for unique IDs. It will generate both the unique IDs (in a text file) and the sample distribution (in csv).

```Bash
segul id --dir [a-path-to-a-directory] --input-format [sequence-format-keyword] --map
```

### Filtering alignments

Segul provide multiple filtering parameters.

```Bash
segul filter --dir [a-path-to-a-directory] --input-format [sequence-format-keyword] <parameters>
```

For example, to filter based on taxon completeness:

```Bash
segul filter --dir [a-path-to-a-directory] --input-format [sequence-format-keyword] --percent [percentages-of-minimal-taxa]
```

Other available parameters are multiple minimal taxon completeness `--npercent`, alignment length `--len`, numbers of minimal parsimony informative sites `--pinf`, and percent of minimal parsimony informative sites `--percent-inf`.

By default, the app will copy files that are match with the parameter to a new folder. If you would like to concat the results instead, you can specify it by passing `--concat` flags. All the options available for the concat function above also available for concatenating filtered alignments.

### Extracting sequences in alignments

You can also extract sequences from a collection of alignments. It can be done by supplying a list of IDs directly on the command line or in text file. The app finds for the exact match. You can also use regular expression to search for matching IDs.

To extract sequences by inputing the IDs in the command line:

```bash
segul extract --dir [path-to-alignment-directory] --input-format [sequence-format-keyword] --id [id_1] [id_2] [id_3]
```

You can specify as many id as you would like. However, for long list of IDs, it may be better to input it using a text file. In the file it should be only the ID list, one ID each line:

```bash
sequence_1
sequence_2
sequence_3
sequence_4
```

The the command will be:

```bash
segul extract --dir [path-to-alignment-directory] --input-format [sequence-format-keyword] --file [path-to-text-file]
```

For using regular expression:

```bash
segul extract -d gblock_trimmed_80p/ -f nexus --re="regex-syntax"
```

The app uses the rust [regex library](https://docs.rs/regex/1.5.4/regex/) to parse regular expression. The syntax is similar to Perl regular expression (find out more [here](https://docs.rs/regex/1.5.4/regex/)).

### Batch renaming sequence IDs

To rename sequence IDs in multiple alignments, you need to input the sequence IDs in tsv or csv format with header. For example:

| Original_names        | New_names                |
| --------------------- | ------------------------ |
| Genus_species1_random | Genus_species1_voucherID |
| Genus_species2_random | Genus_species2_voucherID |

To simplify this process, you can generate unique IDs for all of your alignments using the `id` sub-command.

```Bash
segul id -d [alignment-dir] -f [sequence-format-keyword]
```

Copy the IDs to Excel and then write a new names and the header names as above. Save the file as csv or tsv. The program will infer the file format based on the file extension. Use it as an `--names` or `-n` input for renaming the sequence IDs using `rename` sub-command:

```Bash
segul rename -d [alignment-dir] -f [sequence-format-keyword] -n [file-path-to-IDs]
```

Example:

```Bash
segul rename -d uce-loci/ -f nexus -n new_names.csv
```

You can also change the output format by using `--output-format` or `-F` option.

### Translating DNA sequences

List of supported [NCBI Genetic Code Tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#top) is available [here](https://github.com/hhandika/segul/wiki/5.-Usages#translating-dna-sequences).

To translate dna alignment to amino acid:

```Bash
segul translate -d [path-to-alignment-files] -f [sequence-format-keyword]
```

By default, the app will use the standard code table (NCBI Table 1). To set the translation table, use the `--table` option. For example, to translate dna sequences using NCBI Table 2 (vertebrate MtDNA):

```Bash
segul translate -d loci/ -f fasta --table 2
```

You can also set the reading frame using the `--rf` option:

```Bash
segul translate -d loci/ -f fasta --table 2 --rf 2
```

To show all the table options, use the `--show-tables` flag:

```Bash
segul translate --show-tables
```

## Logging

Most information that is printed to the terminal is written to the log file (named `segul.log`). It is written to the current working directoy. Unlike the terminal output that we try to keep it clean and only show the most important information, the log file will also contain the dates, times, and the log level status. Each time you run the app, if the log file exists in the same directory, the app will append the log output to the same log file. Rename this file or move it to a different directory if you would like to keep a different log file for each task.

Learn more about using SEGUL [here](https://github.com/hhandika/segul/wiki/5.-Usages).

## Contribution

We welcome any kind of contribution, from issue reporting, ideas to improve the app, to code contribution. For ideas and issue reporting please post in [the Github issues page](https://github.com/hhandika/segul/issues). For code contribution, please fork the repository and send pull requests to this repo

<!-- ## Acknowledgment

We thank Giovani for testing the earlier version of the app and provide some ideas to further develop it. Some `SEGUL` features are inspired by Phyluce pipelines.  -->
