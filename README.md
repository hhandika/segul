# SEGUL

![Segul-Tests](https://github.com/hhandika/segul/workflows/Segul-Tests/badge.svg)
![Crate-IO](https://img.shields.io/crates/v/segul)
![Crates-Download](https://img.shields.io/crates/d/segul?color=orange&label=crates.io-downloads)
![GH-Release](https://img.shields.io/github/v/tag/hhandika/segul?label=gh-releases)
![GH-Downloads](https://img.shields.io/github/downloads/hhandika/segul/total?color=blue&label=gh-release-downloads)
[![LoC](https://tokei.rs/b1/github/hhandika/segul?category=code)](https://github.com/XAMPPRocky/tokei)
![last-commit](https://img.shields.io/github/last-commit/hhandika/segul)
![License](https://img.shields.io/github/license/hhandika/segul)

SEGUL is an ultrafast, memory-efficient command-line (cli) application for working with sequence alignments. It is a cross-platform, single executable app, with zero runtime dependency on MacOS and Windows. On Linux, it relies on only a library provided by the OS. In simple words, you need not worry about dependencies. It is designed to handle operations on large genomic datasets, while using minimal computational resources. However, it also provides convenient features for working on smaller datasets (e.g., Sanger datasets). In our tests, it consistently offers a faster and more efficient (low memory footprint) alternative to existing applications for a variety of sequence alignment manipulations ([see benchmark](https://github.com/hhandika/segul-bench)).

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

<!-- Upcoming features and bug fixes ([try](https://github.com/hhandika/segul/wiki/2.-Installation#try-development-features)):

1. A feature to generate summary sequence statistics for each locus.
2. Taxon summary statistics now include proportion of missing data, nucleotide counts, GC content and AT content. -->

Supported sequence formats:

1. NEXUS
2. Relaxed PHYLIP
3. FASTA

All of the formats are supported in interleave and sequential versions. The app supports both DNA and amino acid sequences.

Supported partition formats:

1. RaXML
2. NEXUS

The NEXUS partition can be written as a charset block embedded in NEXUS formatted sequences or a separate file.

Documentation: [GitHub Wiki](https://github.com/hhandika/segul/wiki)

Citation: [pre-print](https://doi.org/10.22541/au.165167823.30911834/v1)

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
    - [Usages](#usages)
      - [Computing sequence summary statistics](#computing-sequence-summary-statistics)
      - [Concatenating alignments](#concatenating-alignments)
      - [Converting sequence formats](#converting-sequence-formats)
      - [Converting partition formats](#converting-partition-formats)
      - [Extracting sequences in alignments](#extracting-sequences-in-alignments)
      - [Filtering alignments](#filtering-alignments)
      - [Getting sample IDs from a collection of alignments](#getting-sample-ids-from-a-collection-of-alignments)
      - [Map sample distribution in a collection of alignments](#map-sample-distribution-in-a-collection-of-alignments)
      - [Splitting concatenated alignments by partitions](#splitting-concatenated-alignments-by-partitions)
      - [Batch removing sequence based on user-defined IDs](#batch-removing-sequence-based-on-user-defined-ids)
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

We offer multiple options to install `segul`. The easiest and quickest way is using the pre-compiled binary. You won't have to worry about having the Rust compiler tool-chain installed in your system. If you are already using applications written in Rust, or if none of the available compiled binaries work for your system, we recommend installing `segul` using the Rust Package Manager.

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

Learn more about SEGUL command structure and expected behaviors for each argument [here](https://github.com/hhandika/segul/wiki/4.-Command-Options).

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

Both of the input options are available in all subcommands. Learn more about input options ([here](https://github.com/hhandika/segul/wiki/4.-Command-Options#input-options)).

### Datatype

The app support both DNA and amino acid sequences. It will check whether the sequences contain only valid IUPAC characters of the datatype. By default, it sets to DNA sequences. Use the option `--datatype aa` if your input is amino acid sequences. For example:

```Bash
segul convert -d /alignments -f nexus --datatype aa
```

Learn more about supported data type [here](https://github.com/hhandika/segul/wiki/4.-Command-Options#data-types).

### Output

Most functions will save into their default directory. For example, the concat function will default to create `SEGUL-concat` directory and will save its output files into the directory. To specify the output directory, use the `--output` or `-o` option. For example:

```Bash
segul convert -d /alignments -f nexus -o alignments_concat
```

By default, the app avoids over-writing files with similar names. The app will check if a such file or directory exists and will ask if you would like to remove it. The app will exit if you decide to not remove it. You can pass `--overwrite` flag to remove existing output files or directories without asking first.

Learn more about specifying the output [here](https://github.com/hhandika/segul/wiki/4.-Command-Options#output).

### Usages

#### Computing sequence summary statistics

`segul` generate summary statistics for all alignment, each locus, and each taxon. To generate the summary statistics in for alignment files in a directory (usign `--dir` or `-d` option):

```Bash
segul summary --dir [a-path-to-a-directory] --input-format [sequence-format-keyword]
```

Using `--input` (or `-i` in short version) option ([learn the difference](https://github.com/hhandika/segul/wiki/4.-Command-Options#input-options)):

```Bash
segul summary --input [a-path/wilcard-to-files]
```

Learn more about computing sequence summary statistics [here](https://github.com/hhandika/segul/wiki/5.-Usages#computing-sequence-summary-statistics).

#### Concatenating alignments

To concat all alignments in a directory using `--dir` (or `-d` in short version):

```Bash
segul concat --dir [a-path-to-a-directory] --input-format [sequence-format-keyword]
```

Using `--input` (or `-i` in short version) option ([learn the difference](https://github.com/hhandika/segul/wiki/4.-Command-Options#input-options)):

```Bash
segul concat --input [path/wildcard-to-your-files]
```

Learn more about concatenating alignments [here](https://github.com/hhandika/segul/wiki/5.-Usages#concatenating-alignments).

#### Converting sequence formats

Segul can convert a single sequence file or multiple sequence files in a directory. To input the file, use `--dir` or `-d` option:

```Bash
segul convert --dir [path-to-your-repository] --input-format [sequence-format-keyword] --output-format [sequence-format-keyword]
```

Using `--input` (or `-i` in short version) option ([learn the difference](https://github.com/hhandika/segul/wiki/4.-Command-Options#input-options)):

```Bash
segul convert --input [path/wildcard-to-your-files] --output-format [sequence-format-keyword]
```

Learn more about converting sequence formats [here](https://github.com/hhandika/segul/wiki/5.-Usages#converting-sequences-to-a-different-format).

#### Converting partition formats

`segul` can convert a single and multiple partition files in multiple folders. It can extract partition embedded in NEXUS sequence files or merge codon model partition. For this command, the input option is only available using the `--input` option (or `-i` in short version).

```Bash
segul partition --input <a-path/wildcard-to-partition> --input-part <input-partition-format> --output-part<output-partition-format>
```

For example, to extract nexus in-file partitions (called charset format in `segul`):

```Bash
segul partition --input concatenated_alignment.nex --input-part charset --output-part nexus
```

You can also use wildcard to convert multiple concatenated alignment partitions at once:

```Bash
segul partition --input ./*/concatenated_alignment.nex --input-part charset --output-part nexus
```

When merging codon model partition, `segul` can clean up locus names if they are suffixes with subset{codon_pos}, subset\_{codon_pos}, or {codon_pos}\_pos. For example, for this partition format:

```Text
DNA, locus_1_subset1 = 1-10\3
DNA, locus_1_subest2 = 2-10\3
DNA, locus_1_subset3 = 3-10\3
DNA, locus_2_subset1 = 11-20\3
DNA, locus_2_subest2 = 12-20\3
DNA, locus_2_subset3 = 13-20\3
```

or

```Text
DNA, locus_1_1stpos = 1-10\3
DNA, locus_1_2ndpos = 2-10\3
DNA, locus_1_3rdpos = 3-10\3
DNA, locus_2_1stpos = 11-20\3
DNA, locus_2_2ndpos = 12-20\3
DNA, locus_2_3rdpos = 13-20\3
```

Will be merge as:

```Text
DNA, locus_1 = 1-10
DNA, locus_2 = 11-20
```

`segul` will save the resulting partition file in the current directory and the output file name will be the input file suffixed it with `_partition`. For example, if the input file is `alignment.nexus`, the output filename will be `alignment_partition.nex` or`alignment_partition.txt` (if the output is in a RaxML format)

Learn more about converting partition formats [here](https://github.com/hhandika/segul/wiki/5.-Usages#converting-partition-formats).

#### Extracting sequences in alignments

You can also extract sequences from a collection of alignments. It can be done by supplying a list of IDs directly using the command line or in a text file. `segul` will look for the exact match of the IDs and is case sensitive. You can also use regular expression for more flexiple option.

To extract sequences by inputing the IDs in the command line:

```bash
segul extract --dir [path-to-alignment-directory] --input-format [sequence-format-keyword] --id [id_1] [id_2] [id_3]
```

Using `--input` (or `-i` in short version) option ([learn the difference](https://github.com/hhandika/segul/wiki/4.-Command-Options#input-options)):

```bash
segul extract --input [a-path/wilcard-to-files] --id [id_1] [id_2] [id_3]
```

You can specify as many id as you would like. However, for a long list of IDs, it may be more convenient to input a text file. In the file, the content should be only the list of IDs, one ID each line:

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

`segul` uses the rust [regex library](https://docs.rs/regex/1.5.4/regex/) to parse regular expression. The syntax is similar to Perl regular expression (find out more [here](https://docs.rs/regex/1.5.4/regex/)).

Learn more about extracting sequences [here](https://github.com/hhandika/segul/wiki/5.-Usages#extracting-sequences-in-alignments).

#### Filtering alignments

Segul provide multiple filtering parameters (see below).

Using a directory `--dir` (or `-d` in short version):

```Bash
segul filter --dir [a-path-to-a-directory] --input-format [sequence-format-keyword] <parameters>
```

Using `--input` (or `-i` in short version) option ([learn the difference](https://github.com/hhandika/segul/wiki/4.-Command-Options#input-options)):

```Bash
segul filter --input [a-path/wilcard-to-files] <parameters>
```

For example, to filter based on taxon completeness:

```Bash
segul filter --dir [a-path-to-a-directory] --input-format [sequence-format-keyword] --percent [percentages-of-minimal-taxa]
```

Other available parameters are multiple minimal taxon completeness `--npercent`, alignment length `--len`, numbers of minimal parsimony informative sites `--pinf`, and percent of minimal parsimony informative sites `--percent-inf`.

By default, the app will copy files that are match with the parameter to a new folder. If you would like to concat the results instead, you can specify it by passing `--concat` flags. All the options available for the concat function [above](#concatenating-alignments) also available for concatenating filtered alignments.

Learn more about filtering alignments [here](https://github.com/hhandika/segul/wiki/5.-Usages#filtering-alignments).

#### Getting sample IDs from a collection of alignments

You have multiple alignments and want to know what are samples you have in all of those alignment. You can easily do it using `segul`. The app can find all the unique IDs across thousands of alignments within seconds.

The command using a `--dir` (or `-d` in short version):

```Bash
segul id --dir [a-path-to-a-directory] --input-format [sequence-format-keyword]
```

Using `--input` (or `-i` in short version) option ([learn the difference](https://github.com/hhandika/segul/wiki/4.-Command-Options#input-options)):

```Bash
segul id --input [a-path/wilcard-to-files]
```

It will generate a text file that contains all the unique IDs across your alignments.

Learn more about getting sample IDs [here](https://github.com/hhandika/segul/wiki/5.-Usages#finding-unique-ids-in-alignments).

#### Map sample distribution in a collection of alignments

If you would like to know how the samples distributed across your alignments, you only need to add the `--map` flag when searching for unique IDs. It will generate both the unique IDs (in a text file) and the sample distribution in boolean format (save to a csv file).

Using a `--dir` (or `-d` in short version):

```Bash
segul id --dir [a-path-to-a-directory] --input-format [sequence-format-keyword] --map
```

Using `--input` (or `-i` in short version) option ([learn the difference](https://github.com/hhandika/segul/wiki/4.-Command-Options#input-options)):

```Bash
segul id --input [a-path/wilcard-to-files] --map
```

Learn more about mapping sample distribution [here](https://github.com/hhandika/segul/wiki/5.-Usages#mapping-sample-distribution-in-a-collection-of-alignments).

#### Splitting concatenated alignments by partitions

To split alignment by partitions, you need the alignment file and the alignment partition. If the input partition file is not provided, `segul` will assume the partition is in the alignment file. It will fail to run if it could not find the partition in the input alignment.

The input option accepts a single file only. The `--dir` option is not available for this subcommand.

```Bash
segul split --input [a-path-to-a-concatenated-alignment] --input-part [a-path-to-partition-file]
```

By default, `segul` will use the alignment name as an output directory. To provide the output directory name, use the `--output` or `-o` option.

For the resulting alignments, `segul` will use the locus name in the partition file as the output filename. You can prefix this filename by using `--prefix` option. The resulting filenames will be {prefix}\_{locus-name}.{file-extension}.

Learn more about splitting concatenated alignments [here](https://github.com/hhandika/segul/wiki/5.-Usages#splitting-concatenated-alignments-by-partition).

#### Batch removing sequence based on user-defined IDs

Based on a list of IDs, you can remove sequences in a collection of alignments.

```Bash
segul remove --dir [alignment-dir] -f [sequence-format-keyword] --id [list-of-id]
```

Using regular expression:

```Bash
segul remove --dir [alignment-dir] -f [sequence-format-keyword] --re=["regex"]
```

If you remove more than a half of the sequence, `segul`'s extract feature is more efficient.

#### Batch renaming sequence IDs

To rename sequence IDs in multiple alignments, you need to input the sequence IDs in tsv or csv format with header. For example:

| Original_names        | New_names                |
| --------------------- | ------------------------ |
| Genus_species1_random | Genus_species1_voucherID |
| Genus_species2_random | Genus_species2_voucherID |

To simplify this process, you can generate unique IDs for all of your alignments using the `id` sub-command.

```Bash
segul id --dir [alignment-dir] -f [sequence-format-keyword]
```

Copy the IDs to a spreadsheet application, such as Microsoft Excel and then write new names and the header names as above. Save the file as csv or tsv. The program will infer the file format based on the file extension. Use it as an `--names` or `-n` input:

```Bash
segul rename --dir [alignment-dir] -f [sequence-format-keyword] --replace [file-path-to-IDs]
```

Using `--input` (or `-i` in short version) option ([learn the difference](https://github.com/hhandika/segul/wiki/4.-Command-Options#input-options)):

```Bash
segul rename --input [a-path/wilcard-to-files] --replace [file-path-to-IDs]
```

Example:

```Bash
segul rename -d uce-loci/ -f nexus --replace new_names.csv
```

You can also remove a part of the sequence name using a string input. Command below will remove each occurrence of `LOCUS` in the sequence IDs.

```Bash
segul rename -d -d uce-loci/ -f nexus --remove="LOCUS"
```

Using regular expression:

```Bash
segul rename -d -d uce-loci/ -f nexus --remove-re="(?i)^LOCUS-\d{3}"
```

The regex above match the first occurence of `locus-` (case insensitive) followed by 3 digits. For multiple occurrences in a single sequence name, use the option `--remove-re-all`. If your sequence names as below:

```Text
species_epithet_locus-100
species_epithet_locus-102
species_epithet_locus-103
```

They will change to:

```Text
species_epithet
species_epithet
species_epithet
```

You can also replace a part of the sequence name with other string:

```Bash
segul rename -d -d uce-loci/ -f nexus --replace-from="GENE" --replace-to="LOCUS
```

Using regular expression:

```Bash
segul rename -d -d uce-loci/ -f nexus --replace-from-re="(?i)GENE" --replace-to="LOCUS
```

You can also change the output format by using `--output-format` or `-F` option.

Learn more about batch renaming sequence IDs [here](https://github.com/hhandika/segul/wiki/5.-Usages#batch-renaming-sequence-ids).

#### Translating DNA sequences

`segul` translate DNA sequence based on [NCBI Genetic Code Tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#top). List of the supported tables is available [here](https://github.com/hhandika/segul/wiki/5.-Usages#translating-dna-sequences).

To translate dna alignment to amino acid using `--dir` (or `-d` in short version):

```Bash
segul translate --dir [path-to-alignment-files] -f [sequence-format-keyword]
```

Using `--input` (or `-i` in short version) option ([learn the difference](https://github.com/hhandika/segul/wiki/4.-Command-Options#input-options)):

```Bash
segul translate --input [a-path/wilcard-to-files]
```

By default, the app will use the NCBI standard code table (Table 1). To set the translation table, use the `--table` option. For example, to translate dna sequences using NCBI Table 2 (vertebrate MtDNA):

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

Learn more about translating DNA sequences to amino acid sequences [here](https://github.com/hhandika/segul/wiki/5.-Usages#translating-dna-sequences).

## Logging

Most information that is printed to the terminal is written to the log file (named `segul.log`). It is written to the current working directoy. Unlike the terminal output that we try to keep it clean and only show the most important information, the log file will also contain the dates, times, and the log level status. Each time you run the app, if the log file exists in the same directory, the app will append the log output to the same log file. Rename this file or move it to a different directory if you would like to keep a different log file for each task.

Learn more about using SEGUL [here](https://github.com/hhandika/segul/wiki/5.-Usages).

## Contribution

We welcome any kind of contribution, from issue reporting, ideas to improve the app, to code contribution. For ideas and issue reporting please post in [the Github issues page](https://github.com/hhandika/segul/issues). For code contribution, please fork the repository and send pull requests to this repo

<!-- ## Acknowledgment

We thank Giovani for testing the earlier version of the app and provide some ideas to further develop it. Some `SEGUL` features are inspired by Phyluce pipelines.  -->
