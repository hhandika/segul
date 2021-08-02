# SEGUL

![Segul-Tests](https://github.com/hhandika/segul/workflows/Segul-Tests/badge.svg)
![Crate-IO](https://img.shields.io/crates/v/segul)
![GH-Release](https://img.shields.io/github/v/release/hhandika/segul)
![License](https://img.shields.io/github/license/hhandika/segul)

SEGUL is an ultrafast and memory efficient command-line (cli) application for working with sequence alignments that typically done using interpreted programming languages, such as Python, R, or Perl. It is designed to handle genomic datasets, but just as capable for Sanger datasets. In our test using a dataset with 4060 UCE loci, for instance, compare to a program written using biopython library, SEGUL is >40x faster for alignment concatenation while using 3x less RAM space.

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

## Table of Contents

- [Supported Platforms](#supported-platforms)
- [Quick Start](#quick-start)
- [Installation](#installation)
  - [Using a pre-compiled library](#using-a-pre-compiled-library)
  - [Installing using the Cargo Package Manager](#installing-using-the-cargo-package-manager)
  - [Installing from a GitHub repository](#installing-from-a-github-repository)
- [Command Structure](#command-structure)
  - [Available subcommands](#available-subcommands)
- [Usages](#usages)
  - [Converting sequences to a different format](#converting-sequences-to-a-different-format)
  - [Concatenating sequences](#concatenating-sequences)

## Supported Platforms

The program may work in any rust supported [platform](https://doc.rust-lang.org/nightly/rustc/platform-support.html). Below is a list of operating system that we tested and is guaranteed to work:

- Linux
- MacOS
- Windows
- Windows Subsystem for Linux (WSL)

## Quick Start

If you already using a rust program and familiar with its [toolchain](https://www.rust-lang.org/learn/get-started), the best option is to install the app using cargo. If you are new to using a command line application, installing through cargo is also the easiest route ([see details in the installation instruction](#installing-through-cargo)). To install through cargo:

```Bash
cargo install segul
```

You can also use the pre-compiled binary available in [the release page](https://github.com/hhandika/segul/releases/). The installation is similar to any other single executable command line app, such as the phylogenetic programs IQ-Tree and RaXML. You only need to make sure the path to the app is registered in your environment variable, so that the app can be called from anywhere in your system ([see instructions](#using-a-pre-compiled-library)).

The program command structure is similar to git, gh-cli, or any other program that use subcommands. The program file name will be `segul` for Linux/MacOS/WSL and `segul.exe` for Windows.

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

> :warning: **Unlike the log file, for the other outputs, the app will over-write existing files with the same names**: Careful in specifying the output file names. Future updates will prevent it.

## Installation

We want the installation as flexible as possible. We offer three ways to install the program. Each of the options has pros and cons.

1. Using a pre-compiled library. The quickest and the most straigtforward installation route, but the app may not be fine-tuned for your specific hardware.
2. Using the rust package manager, cargo. It is similar to pip for python, gem for ruby, npm for JavaScript, etc. This is the recommended option. Cargo is a part of the Rust programming language toolchain. The toolchain is small and easy to install ([see below](#installing-through-cargo)). Cargo will also help to manage the app and allow for a quick update whenever the new version is released.
3. Compiling it from the source in the Github repository. It gives you the most up to date version, but the resulting binary may be less stable than the binary installed through the other routes.

### Using a pre-compiled library

The pre-compiled library is available in [the release page](https://github.com/hhandika/segul/releases/). The typical workflow is as follow:

1. Download the zip file in [the release page](https://github.com/hhandika/segul/releases/) to your computer, using a browser or using a command line app.
2. Extract the zip file.
3. Make the binary executable (Linux and MacOS only).
4. Put it in a folder registered in your environment variable. It can be run from the folder where you extract the app too, but this option is less convenient when dealing with datasets in different folders or storage partitions.

See specific details below:

#### Linux/WSL/MacOS

First, copy the link to the zip file in [the release page](https://github.com/hhandika/segul/releases/). We provide two versions of the app for Linux. The zip file labeled with HPC is compiled using Red Hat Enterprise Linux Server 7.9 (Kernel Version 3.10). If you are running the program in HPC, you should use this version. The other version (labeled Linux only) is compiled using [Ubuntu 20.04 LTS (Kernel version 5.8)](https://github.com/actions/virtual-environments/blob/main/images/linux/Ubuntu2004-README.md). You should use this if you are using WSL or more up to date native Linux distros. Simply put, if you encounter [GLIBC](https://www.gnu.org/software/libc/) error, try using the HPC version. If the issue still persists, try to [install the app using cargo](#installing-through-cargo).

For MacOS, the executable is available for an Intel Mac. If you are using Apple silicon Macs (Apple M1), we recommend installing it using cargo.

Here, we use the version 0.3.1 as an example. You should replace the link with the most up to date version available in the release page.

- Download the library.

```Bash

wget https://github.com/hhandika/segul/releases/download/v0.3.1/segul-MacOS-x86_64.zip
```

- Unzip the file.

```Bash
unzip segul-MacOS-x86_64.zip
```

- Make it executable

```Bash
chmod +x segul
```

If you would like the binary executable for all users:

```Bash
chmod a+x segul
```

The next step is putting the binary in a folder registered in your path variable. It is always best to avoid registering too many paths in your environment variable. It will slow down your terminal startup if you do. If you already used a single executable cli app, the chance is that you may already have a folder registered in your path variable. Copy SEGUL executable to the folder. Then, try call SEGUL from anywhere in your system:

```Bash
segul --version
```

It should show the SEGUL version number.

If you would like to setup a folder in your environment variable, take a look at simple-qc [installation instruction](https://github.com/hhandika/simple-qc).

#### Windows

The installation procedure is similar to the MacOS or Linux. After downloading the zip file for Windows and extracting it, you will setup your environment variable that point to the path where you will put the executable. In Windows, this is usually done using GUI. Follow this amazing guideline from the stakoverflow [to setup the environment variable](https://stackoverflow.com/questions/44272416/how-to-add-a-folder-to-path-environment-variable-in-windows-10-with-screensho). After setup, copy the segul.exe file to the folder.

### Installing using the Cargo Package Manager

This is the recommended option. Cargo will compile the app, manage its dependencies, and fine-tuned it for your specific hardware. It also allows to easily updating the app.

First, download and install [the rust compiler toolchain](https://www.rust-lang.org/learn/get-started). It requires rust version 1.5 or higher. Then, check if the toolchain installation successful:

```Bash
cargo --version
```

It should show the cargo version number. Then, install the app:

```Bash
cargo install segul
```

If you encounter a compiling  issue (usually happens on Linux or Windows), you may need to install the C-development toolkit. For Debian-based Linux distribution, such as Debian, Ubuntu, PopOS, etc., the easiest way is to install build-essential:

```Bash
sudo apt install build-essential
```

For OpenSUSE:

```Bash
zypper install -t pattern devel_basis
```

For Fedora:

```Bash
sudo dnf groupinstall "Development Tools" "Development Libraries"
```

For Windows, you only need to install the GNU toolchain for rust. The installation should be straighforward using rustup. Rustup comes as a part of the rust-compiler toolchain. It should be available in your system at the same time as you install cargo.

```Bash
rustup toolchain install stable-x86_64-pc-windows-gnu
```

Then set the GNU toolchain as the default compiler

```Bash
rustup default stable-x86_64-pc-windows-gnu
```

Try to install SEGUL again:

```Bash
cargo install segul
```

### Installing from a GitHub Repository

You will need [rust compiler toolchain](https://www.rust-lang.org/learn/get-started). The setup procedure is similar to installing the app using cargo. To install the development version for any supported platform:

```Bash
cargo install --git https://github.com/hhandika/segul.git
```

You should have SEGUL ready to use.

It is equivalent to:

```Bash
git clone https://github.com/hhandika/segul

cd segul/

cargo build --release
```

The different is that, for the latter, the executable will be in the `segul` repository: `/target/release/segul`. Copy the `segul` binary and then add it to your environment path folder.

Then, try to call SEGUL:

```Bash
segul --version
```

## Updating the app

If you install the app using cargo, updating the app is the same as installing it:

```Bash
cargo install segul
```

Cargo will check whether the version of the app in your computer different from the version in the rust package repository ([crates.io](https://crates.io/crates/segul)) and will install the newer version if it is available. Similar procedure is also applied for installing from the GitHub repository:

```Bash
cargo install --git https://github.com/hhandika/segul.git
```

If you used the pre-compiled binary, replace the old binary with the newer version manually.

## Uninstalling the app

It is also easy to do if you install the app using cargo:

```Bash
cargo uninstall segul
```

Rust toolchain, including cargo, can be uninstall easily too:

```Bash
rustup self uninstall
```

Remove the app manually if you use the pre-compiled binary.

## Command Structure

### Available subcommands

```Bash
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

### Input options

- Option `-i` or `input`: Use for a single file input. Only available for convert and summary subcommands.
- Option `-d` or `dir`: If your input is a path to a directory. The directory input requires users to specify the input format. Available for all subcommands.
- Option `w` or `wildcard`: If your input is wilcards. This is more flexible than the other two input options and can accept multiple values. Available for all subcommands.

### Input format

Arguments: `-f` or `--input-format`

Availabilities: all subcommands

It is used to specify the input format. For a sinlge input `-i` or `--input` and `-w` or `--wildcard`, this is not required.

Input format options (all in lowercase):

- `auto` (default)
- `nexus`
- `phylip`
- `fasta`

### Output

Arguments: `-o` or `--output`

Availabilities: all subcommands

For a single output task, such as converting a single file, or concatenating alignment, the output will be the file name for the output. For a multiple output task, such as converting multiple files to a different format, the output will be the directory name for the output. The program will use the input file name for the each output file.

The program by default write to the current working directory.

### Output format

Arguments: `-F` or `--output-format`

Availabilities: all subcommands

By default the output format is `nexus`. Use this option to specify the output format. Below is the available output formats.

Sequential formats:

- `nexus`
- `phylip`
- `fasta`

Interleaved formats:

- `fasta-int`
- `nexus-int`
- `phylip-int`

### Data type

Argument: `--datatype`

Availabilities: all subcommands

The app support both DNA and amino acid sequences. By default the datatype is set for DNA. If your input file is amino acid sequences, you will need to change the data type to `aa`. By specifying the data type, the app will check if your sequence files contain only IUPAC characters. Except for computing summary statistics, you can set data type to `ignore` to skip checking the IUPAC characters. This usually speed app the computation for about 40%. Use this option when you are sure your sequences contain only IUPAC characters.

### Special options

#### Partition format

Arguments: `-p` or `--part`

Availabilities: concat and filter subcommands

This option is used to specify the partition format. By default the format is nexus. Available options:

- `charset` (embedded in a nexus sequence)
- `nexus`
- `raxml`

#### Percentage interval

Arguments: `-i` or `--interval`

Availability: summary subcommand

This option is to specify the percentage decrement interval for computing data matrix completeness in summary statistics. Available interval: `1`, `2`, `5`, `10`.

### Filtering options

Only available for filter subcommands. Available options:

- `-l` or `--len`: To filter alignments based on minimal alignment length
- `--percent`: To filter based on percentage of data matrix completeness.
- `--npercent`: The same as `--percent`, but accept multiple values.
- `--pinf`: To filter based on the number of parsimony informative sites.
- `--ntax`: To defined the total number of taxa. By default the app determines the number of taxa in all the alignments based on the numbers of unique IDs.

### Special Flags

`--codon`: Use to set the partition format to codon model. Available in concat subcommand, and filtering subcommand (if you choose to concatenate the result).

`--concat`: Available for filtering to concat the filtering results in lieu to copying the files.

`--sort`: Available for convert to sort the sequences based on their IDs in alphabetical order.

`-h` or `--help`: To display help information.

`-v` or `--version`: To display the app version information.

## Usages

### Converting sequences to a different format

Segul can convert from a single file input (`-i` or `--input`), a directory input (`-d` or `--dir`), or using a wildcard input (`-c` `--wildcard`).

#### Using a single file input

```Bash
segul convert --input [path-to-your-repository] --input-format [choose-one]
```

In short format:

```Bash
segul convert -i [path-to-your-repository] -f [sequence-format]
```

By default it converts to nexus. To choose a different output format, use the output format option (`-F` or `--output-format`):

```Bash
segul convert --input [path-to-your-repository] --input-format [sequence-format] --output-format [sequence-format]
```

In short format, notice the uppercase 'F' for the output format:

```Bash
segul convert -i [path-to-your-repository] -f [sequence-format] -F [sequence-format]
```

You can also skip specifying the input format and the program will infer it based on the file extension:

```Bash
segul convert -i [path-to-your-repository]
```

By default the program will use the input file name for the output. To specify, the output name use the `-o` or `--output` option. There is no need to include the extension for the output name.

Using the `--sort` flag, you can also sort the sequence based on their IDs in alphabetical order.

For example, to convert a file name `sequence.fasta` to a phylip format and we will sort the result.

```Bash
segul convert -i sequence.phy -f fasta -F phylip -o new_sequence --sort
```

#### Batch converting sequence files in a directory

The conversion command for a directory input is similar to converting a single file. Unlike the single file, the app require you to specify the input format and the output name. The output name will be the directory name for the output files, whereas the output file name will be the same as the input file.

```Bash
segul convert --dir [path-to-your-repository] --input-format [choose-one] --output [your-output-dir-name]
```

In shortened format

```Bash
segul convert -d [path-to-your-repository] -f [sequence-format] -o [your-output-dir-name]
```

For example, suppose we want to convert all the fasta files in the directory below to a phylip format and name the output directory `new_sequences`:

```Bash
sequences/
├── seq_1.fasta
├── seq_2.fasta
└── seq_3.fasta
```

The command will be:

```Bash
segul convert -d sequences/ -f fasta -F phylip -o new_sequences
```

The resulting directory will be:

```Bash
new_sequences/
├── seq_1.phy
├── seq_2.phy
└── seq_3.phy
```

#### Batch converting sequence files using wildcard

All the options for a single input or a directory is also available for a wildcard. The program can also infer the input format. Unlike any other input, the wildcard can take multiple values. This allow you to batch converting files in different folders. The ouput will be in a single directory. It is required to specify the output name and will be used as a name for the output directory.

```Bash
segul convert -c [wildcard-1] [wildcard-2] [wildcard-3] -f [sequence-format] -o [your-output-dir-name]
```

### Concatenating sequences

The app concat multiple alignments and write the partition setting for the resulting files. The input options are `-d` or `--dir` and `-c` or `--wildcard`. To specify the partition format, you will use the `-p` or `--part` option. You can also write the partition to a codon model format by using the flag `--codon`.

### Filtering sequences

### Computing sequence summary statistics

### Finding all IDs in a collection of alignments
