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

All of the formats are supported in interleave and sequential.

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
  - [Installing through Cargo](#installing-through-cargo)
  - [Installing from a GitHub repository](#installing-from-a-github-repository)
- [Usages](#usages)
  - [Available command options](#available-command-options)

## Supported Platforms

The program may work in any rust supported [platform](https://doc.rust-lang.org/nightly/rustc/platform-support.html). Below is a list of operating system that we tested and is guaranteed to work:

- Linux
- MacOS
- Windows
- Windows Sub-System for Linux

## Quick Start

If you already using a rust program and familiar with its [toolchain](https://www.rust-lang.org/learn/get-started), the best option is to install the app using cargo. If you are new to using a command line application, installing through cargo is also the easiest route ([see details in the installation instruction](#installation)). To install through cargo:

```{Bash}
cargo install segul
```

You can also use the pre-compiled binary available in [the release page](https://github.com/hhandika/segul/releases/). The installation is similar to any other single executable command line app, such as the phylogenetic programs IQ-Tree and RaXML. You only need to make sure the path to the app is registered in your environment variable, so that the app can be called from anywhere in your system. If you are not familiar with it, please see the detailed instruction below.

The program command structure is similar to git, gh-cli, or any other program that use subcommands. The program file name will be `segul` for Linux/MacOS and `segul.exe` for Windows.

```{Bash}
[THE-PROGRAM-FILENAME] <SUBCOMMAND> [OPTIONS] <VALUES> <A-FLAG-IF-APPLICABLE>
```

To check for available subcommand:

```{Bash}
segul --help
```

To check for available options and flags for each sub-command:

```{Bash}
segul <SUBCOMMAND> --help
```

For example, to concat all the alignments in a directory named `nexus-alignments`:

```{Bash}
segul concat --dir nexus-alignments --format nexus
```

It is also available in short options:

```{Bash}
segul concat -d nexus-alignments -f nexus
```

The program outputs are the resulting files from each task and a log file. Most information that is printed to the terminal is written to the log file. Unlike the terminal output that we try to keep it clean and only show the most important information, the log file will also contain the dates, times, and the log level status. The program will append the log outputs to the same log file (named `segul.log`) each time you run the program if the file exist in the same directory. Rename this file or move it to a different folder if you would like to keep a different log file for each task.

## Installation

We want the installation as flexible as possible. We offer three ways to install the program. Each of the options has its own advantages and disavantages.

1. Using a pre-compiled library. The quickest and the most straigtforward installation route, but the app may not fine-tuned for your specific hardware.
2. Through the rust package manager, cargo. It is similar to pip for python, gem for ruby, npm for JavaScript, etc. This is the recommended option. Cargo is a part of the Rust programming language toolchain. The toolchain is small and easy to install ([see below](#installing-through-cargo)). Cargo will also help to manage the app and allow for a quick update whenever the new version is released.
3. Compiling it from the source in the Github repository. It gives you the most up to date version, but the resulting binary may be less stable than the binary installed through the other routes.

### Using a pre-compiled library

The pre-compiled library is available at [the release page](https://github.com/hhandika/segul/releases/). The typicall workflow is as folow:

1. Copy the link to the zip file in [the release page](https://github.com/hhandika/segul/releases/) or download it to your computer.
2. Extract the zip file.
3. Make the binary executable (Linux and MacOS only).
4. Put it in a folder registered in your environment variable. It can be run from the folder where you extract the app too, but you will need your datasets in the same folder.

See specific details below:

#### Linux/MacOS

First, copy the link to the zip file in [the release page](https://github.com/hhandika/segul/releases/). We provide two versions of the app for Linux. The zip file labeled with HPC is compiled using Red Hat Enterprise Linux Server 7.9 (Kernel Version 3.10). If you are running the program in HPC, you should use this version. The other version (labeled Linux only) is compiled using [Ubuntu 20.04 LTS (Kernel version 5.8)](https://github.com/actions/virtual-environments/blob/main/images/linux/Ubuntu2004-README.md). You should use this if you are using more up to date Linux. Simply put, if you encounter [GLIBC](https://www.gnu.org/software/libc/) error, try using the HPC version. If the issue still persists, try install using cargo.

For MacOS, the executable is available for an Intel Mac. If you are using Apple silicon Macs (M1), may be best to install it through cargo.

Here, we use the version 0.3.1 as an example. You should replace the link with the most up to date version available in the release page.

- Download the library.

```{Bash}

wget https://github.com/hhandika/segul/releases/download/v0.3.1/segul-MacOS-x86_64.zip
```

- Unzip the file.

```{Bash}
unzip segul-MacOS-x86_64.zip
```

- Make it executable

```{Bash}
chmod +x segul
```

If you would like the binary executable for all users:

```{Bash}
chmod a+x segul
```

The next step is putting the binary in a folder registered in your path variable. It is always best to avoid registering too many paths in your environment variable. It will slow down your terminal startup if you do. If you already used a single executable cli app, the chance is that you may already have a folder registered in your path variable. Copy SEGUL executable to the folder. Then, try call SEGUL from anywhere in your system:

```{Bash}
segul --version
```

It should show the SEGUL version number.

If you would like to setup a folder in your environment variable, take a look at simple-qc [installation instruction](https://github.com/hhandika/simple-qc).

#### Windows

The installation procedure is similar to the MacOS or Linux. After downloading the zip file for Windows and extracting it, you will setup your environment variable that point to the path where you will put the executable. Follow this amazing guideline from the stakoverflow [to setup the environment variable](https://stackoverflow.com/questions/44272416/how-to-add-a-folder-to-path-environment-variable-in-windows-10-with-screensho). After setup, copy the segul.exe file to the folder.

### Installing through Cargo

This is the recommended option. Cargo will compile the program and fine-tuned it for your specific hardware. It also allows to easily updating the app. You will need the [rust compiler toolchain](https://www.rust-lang.org/learn/get-started). It requires rust version 1.5 or higher. Then, check if the toolchain installation successful:

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

For OpenSUSE:

```{Bash}
zypper install -t pattern devel_basis
```

For Fedora:

```{Bash}
sudo dnf groupinstall "Development Tools" "Development Libraries"
```

For Windows, you only need to install the GNU toolchain for rust. The installation should be straighforward using rustup. Rustup comes as a part of the rust-compiler toolchain. It should be available in your system at the same time as you install cargo.

```{Bash}
rustup toolchain install stable-x86_64-pc-windows-gnu
```

Then set the GNU toolchain as the default compiler

```{Bash}
rustup default stable-x86_64-pc-windows-gnu
```

### Installing from a GitHub Repository

You will need [rust compiler toolchain](https://www.rust-lang.org/learn/get-started). To install the development version for any supported platform:

```{Bash}
cargo install --git https://github.com/hhandika/segul.git
```

You should have segul ready to use.

It is equivalent to:

```{Bash}
git clone https://github.com/hhandika/segul

cd segul/

cargo build --release
```

The different is, for the latter, the executable will be in the `segul` repository: /target/release/segul. Copy the `segul` binary and then add it to your environment path folder.

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
