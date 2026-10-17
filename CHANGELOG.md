# Release Notes

## v0.24.0

Library Changes:

- The library code now lives in a new crate, [`segul-core`](https://crates.io/crates/segul-core). The `segul` crate is the command-line app.
- `segul-core` uses Cargo features to enable only what you need: `align`, `sequence`, `read`, `contig`, and `maf`. All features are on by default. For example:

  ```toml
  segul-core = { version = "0.24", default-features = false, features = ["read"] }
  ```

- Library code no longer depends on `clap` or `dialoguer`. `noodles` is only needed by the `read` feature.

Migration:

- Replace `segul::` with `segul_core::` in your imports.
- `segul` still re-exports `core`, `helper`, `parser`, `stats`, and `writer` from `segul-core`, so existing code keeps compiling. These re-exports will be removed in a future release.

Breaking Changes:

- MAF conversion is now `segul maf convert` (was `segul maf maf`), matching the documentation.
- `--bed <FILE>` replaces `--reference <FILE> --from-bed`. The BED file is optional.
- The MAF output format option is now `-F/--output-format`, like other commands (was `-t`). NEXUS output is supported. The default output directory is `Maf-Convert`.
- Sample names drop the assembly suffix of MAF sources (`hg18.chr7` becomes `hg18`), so each species is one sample across blocks. Use `--full-source-name` to keep full names.
- Blocks that match no BED region are no longer written by default. Use `--write-unmatched` to write them to `unmatched/` (was `missing-refs/`).
- `-d` finds `.maf` and `.maf.gz` files directly in the directory, without searching subdirectories, like other commands.
- Library: `MafReader` yields `anyhow::Result<MafParagraph>`, `MafConverter` uses a builder and returns a `MafConvertReport`, and `BedParser::parse` takes `&self` and returns `anyhow::Result`.

New Features:

- `segul maf summary`: per-file and per-species MAF statistics (blocks, alignment length, reference coverage, gaps) written to CSV.
- `segul maf convert` with a BED file assigns each block to every locus it overlaps on the same chromosome, clips it to the locus (`--no-clip` keeps whole blocks), and merges blocks of a locus across all input files, padding missing species with gaps. A `maf_convert.csv` report lists every locus.
- `segul maf convert` without a BED file writes one file per block, or one per reference chromosome with `--group-by chrom`.
- Block filters: `--min-block-len`, `--min-species`, `--min-score`, and `--species` or `--species-file` to keep only listed species.
- `--strand-aware` reverse complements loci on the minus strand of the BED file.
- `--ref-species` sets the reference species. It defaults to the first species of the first block.
- Gzipped MAF (`.maf.gz`) and BED (`.bed.gz`) files are read directly.
- `segul sequence transpose`: transpose locus alignments into per-sample sequence files. Each sample gets its own file with one sequence per locus, labeled by locus name, with gaps and missing data removed.
- `segul align distance`: compute pairwise genetic distances for each alignment. Supports uncorrected p-distance, Jukes-Cantor (JC69), Kimura 2-parameter (K2P), Tamura-Nei (TN93), and Poisson (amino acid) models, using pairwise deletion for gaps, missing data, and ambiguous characters. Writes a square distance matrix for each locus to `per_locus/` (use `--triangle lower` or `--triangle upper` to leave the other triangle empty, and `--diagonal zero`, `blank`, or `max` to set the diagonal) and a summary of each locus's distances (min, max, mean, median, SD) to the output directory.
- `segul align distance` group means: group samples into populations or species with `--group-delimiter` and `--group-fields` (e.g. `_` and 2 turn `Rattus_rattus_FMNH1` into `Rattus_rattus`), `--group-regex`, `--group-file` (a CSV/TSV with sample and group columns), or `--group-auto` to detect the delimiter and field count from sample IDs. Writes within- and between-group mean distances for each locus to `per_locus_group/`, and pooled over all sample pairs of all loci to `<prefix>_group_mean.csv`, with pair counts in `<prefix>_group_pairs.csv`. Use `--group-diagonal` to write the within-group mean (default), `zero`, `blank`, or `max` (each group's highest mean to another group) on the diagonal of the group mean matrices.
- `segul align distance --long`: write the pairwise distances of all loci to one long table, `<prefix>_pairwise_long.csv` (locus, samples, groups, compared sites, distance), and with grouping, group statistics (pairs, mean, SD, min, max) per locus and pooled to `<prefix>_group_long.csv`.

Bug Fixes:

- Fix MAF conversion without a BED file failing when a reference has more than one block.
- Fix BED files with four columns failing to parse because of the line ending. BED lines with CRLF endings, blank lines, or too few columns are handled or reported with the line number.
- Fix MAF blocks of the same locus overwriting each other instead of being joined.
- Fix BED matching that used only the start position and ignored the chromosome.
- Fix consecutive MAF blocks merging when a file with CRLF line endings is read on Linux or macOS.
- Fix negative MAF block scores being ignored.

## v0.23.2

- Add Windows ARM64 support.

## v0.23.1

- Fix ARM support for Linux.

## v0.23.0

New Features:

- [Convert alignments to unaligned sequences](https://www.segul.app/docs/cli-usage/align-unalign).
- [Add sequence to an existing sequence files/alignments](https://www.segul.app/docs/cli-usage/sequence-add).
- [Trim alignments](https://www.segul.app/docs/cli-usage/align-trim).
- New options for filtering alignments ([see details](https://www.segul.app/docs/cli-usage/align-filter)):
  - Minimum or maximum sequence length.
  - Minimum or maximum of parsimony informative sites
  - Minimum taxon counts
  - Based on user-defined list of sequence IDs.

Breaking Changes:

- Filtering argument changes for better consistency (see [docs](https://www.segul.app/docs/cli-usage/align-filter)).
- Remove `--ntax` option from [alignment filtering](https://www.segul.app/docs/cli-usage/align-filter) because by default SEGUL automatically and fastly counts the number of unique taxa across all input alignments.

Bug Fixes:

- Fix max-gap filtering issues.
- Fix concatenation leading to missing data when sequence IDs in FASTA files contain trailing whitespace.

Other Changes:

- Migrate to Rust 2024 edition.
- Update dependencies.

## v0.22.1

- Sequence filtering now does not enforce input is an alignment.

## v0.22.0

New Features:

- Filter sequences based on gaps or sequence length.
- Add a proportion of missing data parameter to the alignment filtering.

API Breaking Changes:

- Rename classes and methods to improve clarity and consistency.
- Remove unnecessary mutability for conversion modules.
- Partition conversion moved to align module.
- Output directory and output format arguments for extraction are now required in the new instance of the extraction struct.

## v0.21.3

- Try to fix ARM CI issues.

## v0.21.2

- Fix auto detect sequence format issue does not match the GUI input.

## v0.21.1

- Fix extra commas in the alignment summary CSV output.
- Fix missing proportion of missing data in alignment summary.
- Fix issues when inputting list of IDs for sequence extraction and removal. Now, require equal sign and semi-colon to separate the IDs.
- Fix spacing and incorrect symbol issues in the terminal output.

### v0.21.0

- Alignment summary now is written to file.
- Read summary with `complete` option now writes per read position results to compressed zip file.
- Prefix option for read, contig summary, and ID extraction.
- More consistent CSV header names for all summary statistics, and sequence ID map.

## v0.20.2

- Fix CI issues.

## v0.20.1

- Fix empty result when summarizing a single contig file.
- Fix CI issues.

## v0.20.0

- ID parsing result is now in a directory consistent with other subcommands.

## v0.19.2

- Fix CI issues.

## v0.19.1

- Fix CI issue when building windows binary.

## v0.19.0

New features:

- Directory input auto detect sequence format based on the file extension.
- New summary statistics for raw reads and contigs.
- New command line structure.

Bug fixes:

- Fix translation table errors.
- Update deprecated dependencies.
- Fix extract output issues.

## v0.18.1

- Fix colored terminal dependency deprecation issues.
- Fix inconsistent header alignment when printing summary statistics to the terminal.

## v0.18.0

New features:

- Add an option to generate sequence statistics for each locus.
- More taxon statistics: proportion of missing data, nucleotide counts, GC content and AT content.
- Add message when choosing "no" folder exists.

Bug Fixes:

- Fix compiling issues.
- Fix command conflicts when using the option `remove-re-all`.

## v0.17.0

New Features:

- Add a feature to remove sequence based on user-defined IDs.
- Add regular expression support to rename sequence IDs.
- Add string replacement option to rename sequence IDs.
- Add option to filter based on the present of all user-defined taxa.
- Replace site counts with character counts in the taxon summary.
- Add integration tests.

Bug Fixes:

- Fix error messages when finding empty matrix.
- Fix missing character counts in amino acid alignment summary.
- Fix input and output format character cases.
- Fix conflict with percent when filtering alignments.
- Fix command conflict errors.
- Make output name consistent with subcommand name for statistics function.

## v0.16.3

- Fix output option show up when it does not use by the partition subcommand.
- Fix unclear error messages when failing to infer partition based on the file extension.

## v0.16.2

- Fix version number errors.

## v0.16.1

- Fix conflicting short commands.

## v0.16.0

- New feature to convert partition formats.
- Now, can parse codon partition.
- Check partition format in the sequence file when splitting alignments.
- Fix incorrect data type in raxml partition when the input is amino acid sequences.
- Clearer command help output (now use [clap](https://docs.rs/clap/latest/clap/) v3).

## v0.15.0

- Add an overwriting flag to all subcommands to allow overwriting existing output files or directories.
- Fix typos in error messages.
- Check sequences for duplicate IDs are the same. Let users know in the error messages.

## v0.14.0

- Add support for RaXML partition that does not contain datatype.
- Add an option to prefix output filename when splitting alignments.
- Replace dot with underscore in the partition locus name.
- Print partition information when splitting alignments.

## v0.13.1

- Fix issues including empty sequences in split alignments.
- Fix unmatched number of character counts in nexus and phylip output when splitting alignments.

## v0.13.0

- Fix errors in the nexus output datatype for amino acid sequence.
- Add a feature to split alignment given partition information.

## v0.12.0

- Add taxon summary to the summary stats feature.
- Summary stats now write the csv output files to a directory.
- Allow to add prefix for summary stats csv output files.
- Fix spaces issue in showing filtering parameters in the log output.

## v0.11.2

- Fix confusing input name in the terminal output.
- Improve fasta parsing performance.
- Fix issues when the app keeps processing when no alignments left after filtering.

## v0.11.1

- Fix file count formatting for sequence extraction terminal output.
- Fix outdated error messages.

## v0.11.0

- Add feature to map sample distribution across a collection of alignments.
- Dir and file exist prompts now show the file names.
- Simplify input option. Wildcard option is now a part of input option.
- Print output format.

## v0.10.1

- Fix spinner messages.
- Fix dry-run terminal output.

## v0.10.0

- Add feature to batch renaming sequence IDs.

## v0.9.1

- Fix missing file extension issues when converting multiple sequences

## v0.9.0

- Add a DNA to Amino Acid translation feature.

## v0.8.6

- Fix inconsistent DNA character ordering when using a single input for summary statistics.

## v0.8.5

- Fix missing wildcard options for filter and id subcommands.

## v0.8.4

- Fix errors in displaying proportion of parsimony informative sites for a single file input.

## v0.8.3

- Minor performance improvement for sequence extraction.
- Fix multiple newlines after spinners when using ID subcommands.

## v0.8.2

- Fix wrong help messages in the extract command arguments.
- Shows output directory in the stdout after extracting sequences.

## v0.8.1

- Fix issues missing summary subcommands.

## v0.8.0

- Add new feature to extract sequence file based on user-defined IDs.

## v0.7.1

- Fix version unmatched version numbers.

## v0.7.0

- Add spinner for converting a single file.
- Output names for concatenating filtered sequences is now default to the directory name with parameters as suffixes.
- Allow "ignore" in the summary data types.
- Fix issues on specifying a directory for summary stats.
- Help info now shows the type of values required for the arguments.
- Fix typos in panic messages.

## v0.6.6

- Try to fix Windows CI issues.

## v0.6.5

- Fix CI issues

## v0.6.4

- Update release rules.

## v0.6.3

- CI release changes.

## v0.6.2

- Clearer terminal messages.

## v0.6.1

- Fix CI token access.

## v0.6.0

- Add support for percent parsimony informative.

## v0.5.0

- Allow to specify both the directory and file name for concatenating alignments.
- Clearer prompt messages.
- More consistent letter cases and spacing for spinner.
- Concat filtered alignment now default to nexus. No need to specify.

## v0.4.0

- Avoid overwriting existing files or directory. Now, ask user to remove or abort.
- Concatenate results now in a directory.
- Better error handling
- Fix stack overflow when getting unique ids.

## v0.3.16

- Change panic implementation to C-style abort.
- Fix issues with too many redundant error messages when panic occurs in multi-threading tasks.

## v0.3.15

- Print data type in the input information.

## v0.3.14

- Fix another CI issue.

## v0.3.13

- Update release steps.

## v0.3.12

- CI fix with newer container.

## v0.3.11

- Add support for older linux.

## v0.3.10

- Change CI implementation.

## v0.3.9

- Fix CI errors.

## v0.3.8

- More consistent error messages for the filter subcommand.

## v0.3.7

- Better error handling.

## v0.3.6

- Fix issues when concatenating filtered alignments. Now, enforce users to input required arguments.

## v0.3.5

- Change option "format" to "input-format" for consistency with "output-format".
- Fix missing a new line after spinner message in converting sequences in a directory.

## v0.3.4

- Fix wildcard issues on windows.

## v0.3.3

- Fix issue with summary printing.

## v0.3.2

- Fix error prematurely writing a log file.

## v0.3.1

- Fix welcome screen showing up when checking the app version.

## v0.3.0

- Add new features: alignment filtering, id finding, and summary statistics.
- More consistent terminal output.
- Add sequence checking for fasta input.
- Faster parser.
- Improvement in memory allocation.

## v0.2.0

- Add support for interleaved outputs.

## v0.1.1

- Faster fasta concat.

## v0.1.0

- First release.
