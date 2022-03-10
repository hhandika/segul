# Release Notes

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

## vO.13.1

- Fix issues including empty sequences in splitted alignments.
- Fix unmatched nchar counts in nexus and phylip output when splitting alignments.

## v0.13.0

- Fix errors in the nexus output datatype for amino acid sequence.
- Add a feature to split alignment given partition information.

## v0.12.0

- Add taxon summary to the summary stats feature.
- Summary stats now write the csv ouput files to a directory.
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

- Add a DNA to Amino Acic translation feature.

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
- Output names for concatenating filtered sequences is now default to the directory name with parameters as sufixes.
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
