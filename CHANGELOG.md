# Release Notes

## Upcoming

- Add spinner for converting a single file.
- Output names for concatenating filtered sequences is now default to the directory name with parameters as sufixes.
- Allow "ignore" in the summary data types.

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

- First stable release.
