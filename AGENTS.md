# Agent instructions

SEGUL is an ultrafast and memory-efficient tool for phylogenomics, written in
Rust. This file tells coding agents how the repository is organized and how
to work in it.

## Hard rules

- Never commit. Do not run `git commit`, `git push`, or anything that creates
  commits or tags. Leave all changes uncommitted for the maintainer to review
  and commit.
- Do not publish crates (`cargo publish`), create GitHub releases, or push
  tags. Dry runs (`cargo publish --dry-run`) are fine.
- Do not delete or overwrite files you did not create in this task, including
  benchmark results and datasets outside the repo.

## Repository layout

The repo is a Cargo workspace (edition 2024, version set once in the root
`Cargo.toml` under `[workspace.package]`).

- `crates/segul-core/`: the library. All analysis code lives here.
  - `src/core/`: one module per command group (`align`, `sequence`, `read`,
    `contig`, `maf`), one file per task (`align/split.rs`,
    `sequence/rename.rs`, ...).
  - `src/parser/`: file format readers (FASTA, NEXUS, PHYLIP, partitions,
    MAF, BED, FASTQ quality scores, delimited text).
  - `src/writer/`: sequence, partition, summary, and archive writers.
  - `src/stats/`: summary statistics.
  - `src/helper/`: shared types (`types.rs`: `InputFmt`, `OutputFmt`,
    `DataType`, `Header`, `SeqMatrix`, ...), file finders, logging, and
    utilities.
  - `tests/files/`: small test inputs used by unit tests, doctests, and the
    CLI integration tests.
- `src/`: the `segul` command-line app.
  - `src/cli/args/`: clap argument definitions (`align.rs`, `sequence.rs`,
    `genomics.rs`), shared argument groups in `args/mod.rs`.
  - `src/cli/<group>/<task>.rs`: a parser struct per command that turns args
    into a `segul-core` call.
  - `src/cli/commands.rs`: dispatches subcommands.
  - `src/lib.rs`: re-exports `core`, `helper`, `parser`, `stats`, `writer`
    from `segul-core` for backward compatibility. Do not add new code here.
- `tests/`: CLI integration tests, one file per command, sharing macros in
  `tests/utils.rs`.
- `bench/`: Python benchmark harness (`segul-bench`), managed with `uv`. See
  `bench/README.md`.
- `scripts/`: release helper scripts (`gh-newrelease`, `gh-upload`, ...).
- `.github/workflows/`: `test.yml` (check, test, feature matrix, install
  smoke test on Linux, macOS, Windows, x86_64 and arm64) and `release.yml`
  (build and upload binaries on tag push).
- `CHANGELOG.md`: release notes, newest first.

Downstream projects build against the `dev` branch: `pysegul`
(`segul-core`) and the Flutter GUI `segui` (`segul`). Treat public
`segul-core` items as API and call out breaking changes.

## Where new code goes

- Analysis logic goes in `crates/segul-core/src/core/<group>/<task>.rs`, with
  parsing, writing, and stats in the matching `parser/`, `writer/`, `stats/`
  modules. Library code must not depend on `clap` or `dialoguer`.
- Gate modules with the right Cargo feature (`#[cfg(feature = "align")]`,
  etc.) in the parent `mod.rs`. Shared code used by several features uses
  `#[cfg(any(feature = "...", ...))]`. Each feature must build on its own:
  `cargo check -p segul-core --all-targets --no-default-features --features <name>`.
- Add the CLI wiring in `src/cli/args/<group>.rs` (args), a
  `src/cli/<group>/<task>.rs` parser, and the match arm in
  `src/cli/commands.rs`.
- New dependencies go in `[workspace.dependencies]` in the root
  `Cargo.toml` and are referenced with `dep.workspace = true`. Optional ones
  belong behind the feature that needs them (see `read` and `noodles`).

## Coding style

Match the code around you. In particular:

- Format with `cargo fmt`, and keep `cargo clippy --workspace --all-targets`
  free of warnings.
- Import order: `std`, then external crates, then `crate::` paths, with a
  blank line between groups. Items from the same module are grouped in
  braces (`use std::path::{Path, PathBuf};`).
- Tasks are structs holding borrowed inputs (`&'a Path`, `&'a DataType`,
  `&'a InputFmt`, ...) with a `new` constructor, plus a method that runs the
  task (`split`, `concat_alignment`, ...). Long argument lists use a
  constructor, not a builder, unless the module already uses a builder (MAF
  conversion does).
- Shared behavior comes from small traits with default methods, implemented
  with empty impls at the top of the file (`impl OutputPrint for X<'_> {}`).
- Visibility is as narrow as possible: `pub(crate)` or `pub(in crate::cli)`
  in the CLI, `pub` only for library API.
- Parallelism uses `rayon` (`par_iter`). Counters across threads use
  atomics.
- User-facing progress uses `utils::set_spinner()` and `log::info!` with the
  aligned `"{:18}: {}"` layout and `colored` section headers.
- Errors: newer library code returns `anyhow::Result` (MAF, BED, writers).
  Existing code often uses `.expect("Failed ...")` with a short message
  starting with "Failed". Prefer `anyhow::Result` for new library functions
  that can fail on user input; keep `expect` messages descriptive.
- Doc comments: `//!` module headers describing what the module does, `///`
  on public items. Comments explain why, briefly. Keep comment density similar
  to neighboring code.
- Names are descriptive and spelled out (`input_fmt`, `output_fmt`,
  `datatype`, `part_path`). Formats use the existing enums, not strings.

## Tests

- Unit tests live in a `#[cfg(test)] mod test` at the bottom of each file and
  read inputs from `crates/segul-core/tests/files/`. Add a small fixture file
  there rather than generating data inline when a test needs a real format.
- CLI tests in `tests/<command>.rs` run the binary with `assert_cmd` in a temp
  dir, using the macros in `tests/utils.rs` (`initiate_cmd!`,
  `test_results!`). Paths are relative to the repo root.
- Public examples in doc comments are doctests; keep them compiling.
- Run before finishing:

  ```bash
  cargo fmt --all
  cargo clippy --workspace --all-targets
  cargo test --workspace
  ```

  When touching feature gates, also run the feature matrix from
  `.github/workflows/test.yml`.

## Benchmarks

- `bench/` compares release versions with the working tree. Config is
  `bench/bench.toml`; real datasets come from
  `~/Documents/segul-test-datasets` and follow the segul-bench paper setup
  (loci folders for concat/convert/summary, concatenated alignments for
  split).
- Run with `uv run segul-bench run` from `bench/`. Full runs take a long time;
  run them in the background. `bench/figures/` and `bench/.cache/` are not
  tracked.
- Python code there uses type hints, short docstrings, and only the standard
  library plus matplotlib.

## Changelog and releases

- User-visible changes go in `CHANGELOG.md` under the upcoming version, in
  the existing sections (Installation, Library Changes, Migration, Breaking
  Changes, ...). Write plain sentences describing the effect for users.
- The workspace version in the root `Cargo.toml` is the single source of
  truth for both crates. `segul-core` must be published before `segul`.
- Releases are titled with the tag name. Binaries are built by
  `release.yml` when a tag is pushed.
