# segul-bench

Track SEGUL execution time and peak memory across versions. `segul-bench` runs
the same tasks with previous releases and the current working tree, then
writes:

- **JSON** results in `results/` (committed), used by the website to show
  performance across versions.
- **PNG** charts in `figures/` (ignored by git), for checking results by eye.

## Requirements

- [uv](https://docs.astral.sh/uv/)
- Rust toolchain (`cargo`), to build the current version and the timer helper
- Linux or macOS (x86_64 or arm64), the platforms with release binaries

## Usage

Run from this directory, or from anywhere with `uv run --project <path-to-segul>/bench`.

```bash
uv run segul-bench run
```

This builds the working tree, downloads the release binaries listed in
`bench.toml`, generates the synthetic datasets on first use, runs every task,
prints a comparison table, writes `results/<date>-<os>-<arch>.json` and
`results/latest.json`, and draws the charts.

Common options:

```bash
# Quick check: small datasets, fewer replicates, existing build
uv run segul-bench run --scale small --replicates 2 --skip-build

# Choose versions (oldest first) and tasks
uv run segul-bench run --versions v0.21.3,v0.23.2,current --tasks concat-nexus,contig-summary

# Exit with 1 if current is more than 5% slower or larger than the previous version
uv run segul-bench run --fail-on-regression
```

Other commands:

```bash
uv run segul-bench data --scale full        # generate the synthetic datasets only
uv run segul-bench fetch v0.22.1 v0.23.2    # download release binaries only
uv run segul-bench report                   # print the table of results/latest.json
uv run segul-bench plot                     # draw charts of results/latest.json
uv run segul-bench export <path-to-segui>/website/src/data/segul-versions.json
```

## How it measures

- **Time:** wall time, plus user and system CPU time.
- **Memory:** peak resident memory (`ru_maxrss`) of the segul process.
- Commands run through `segul-bench-timer` (`timer/`), a small Rust helper
  built automatically with cargo. Measuring from Python directly would be wrong:
  Linux keeps the memory high-water mark of the spawning process across `exec`,
  so every command would report at least the memory of the Python runner.
- Each version gets one warm-up run (not recorded, it warms the file cache),
  then `replicates` recorded runs. Version order alternates between replicates
  to cancel drift such as CPU throttling.
- **Same output:** each version's output files are hashed per line and
  ignoring line order, so the table shows whether versions produce the same
  results even if they order taxa differently. Log files are ignored.

Results depend on the machine. Compare versions within one results file, not
across machines. Close other heavy programs while running.

## Configuration

`bench.toml` sets the default versions, replicates, scale, datasets and tasks.

Synthetic datasets are generated with a fixed seed into
`.cache/data/<scale>/<id>/`, so they are the same on every machine. Use
`--scale small` for quick checks and `full` (0.3–1 GB per dataset) for results.

To benchmark real data, add a dataset with a local path. Tasks using it are
skipped when the path does not exist:

```toml
[[datasets]]
id = "oliveros-2019"
kind = "real"
path = "~/data/oliveros-2019/nexus"
description = "Oliveros et al. 2019 (522.5 MBases, DNA)"

[[tasks]]
id = "concat-oliveros"
label = "Concatenate Oliveros et al. 2019"
dataset = "oliveros-2019"
args = ["align", "concat", "-d", "{input}", "-f", "nexus", "-o", "{output}"]
```

In `args`, `{input}` is the dataset directory and `{output}` a fresh output
directory for each run. Add `min_version = "0.24.0"` to skip releases that
do not have a command.

## Results format

```json
{
  "schema_version": 1,
  "created": "2026-10-20T10:00:00Z",
  "machine": {"os": "Linux", "arch": "x86_64", "cpu": "...", "cores": 24, "ram_gb": 62.7},
  "config": {"replicates": 5, "warmup": 1, "scale": "full"},
  "versions": [
    {"id": "v0.23.2", "source": "release", "version": "0.23.2"},
    {"id": "current", "source": "build", "version": "0.24.0", "git": "v0.23.2-45-g4b30607-dirty"}
  ],
  "datasets": [{"id": "loci-nexus", "kind": "synthetic", "description": "...", "size_mb": 276.1, "files": 2000}],
  "tasks": [
    {
      "id": "concat-nexus",
      "label": "Concatenate NEXUS loci",
      "dataset": "loci-nexus",
      "args": ["align", "concat", "..."],
      "results": [
        {
          "version": "v0.23.2",
          "replicates": [{"wall_s": 1.23, "user_s": 1.51, "sys_s": 0.2, "max_rss_mb": 374.1, "ok": true, "exit_code": 0}],
          "summary": {
            "wall_s": {"mean": 1.23, "sd": 0.01, "min": 1.22, "median": 1.23, "n": 5},
            "max_rss_mb": {"mean": 374.3, "sd": 0.3, "min": 374.1, "median": 374.2, "n": 5}
          },
          "output_hash": "674cd834a11f1e2b"
        },
        {"version": "v0.21.3", "skipped": "requires 0.22.0 or later"}
      ]
    }
  ]
}
```

`summary` only includes successful replicates. A result has `error` with the
last lines of stderr when a version fails, and `skipped` instead of
measurements when `min_version` excludes it. Bump `schema_version` when the
format changes.
