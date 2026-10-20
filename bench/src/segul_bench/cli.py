"""segul-bench: track SEGUL execution time and memory usage across versions."""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

from . import datasets as ds
from .config import DEFAULT_CONFIG, FIGURES_DIR, RESULTS_DIR, SCALES, load_config
from .results import (
    build_report,
    load_report,
    machine_info,
    print_table,
    result_to_json,
    write_report,
)
from .runner import run_task, timer_binary
from .versions import fetch_release, resolve

LATEST = RESULTS_DIR / "latest.json"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="segul-bench", description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help="benchmark config (default: bench.toml)")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="run the benchmark and write JSON results")
    run.add_argument("--versions", help="comma-separated, oldest first, e.g. v0.22.1,v0.23.2,current")
    run.add_argument("--tasks", help="comma-separated task ids (default: all)")
    run.add_argument("--replicates", type=int, help="recorded runs per version")
    run.add_argument("--warmup", type=int, help="unrecorded runs per version before measuring")
    run.add_argument("--scale", choices=SCALES, help="synthetic dataset size")
    run.add_argument("--label", help="machine label stored in the results, e.g. 'workstation'")
    run.add_argument("--skip-build", action="store_true", help="use the existing target/release/segul")
    run.add_argument("--keep-output", action="store_true", help="keep task outputs in .cache/work")
    run.add_argument("--no-plot", action="store_true", help="do not write PNG charts")
    run.add_argument("--fail-on-regression", action="store_true", help="exit with 1 if current regresses")
    run.add_argument("--output-dir", type=Path, default=RESULTS_DIR, help="where to write JSON results")

    data = sub.add_parser("data", help="generate the synthetic datasets")
    data.add_argument("--scale", choices=SCALES)
    data.add_argument("--datasets", help="comma-separated dataset ids (default: all)")

    fetch = sub.add_parser("fetch", help="download release binaries")
    fetch.add_argument("tags", nargs="+", help="release tags, e.g. v0.23.2")
    fetch.add_argument("--refresh", action="store_true", help="download again even if cached")

    report = sub.add_parser("report", help="print the comparison table of a results file")
    report.add_argument("results", type=Path, nargs="?", default=LATEST)

    plot = sub.add_parser("plot", help="write PNG charts of a results file")
    plot.add_argument("results", type=Path, nargs="?", default=LATEST)
    plot.add_argument("--output-dir", type=Path, default=FIGURES_DIR)

    export = sub.add_parser("export", help="copy a results file for the website")
    export.add_argument("destination", type=Path, help="e.g. ../segui/website/src/data/segul-versions.json")
    export.add_argument("--results", type=Path, default=LATEST)

    args = parser.parse_args(argv)
    # Show progress right away when output goes to a file or pipe.
    sys.stdout.reconfigure(line_buffering=True)
    config = load_config(args.config)
    handlers = {"run": cmd_run, "data": cmd_data, "fetch": cmd_fetch, "report": cmd_report, "plot": cmd_plot, "export": cmd_export}
    return handlers[args.command](args, config) or 0


def split(value: str | None) -> list[str] | None:
    return [v.strip() for v in value.split(",") if v.strip()] if value else None


def cmd_run(args, config) -> int:
    if args.replicates is not None:
        config.replicates = args.replicates
    if args.warmup is not None:
        config.warmup = args.warmup
    scale = args.scale or config.scale
    version_ids = split(args.versions) or config.versions
    task_ids = split(args.tasks)
    tasks = [t for t in config.tasks if task_ids is None or t.id in task_ids]
    if task_ids:
        unknown = set(task_ids) - {t.id for t in tasks}
        if unknown:
            raise SystemExit(f"Unknown task ids: {', '.join(sorted(unknown))}")

    versions = resolve(version_ids, skip_build=args.skip_build)
    timer = timer_binary()
    task_reports, dataset_reports = [], {}
    for idx, task in enumerate(tasks, 1):
        dataset = config.datasets[task.dataset]
        input_dir = ds.ensure(dataset, scale)
        print(f"\n[{idx}/{len(tasks)}] {task.label} ({task.id})")
        if input_dir is None:
            print(f"  skipped: dataset {dataset.id} not found at {dataset.path}")
            continue
        dataset_reports.setdefault(dataset.id, ds.describe(dataset, input_dir))
        results = run_task(timer, task, dataset, versions, input_dir, config, keep_output=args.keep_output)
        task_reports.append(
            {
                "id": task.id,
                "label": task.label,
                "dataset": task.dataset,
                "args": task.args,
                "results": [result_to_json(r) for r in results],
            }
        )

    report = build_report(
        machine=machine_info(args.label),
        config={"replicates": config.replicates, "warmup": config.warmup, "scale": scale},
        versions=[v.to_json() for v in versions],
        datasets=list(dataset_reports.values()),
        tasks=task_reports,
    )
    path = write_report(report, args.output_dir)
    regressions = print_table(report, config.regression_threshold)
    print(f"\nResults: {path}")
    if not args.no_plot:
        plot_files(report, FIGURES_DIR)
    return 1 if regressions and args.fail_on_regression else 0


def cmd_data(args, config) -> int:
    scale = args.scale or config.scale
    ids = split(args.datasets)
    for dataset in config.datasets.values():
        if ids and dataset.id not in ids:
            continue
        path = ds.ensure(dataset, scale)
        if path is None:
            print(f"{dataset.id}: not found at {dataset.path}")
        else:
            info = ds.describe(dataset, path)
            print(f"{dataset.id}: {info['files']} files, {info['size_mb']} MB ({path})")
    return 0


def cmd_fetch(args, config) -> int:
    for tag in args.tags:
        version = fetch_release(tag, refresh=args.refresh)
        print(f"{tag}: {version.binary}")
    return 0


def cmd_report(args, config) -> int:
    print_table(load_report(args.results), config.regression_threshold)
    return 0


def cmd_plot(args, config) -> int:
    plot_files(load_report(args.results), args.output_dir)
    return 0


def cmd_export(args, config) -> int:
    load_report(args.results)  # checks the schema version
    args.destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(args.results, args.destination)
    print(f"Copied {args.results} to {args.destination}")
    return 0


def plot_files(report: dict, out_dir: Path) -> None:
    from .plot import plot_report  # matplotlib is only needed here

    files = plot_report(report, out_dir)
    print(f"Charts: {len(files)} PNG files in {out_dir}")


if __name__ == "__main__":
    sys.exit(main())
