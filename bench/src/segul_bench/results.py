"""Benchmark results: JSON output, summaries and the comparison table."""

from __future__ import annotations

import json
import os
import platform
import statistics
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from .runner import TaskVersionResult

SCHEMA_VERSION = 1
METRICS = ("wall_s", "max_rss_mb")


def summarize(values: list[float]) -> dict:
    if not values:
        return {}
    return {
        "mean": round(statistics.fmean(values), 4),
        "sd": round(statistics.stdev(values), 4) if len(values) > 1 else 0.0,
        "min": round(min(values), 4),
        "median": round(statistics.median(values), 4),
        "n": len(values),
    }


def result_to_json(result: TaskVersionResult) -> dict:
    data: dict = {"version": result.version}
    if result.skipped:
        data["skipped"] = result.skipped
        return data
    ok = [r for r in result.replicates if r.ok]
    data["replicates"] = [r.to_json() for r in result.replicates]
    data["summary"] = {metric: summarize([getattr(r, metric) for r in ok]) for metric in METRICS}
    data["output_hash"] = result.output_hash
    if result.error:
        data["error"] = result.error
    return data


def machine_info(label: str | None) -> dict:
    info = {
        "os": platform.system(),
        "os_release": platform.release(),
        "arch": platform.machine(),
        "cpu": cpu_model(),
        "cores": os.cpu_count(),
        "ram_gb": ram_gb(),
    }
    if label:
        info["label"] = label
    return info


def cpu_model() -> str:
    if platform.system() == "Darwin":
        return sysctl("machdep.cpu.brand_string") or platform.processor()
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor()


def ram_gb() -> float | None:
    if platform.system() == "Darwin":
        value = sysctl("hw.memsize")
        return round(int(value) / 1024**3, 1) if value else None
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemTotal"):
                return round(int(line.split()[1]) / 1024**2, 1)
    except OSError:
        pass
    return None


def sysctl(name: str) -> str | None:
    try:
        return subprocess.run(["sysctl", "-n", name], capture_output=True, text=True).stdout.strip()
    except OSError:
        return None


def build_report(machine, config, versions, datasets, tasks) -> dict:
    return {
        "schema_version": SCHEMA_VERSION,
        "created": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "machine": machine,
        "config": config,
        "versions": versions,
        "datasets": datasets,
        "tasks": tasks,
    }


def write_report(report: dict, results_dir: Path) -> Path:
    results_dir.mkdir(parents=True, exist_ok=True)
    machine = report["machine"]
    stamp = report["created"][:10]
    name = f"{stamp}-{machine['os'].lower()}-{machine['arch']}.json"
    path = results_dir / name
    text = json.dumps(report, indent=2) + "\n"
    path.write_text(text)
    (results_dir / "latest.json").write_text(text)
    return path


def load_report(path: Path) -> dict:
    report = json.loads(path.read_text())
    if report.get("schema_version") != SCHEMA_VERSION:
        raise SystemExit(f"{path} has schema version {report.get('schema_version')}, expected {SCHEMA_VERSION}")
    return report


def print_table(report: dict, threshold: float = 0.05) -> list[str]:
    """Print each task's results per version, relative to the first version.
    Returns the regressions of `current` against the version before it."""
    regressions = []
    for task in report["tasks"]:
        print(f"\n{task['label']}  [{task['dataset']}]")
        print(f"  {'version':>10} {'time (s)':>16} {'speedup':>8} {'peak RSS (MB)':>18} {'memory':>7}  same output")
        results = [r for r in task["results"] if not r.get("skipped") and r.get("summary", {}).get("wall_s")]
        if not results:
            print("  no results")
            continue
        base = results[0]
        for r in task["results"]:
            if r.get("skipped"):
                print(f"  {r['version']:>10}  skipped: {r['skipped']}")
                continue
            wall = r["summary"].get("wall_s")
            rss = r["summary"].get("max_rss_mb")
            if not wall:
                print(f"  {r['version']:>10}  failed: {r.get('error', '').splitlines()[-1:] or ''}")
                continue
            speedup = base["summary"]["wall_s"]["mean"] / wall["mean"]
            mem_ratio = rss["mean"] / base["summary"]["max_rss_mb"]["mean"]
            same = "yes" if r.get("output_hash") == base.get("output_hash") else "NO"
            print(
                f"  {r['version']:>10} {wall['mean']:>9.3f} ±{wall['sd']:<6.3f}"
                f"{speedup:>7.2f}× {rss['mean']:>11.1f} ±{rss['sd']:<6.1f}"
                f"{mem_ratio:>6.2f}×  {same}"
            )
        current = next((r for r in results if r["version"] == "current"), None)
        previous = [r for r in results if r["version"] != "current"]
        if current and previous:
            prev = previous[-1]
            for metric, name in (("wall_s", "time"), ("max_rss_mb", "memory")):
                ratio = current["summary"][metric]["mean"] / prev["summary"][metric]["mean"]
                if ratio > 1 + threshold:
                    msg = f"{task['id']}: {name} {ratio:.2f}× of {prev['version']}"
                    regressions.append(msg)
    if regressions:
        print("\nRegressions in current (more than {:.0%} worse than the previous version):".format(threshold))
        for msg in regressions:
            print(f"  - {msg}")
    return regressions
