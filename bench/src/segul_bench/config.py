"""Paths and the benchmark configuration (bench.toml)."""

from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from pathlib import Path

# src/segul_bench/config.py -> bench/
BENCH_DIR = Path(__file__).resolve().parents[2]
REPO_DIR = BENCH_DIR.parent
CACHE_DIR = BENCH_DIR / ".cache"
RESULTS_DIR = BENCH_DIR / "results"
FIGURES_DIR = BENCH_DIR / "figures"
DEFAULT_CONFIG = BENCH_DIR / "bench.toml"

SCALES = ("small", "full")


@dataclass
class Dataset:
    id: str
    kind: str  # "synthetic" or "real"
    description: str
    generator: str | None = None
    params: dict[str, dict] = field(default_factory=dict)
    path: Path | None = None


@dataclass
class Task:
    id: str
    label: str
    dataset: str
    args: list[str]
    min_version: str | None = None


@dataclass
class Config:
    versions: list[str]
    replicates: int
    warmup: int
    scale: str
    regression_threshold: float
    datasets: dict[str, Dataset]
    tasks: list[Task]


def load_config(path: Path = DEFAULT_CONFIG) -> Config:
    with open(path, "rb") as f:
        raw = tomllib.load(f)
    defaults = raw.get("defaults", {})
    datasets = {}
    for d in raw.get("datasets", []):
        dataset_path = Path(d["path"]).expanduser() if "path" in d else None
        datasets[d["id"]] = Dataset(
            id=d["id"],
            kind=d.get("kind", "synthetic"),
            description=d.get("description", ""),
            generator=d.get("generator"),
            params=d.get("params", {}),
            path=dataset_path,
        )
    tasks = [
        Task(
            id=t["id"],
            label=t["label"],
            dataset=t["dataset"],
            args=list(t["args"]),
            min_version=t.get("min_version"),
        )
        for t in raw.get("tasks", [])
    ]
    for task in tasks:
        if task.dataset not in datasets:
            raise SystemExit(f"Task {task.id} uses unknown dataset {task.dataset}")
    return Config(
        versions=list(defaults.get("versions", ["current"])),
        replicates=int(defaults.get("replicates", 5)),
        warmup=int(defaults.get("warmup", 1)),
        scale=defaults.get("scale", "full"),
        regression_threshold=float(defaults.get("regression_threshold", 0.05)),
        datasets=datasets,
        tasks=tasks,
    )


def workspace_version() -> str:
    """The version in the root Cargo.toml, used as the version of `current`."""
    with open(REPO_DIR / "Cargo.toml", "rb") as f:
        cargo = tomllib.load(f)
    return cargo["workspace"]["package"]["version"]
