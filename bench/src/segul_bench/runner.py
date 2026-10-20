"""Run tasks and measure wall time, CPU time and peak memory.

Commands run through `segul-bench-timer` (bench/timer), a small native helper
that spawns the command and reports its resource usage. Measuring from
Python directly would be wrong: Linux keeps the memory high-water mark of
the spawning process across `exec`, so every command would report at least
the memory of the Python runner.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from .config import BENCH_DIR, CACHE_DIR, Config, Dataset, Task
from .versions import Version, parse_version

# Output files that differ between runs for reasons unrelated to results.
IGNORED_OUTPUTS = {".log"}


@dataclass
class Replicate:
    wall_s: float
    user_s: float
    sys_s: float
    max_rss_mb: float
    ok: bool
    exit_code: int

    def to_json(self) -> dict:
        return {
            "wall_s": round(self.wall_s, 4),
            "user_s": round(self.user_s, 4),
            "sys_s": round(self.sys_s, 4),
            "max_rss_mb": round(self.max_rss_mb, 2),
            "ok": self.ok,
            "exit_code": self.exit_code,
        }


@dataclass
class TaskVersionResult:
    version: str
    replicates: list[Replicate] = field(default_factory=list)
    output_hash: str | None = None
    skipped: str | None = None
    error: str | None = None


def timer_binary() -> Path:
    """Build the timer helper once (it is tiny) and return its path."""
    target_dir = CACHE_DIR / "timer-target"
    binary = target_dir / "release" / "segul-bench-timer"
    subprocess.run(
        [
            "cargo", "build", "--release", "--quiet",
            "--manifest-path", str(BENCH_DIR / "timer" / "Cargo.toml"),
            "--target-dir", str(target_dir),
        ],
        check=True,
    )
    return binary


def measure(timer: Path, args: list[str], cwd: Path, stderr_path: Path) -> Replicate:
    with open(stderr_path, "wb") as stderr:
        proc = subprocess.run([str(timer), *args], cwd=cwd, stdout=subprocess.PIPE, stderr=stderr)
    if proc.returncode != 0:
        raise SystemExit(f"segul-bench-timer failed to run {args[0]}, see {stderr_path}")
    usage = json.loads(proc.stdout)
    return Replicate(
        wall_s=usage["wall_s"],
        user_s=usage["user_s"],
        sys_s=usage["sys_s"],
        max_rss_mb=usage["max_rss_mb"],
        ok=usage["exit_code"] == 0,
        exit_code=usage["exit_code"],
    )


def output_hash(out_dir: Path) -> str | None:
    """Hash output contents, ignoring line order.

    Each line is hashed on its own and the hashes are summed per file, so
    outputs whose rows only come in a different order (e.g. taxon order)
    hash the same. Files are combined by their relative path.
    """
    if not out_dir.exists():
        return None
    total = hashlib.sha256()
    for path in sorted(p for p in out_dir.rglob("*") if p.is_file()):
        if path.suffix in IGNORED_OUTPUTS:
            continue
        line_sum = 0
        with open(path, "rb") as f:
            for line in f:
                digest = hashlib.blake2b(line.rstrip(b"\r\n"), digest_size=8).digest()
                line_sum = (line_sum + int.from_bytes(digest, "little")) % (1 << 64)
        total.update(str(path.relative_to(out_dir)).encode())
        total.update(line_sum.to_bytes(8, "little"))
    return total.hexdigest()[:16]


def build_args(task: Task, dataset: Dataset, binary: Path, input_dir: Path, out_dir: Path) -> list[str]:
    values = {"input": str(input_dir), "output": str(out_dir), "datatype": dataset.datatype}
    values.update({key: value.format(**values) for key, value in dataset.vars.items()})
    try:
        return [str(binary)] + [arg.format(**values) for arg in task.args]
    except KeyError as err:
        raise SystemExit(f"Task {task.id}: unknown placeholder {err} for dataset {dataset.id}") from err


def skip_reason(task: Task, version: Version) -> str | None:
    if task.min_version and version.number < parse_version(task.min_version):
        return f"requires {task.min_version} or later"
    return None


def run_task(
    timer: Path,
    task: Task,
    dataset: Dataset,
    versions: list[Version],
    input_dir: Path,
    config: Config,
    keep_output: bool = False,
) -> list[TaskVersionResult]:
    results = {v.id: TaskVersionResult(version=v.id) for v in versions}
    active = []
    for version in versions:
        reason = skip_reason(task, version)
        if reason:
            results[version.id].skipped = reason
            print(f"  {version.id}: skipped ({reason})")
        else:
            active.append(version)

    work_root = CACHE_DIR / "work" / task.id

    def run_once(version: Version) -> Replicate:
        work = work_root / version.id
        out_dir = work / "out"
        if out_dir.exists():
            shutil.rmtree(out_dir)
        work.mkdir(parents=True, exist_ok=True)
        args = build_args(task, dataset, version.binary, input_dir, out_dir)
        return measure(timer, args, work, work / "stderr.txt")

    for version in active:
        for _ in range(config.warmup):
            rep = run_once(version)
            if not rep.ok:
                results[version.id].error = read_error(work_root / version.id / "stderr.txt")
                print(f"  {version.id}: failed (exit {rep.exit_code}), see JSON for the error")
                break

    runnable = [v for v in active if results[v.id].error is None]
    for rep_idx in range(config.replicates):
        # Alternate the order to cancel drift such as thermal throttling.
        order = runnable if rep_idx % 2 == 0 else list(reversed(runnable))
        for version in order:
            rep = run_once(version)
            results[version.id].replicates.append(rep)
            status = "" if rep.ok else f" FAILED (exit {rep.exit_code})"
            print(
                f"  {version.id:>10} rep {rep_idx + 1}/{config.replicates}: "
                f"{rep.wall_s:7.3f} s {rep.max_rss_mb:9.1f} MB{status}"
            )
            if not rep.ok and results[version.id].error is None:
                results[version.id].error = read_error(work_root / version.id / "stderr.txt")

    for version in runnable:
        results[version.id].output_hash = output_hash(work_root / version.id / "out")
    if not keep_output and work_root.exists():
        shutil.rmtree(work_root)
    return [results[v.id] for v in versions]


def read_error(path: Path, max_lines: int = 5) -> str:
    try:
        lines = path.read_text(errors="replace").strip().splitlines()
    except OSError:
        return ""
    return "\n".join(lines[-max_lines:])
