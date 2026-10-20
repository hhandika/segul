"""PNG charts of a benchmark report, for checking results by eye.

Versions are ordered, so they use one blue ramp from light (oldest) to dark
(newest). The ramp steps were validated as an ordinal palette on the light
chart surface.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e1e0d9"
# Blue steps 250 -> 700. Ordinal ramps start no lighter than step 250.
BLUE_RAMP = [
    "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6",
    "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b",
]
METRICS = {
    "wall_s": ("Wall time", "s"),
    "max_rss_mb": ("Peak memory (RSS)", "MB"),
}


def version_colors(versions: list[str]) -> dict[str, str]:
    n = len(versions)
    if n == 1:
        return {versions[0]: BLUE_RAMP[-1]}
    last = len(BLUE_RAMP) - 1
    return {v: BLUE_RAMP[round(i * last / (n - 1))] for i, v in enumerate(versions)}


def style_axes(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_SECONDARY, labelsize=9)
    ax.set_axisbelow(True)


def mean_of(result: dict, metric: str) -> float | None:
    summary = result.get("summary", {}).get(metric)
    return summary["mean"] if summary else None


def plot_report(report: dict, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    versions = [v["id"] for v in report["versions"]]
    colors = version_colors(versions)
    written = []
    for metric in METRICS:
        written.append(plot_relative(report, metric, versions, colors, out_dir))
    for task in report["tasks"]:
        written.append(plot_task(task, versions, colors, out_dir))
    return written


def plot_relative(report, metric, versions, colors, out_dir) -> Path:
    """Each task's metric relative to the oldest version that ran it,
    so tasks of very different sizes share one axis."""
    title, _ = METRICS[metric]
    tasks = report["tasks"]
    n_versions = len(versions)
    bar_h = 0.8 / n_versions
    fig, ax = plt.subplots(figsize=(9, 0.5 + 0.55 * len(tasks) * max(1, n_versions / 2.5)))
    fig.patch.set_facecolor(SURFACE)
    style_axes(ax)
    max_ratio = 1.0
    for t_idx, task in enumerate(tasks):
        results = {r["version"]: r for r in task["results"]}
        base = next((mean_of(results[v], metric) for v in versions if v in results and mean_of(results[v], metric)), None)
        if not base:
            continue
        for v_idx, version in enumerate(versions):
            value = mean_of(results.get(version, {}), metric)
            if value is None:
                continue
            ratio = value / base
            max_ratio = max(max_ratio, ratio)
            y = t_idx + (v_idx - (n_versions - 1) / 2) * bar_h
            ax.barh(y, ratio, height=bar_h, color=colors[version], edgecolor=SURFACE, linewidth=2)
            if v_idx == n_versions - 1 and n_versions > 1:
                better = "faster" if metric == "wall_s" else "less memory"
                label = f"{ratio:.2f}×" + (f"  ({1 / ratio:.1f}× {better})" if ratio < 0.95 else "")
                ax.text(ratio, y, "  " + label, va="center", fontsize=8, color=INK)
    ax.axvline(1.0, color=INK_SECONDARY, linewidth=1, linestyle="--")
    ax.set_yticks(range(len(tasks)), [t["label"] for t in tasks], color=INK)
    ax.invert_yaxis()
    ax.set_xlim(0, max_ratio * 1.35)
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_xlabel(f"{title} relative to the oldest version (lower is better)", color=INK_SECONDARY)
    fig.suptitle(f"{title} by version", x=0.01, ha="left", color=INK, fontsize=12)
    add_legend(ax, versions, colors)
    add_caption(fig, report)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    path = out_dir / ("time.png" if metric == "wall_s" else "memory.png")
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return path


def plot_task(task, versions, colors, out_dir) -> Path:
    """Absolute values for one task: time and memory side by side,
    each on its own axis."""
    results = {r["version"]: r for r in task["results"]}
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    fig.patch.set_facecolor(SURFACE)
    for ax, (metric, (title, unit)) in zip(axes, METRICS.items()):
        style_axes(ax)
        xs, heights, errs, bar_colors, labels = [], [], [], [], []
        for idx, version in enumerate(versions):
            summary = results.get(version, {}).get("summary", {}).get(metric)
            if not summary:
                continue
            xs.append(idx)
            heights.append(summary["mean"])
            errs.append(summary["sd"])
            bar_colors.append(colors[version])
            labels.append(version)
        ax.bar(xs, heights, yerr=errs, color=bar_colors, width=0.6, edgecolor=SURFACE, linewidth=2,
               error_kw={"ecolor": INK_SECONDARY, "elinewidth": 1, "capsize": 3})
        for x, h, err in zip(xs, heights, errs):
            ax.annotate(
                f"{h:.3g}", (x, h + err), xytext=(0, 4), textcoords="offset points",
                ha="center", va="bottom", fontsize=8, color=INK,
            )
        ax.set_xticks(xs, labels, color=INK)
        ax.yaxis.grid(True, color=GRID, linewidth=0.8)
        ax.set_ylabel(unit, color=INK_SECONDARY)
        ax.set_title(title, loc="left", color=INK, fontsize=10)
        ax.margins(y=0.15)
    fig.suptitle(f"{task['label']}  [{task['dataset']}]", x=0.01, ha="left", color=INK, fontsize=11)
    fig.tight_layout()
    path = out_dir / f"{slug(task['id'])}.png"
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return path


def add_legend(ax, versions, colors) -> None:
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[v]) for v in versions]
    # Above the plot area, so it never covers bars or labels.
    legend = ax.legend(
        handles, versions, loc="lower right", bbox_to_anchor=(1.0, 1.0),
        ncol=len(versions), frameon=False, fontsize=9,
    )
    for text in legend.get_texts():
        text.set_color(INK)


def add_caption(fig, report) -> None:
    machine = report["machine"]
    caption = (
        f"{machine.get('cpu', '')} · {machine.get('cores')} cores · {machine.get('os')} {machine.get('arch')} · "
        f"{report['config'].get('replicates')} replicates · scale {report['config'].get('scale')} · {report['created'][:10]}"
    )
    fig.text(0.01, 0.005, caption, fontsize=7, color=INK_SECONDARY, ha="left", va="bottom")


def slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "-", text)
