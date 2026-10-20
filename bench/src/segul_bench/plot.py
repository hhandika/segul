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
from matplotlib.ticker import FuncFormatter, LogLocator, MaxNLocator, NullFormatter  # noqa: E402

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
    written.append(plot_combined(report, versions, colors, out_dir))
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


# Marker shapes per version, so versions are not told apart by color alone.
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "<", ">", "h"]


def plot_combined(report, versions, colors, out_dir) -> Path:
    """Time against memory for every task in one figure, as small multiples.

    Each panel is a task with one point per version (mean, with SD bars), and
    a line in version order shows how the task moved. Down and left is better.
    Panels have their own scales because tasks differ in size; time uses a
    log scale like the segul-bench paper figures.
    """
    tasks = report["tasks"]
    ncols = 3
    nrows = -(-len(tasks) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 0.6 + 2.9 * nrows), squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    markers = {v: MARKERS[i % len(MARKERS)] for i, v in enumerate(versions)}
    for ax, task in zip(axes.flat, tasks):
        style_axes(ax)
        results = {r["version"]: r for r in task["results"]}
        points = []
        for version in versions:
            summary = results.get(version, {}).get("summary", {})
            wall, rss = summary.get("wall_s"), summary.get("max_rss_mb")
            if wall and rss:
                points.append((version, wall, rss))
        if len(points) > 1:
            ax.plot([p[1]["mean"] for p in points], [p[2]["mean"] for p in points],
                    color=GRID, linewidth=1.5, zorder=1)
        for version, wall, rss in points:
            # Older versions are drawn larger and underneath,
            # so versions with nearly equal results stay visible.
            size = 8 + 2.5 * (len(versions) - 1 - versions.index(version))
            ax.errorbar(
                wall["mean"], rss["mean"], xerr=wall["sd"], yerr=rss["sd"],
                fmt=markers[version], markersize=size, color=colors[version],
                markeredgecolor=SURFACE, markeredgewidth=1.5,
                ecolor=INK_SECONDARY, elinewidth=0.8, capsize=0, zorder=2,
            )
        ax.set_xscale("log")
        set_log_ticks(ax, [p[1]["mean"] for p in points])
        ax.set_ylim(bottom=0, top=max((p[2]["mean"] + p[2]["sd"] for p in points), default=1) * 1.15)
        ax.margins(x=0.25)
        ax.grid(True, which="major", color=GRID, linewidth=0.8)
        ax.set_title(task["label"], color=INK, fontsize=9.5)
    for ax in list(axes.flat)[len(tasks):]:
        ax.set_visible(False)
    fig.supxlabel("Mean wall time (s, log scale)", color=INK_SECONDARY, fontsize=10, y=0.035)
    fig.supylabel("Mean peak memory (MB)", color=INK_SECONDARY, fontsize=10)
    fig.suptitle("Time and memory by version (lower left is better)", x=0.01, ha="left", color=INK, fontsize=12)
    handles = [
        plt.Line2D([], [], linestyle="none", marker=markers[v], markersize=8,
                   color=colors[v], markeredgecolor=SURFACE)
        for v in versions
    ]
    legend = fig.legend(handles, versions, loc="upper right", ncol=len(versions),
                        frameon=False, fontsize=9, bbox_to_anchor=(0.99, 0.995))
    for text in legend.get_texts():
        text.set_color(INK)
    add_caption(fig, report)
    fig.tight_layout(rect=(0.01, 0.05, 1, 0.96))
    path = out_dir / "combined.png"
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return path


def set_log_ticks(ax, values: list[float]) -> None:
    """Plain-number ticks on a log axis (0.5, 1, 2 instead of 5×10⁻¹).
    Narrow ranges get 1-2-5 steps, wide ones only powers of ten."""
    values = [v for v in values if v > 0]
    if not values:
        return
    spread = max(values) / min(values)
    if spread < 3:
        # Too narrow for log steps: a few evenly spaced ticks.
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=2))
    else:
        ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5) if spread < 100 else (1,)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())


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
