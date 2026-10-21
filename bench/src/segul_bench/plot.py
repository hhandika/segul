"""PNG charts of a benchmark report, for checking results by eye.

- scorecard.png: the newest version relative to the oldest, per task, as a
  table colored with the ColorBrewer RdBu diverging scale (blue is better).
- scatter.png: time against memory, one panel per analysis, one color and
  marker per dataset (ColorBrewer Dark2), with arrows from older versions
  to the newest.
- <task>.png: time and memory of one task per version (ColorBrewer Blues,
  light for older versions and dark for newer).
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import colormaps  # noqa: E402
from matplotlib.colors import TwoSlopeNorm, to_hex, to_rgb  # noqa: E402
from matplotlib.ticker import FuncFormatter, LogLocator, MaxNLocator, NullFormatter  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e1e0d9"
METRICS = {
    "wall_s": ("Wall time", "s"),
    "max_rss_mb": ("Peak memory", "MB"),
}
# A change smaller than this counts as "no change".
NO_CHANGE = 0.05
# Short names for segul subcommands, used to group tasks into analyses.
ANALYSES = {
    "align concat": "Concatenate",
    "align convert": "Convert",
    "align summary": "Summarize",
    "align split": "Split",
    "align filter": "Filter",
    "align partition": "Partition",
    "align translate": "Translate",
    "sequence id": "Sample IDs",
    "sequence remove": "Remove sequences",
    "sequence rename": "Rename sequences",
    "sequence extract": "Extract sequences",
    "contig summary": "Contig summary",
    "read summary": "Read summary",
}
# Marker shapes per dataset, so datasets are not told apart by color alone.
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "<", ">", "h", "p", "*"]


def plot_report(report: dict, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in ("time.png", "memory.png", "combined.png"):  # earlier chart names
        (out_dir / stale).unlink(missing_ok=True)
    versions = [v["id"] for v in report["versions"]]
    written = [plot_scorecard(report, versions, out_dir), plot_scatter(report, versions, out_dir)]
    colors = version_colors(versions)
    labels = {d["id"]: dataset_label(d) for d in report["datasets"]}
    for task in report["tasks"]:
        written.append(plot_task(task, versions, colors, labels.get(task["dataset"], task["dataset"]), out_dir))
    return written


# Data helpers -------------------------------------------------------------


def analysis_of(task: dict) -> str:
    command = " ".join(task.get("args", [])[:2])
    return ANALYSES.get(command, command.capitalize() or task["id"])


def means(task: dict) -> dict[str, tuple[float, float]]:
    """Mean wall time and peak memory per version that ran successfully."""
    values = {}
    for result in task["results"]:
        summary = result.get("summary", {})
        if summary.get("wall_s") and summary.get("max_rss_mb"):
            values[result["version"]] = (summary["wall_s"]["mean"], summary["max_rss_mb"]["mean"])
    return values


def grouped_rows(report: dict, versions: list[str]) -> list[tuple[str, list[tuple[str, dict]]]]:
    """Tasks grouped by analysis, in first-seen order, labeled by dataset.
    Tasks without results for at least one version are left out."""
    names = {d["id"]: dataset_label(d) for d in report["datasets"]}
    groups: dict[str, list[tuple[str, dict, dict]]] = {}
    for task in report["tasks"]:
        values = means(task)
        if values:
            groups.setdefault(analysis_of(task), []).append((names.get(task["dataset"], task["dataset"]), task, values))
    rows = []
    for analysis, items in groups.items():
        labels = [name for name, _, _ in items]
        # Fall back to the task label when one dataset has several tasks in a group.
        rows.append((analysis, [
            (name if labels.count(name) == 1 else task["label"], values) for name, task, values in items
        ]))
    return rows


def dataset_label(dataset: dict) -> str:
    """Dataset name with its total sequence length, e.g.
    "Esselstyn et al. 2021 (358.1 Mbp)"."""
    name = dataset.get("name") or dataset["id"]
    bases = dataset.get("bases")
    if not bases:
        return name
    size = f"{bases / 1e6:.1f} {'M aa' if dataset.get('datatype') == 'aa' else 'Mbp'}"
    return f"{name} ({size})"


def compare(values: dict, versions: list[str]) -> tuple[str, str] | None:
    """The oldest and newest versions with results."""
    present = [v for v in versions if v in values]
    if len(present) < 2:
        return None
    return present[0], present[-1]


# Scorecard ---------------------------------------------------------------


def plot_scorecard(report: dict, versions: list[str], out_dir: Path) -> Path:
    """One row per task: how much faster and smaller the newest version is
    than the oldest, colored on a diverging scale centered on no change."""
    rows = grouped_rows(report, versions)
    cmap = colormaps["RdBu"]
    # log(old / new): positive is better (blue), negative worse (red).
    norm = TwoSlopeNorm(vmin=-math.log(10), vcenter=0, vmax=math.log(10))
    n_lines = sum(len(items) + 1 for _, items in rows)
    fig, ax = plt.subplots(figsize=(7.6, 1.0 + 0.22 * n_lines))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    ax.axis("off")

    # Column x positions: label | time cell | memory cell | time values | memory values
    label_x, time_x, mem_x, cell_w, time_val_x, mem_val_x = 0.0, 2.7, 3.95, 1.2, 5.3, 6.45
    y = 0.0
    ax.text(time_x + cell_w / 2, y, "Time", ha="center", va="center", fontsize=9.5, color=INK_SECONDARY)
    ax.text(mem_x + cell_w / 2, y, "Memory", ha="center", va="center", fontsize=9.5, color=INK_SECONDARY)
    ax.text(time_val_x, y, "Time (s)", va="center", fontsize=9.5, color=INK_SECONDARY)
    ax.text(mem_val_x, y, "Peak memory (MB)", va="center", fontsize=9.5, color=INK_SECONDARY)
    compared = set()
    for analysis, items in rows:
        y += 1.15
        ax.text(label_x, y, analysis, va="center", fontsize=10, fontweight="bold", color=INK)
        for label, values in items:
            y += 1
            ax.text(label_x + 0.15, y, label, va="center", fontsize=9, color=INK)
            pair = compare(values, versions)
            if pair is None:
                ax.text(time_x, y, "needs two versions", va="center", fontsize=8.5, color=INK_SECONDARY)
                continue
            old_v, new_v = pair
            compared.add(pair)
            for idx, (x, val_x, unit) in enumerate(((time_x, time_val_x, "s"), (mem_x, mem_val_x, "MB"))):
                old, new = values[old_v][idx], values[new_v][idx]
                change = math.log(old / new)
                color = cmap(norm(max(min(change, norm.vmax), norm.vmin)))
                ax.add_patch(plt.Rectangle((x, y - 0.42), cell_w, 0.84, color=color, linewidth=0))
                ax.text(x + cell_w / 2, y, change_text(old / new, idx), ha="center", va="center",
                        fontsize=8.5, color=text_on(color))
                ax.text(val_x, y, f"{fmt(old)} → {fmt(new)}", va="center", fontsize=8.5, color=INK)
    ax.set_xlim(-0.1, 7.9)
    ax.set_ylim(y + 0.8, -0.6)

    pairs = ", ".join(f"{new} vs {old}" for old, new in sorted(compared))
    fig.suptitle(f"Scorecard: {pairs or 'no comparison'}", x=0.01, ha="left", color=INK, fontsize=12)
    fig.text(0.01, 0.945, "Blue is faster or smaller, red is slower or larger. Values are means, oldest → newest.",
             fontsize=8, color=INK_SECONDARY, va="top")
    add_caption(fig, report)
    fig.tight_layout(rect=(0, 0.02, 1, 0.93), pad=0.4)
    path = out_dir / "scorecard.png"
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return path


def change_text(ratio: float, metric_idx: int) -> str:
    """ratio = old / new, so above 1 means the newest version is better."""
    if abs(math.log(ratio)) < math.log(1 + NO_CHANGE):
        return "no change"
    if ratio > 1:
        return f"{ratio:.1f}× {'faster' if metric_idx == 0 else 'smaller'}"
    return f"{1 / ratio:.1f}× {'slower' if metric_idx == 0 else 'larger'}"


def text_on(color) -> str:
    """Dark text on light cells and light text on dark cells."""
    r, g, b = to_rgb(color)
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return INK if luminance > 0.5 else SURFACE


def fmt(value: float) -> str:
    if value >= 100:
        return f"{value:,.0f}"
    if value >= 10:
        return f"{value:.1f}"
    return f"{value:.2f}"


# Scatter -----------------------------------------------------------------


def plot_scatter(report: dict, versions: list[str], out_dir: Path) -> Path:
    """Time against memory, one panel per analysis. Each dataset has its own
    color and marker; hollow markers are older versions and the arrow points
    to the newest (filled). Lower left is better."""
    rows = grouped_rows(report, versions)
    datasets = []
    for _, items in rows:
        for label, _ in items:
            if label not in datasets:
                datasets.append(label)
    colors = dataset_colors(datasets)
    markers = {d: MARKERS[i % len(MARKERS)] for i, d in enumerate(datasets)}

    ncols = 2 if len(rows) <= 4 else 3
    nrows = max(1, -(-len(rows) // ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 3.9 * nrows + 1.4), squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    for ax, (analysis, items) in zip(axes.flat, rows):
        style_axes(ax)
        ax.set_xscale("log")
        ax.set_yscale("log")
        xs, ys = [], []
        for label, values in items:
            present = [v for v in versions if v in values]
            points = [values[v] for v in present]
            xs += [p[0] for p in points]
            ys += [p[1] for p in points]
            color = colors[label]
            for (x0, y0), (x1, y1) in zip(points, points[1:]):
                ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(
                    arrowstyle="-|>", color=color, alpha=0.55, lw=1.4, shrinkA=6, shrinkB=6, mutation_scale=11))
            for idx, (x, y) in enumerate(points):
                newest = idx == len(points) - 1
                ax.scatter(x, y, s=58, marker=markers[label], zorder=3,
                           facecolor=color if newest else SURFACE, edgecolor=color if not newest else SURFACE,
                           linewidth=1.6 if not newest else 1.0)
        set_log_ticks(ax.xaxis, xs)
        set_log_ticks(ax.yaxis, ys)
        ax.margins(0.15)
        ax.grid(True, which="major", color=GRID, linewidth=0.8)
        ax.set_title(analysis, loc="left", color=INK, fontsize=10.5)
    for ax in list(axes.flat)[len(rows):]:
        ax.set_visible(False)

    fig.supylabel("Peak memory (MB, log scale)", color=INK_SECONDARY, fontsize=10)
    fig.suptitle("Time and memory by analysis (lower left is better)", x=0.01, ha="left", color=INK, fontsize=12)

    dataset_handles = [
        plt.Line2D([], [], linestyle="none", marker=markers[d], markersize=8,
                   markerfacecolor=colors[d], markeredgecolor=SURFACE)
        for d in datasets
    ]
    older = [v for v in versions[:-1]]
    version_handles = [
        plt.Line2D([], [], linestyle="none", marker="o", markersize=8,
                   markerfacecolor=SURFACE, markeredgecolor=INK_SECONDARY, markeredgewidth=1.6),
        plt.Line2D([], [], linestyle="none", marker="o", markersize=8,
                   markerfacecolor=INK_SECONDARY, markeredgecolor=SURFACE),
    ]
    version_labels = [f"Older: {', '.join(older)}" if older else "Older", f"Newest: {versions[-1]}"]
    caption_h = 0.22 / fig.get_figheight()  # the caption line, in figure fraction
    legend = fig.legend(
        dataset_handles + version_handles, datasets + version_labels,
        loc="lower center", ncol=min(4, len(datasets) + 2), frameon=False, fontsize=8.5,
        bbox_to_anchor=(0.5, caption_h), columnspacing=1.6, handletextpad=0.4,
    )
    for text in legend.get_texts():
        text.set_color(INK)
    add_caption(fig, report)
    # Stack from the bottom: caption, legend, x label, then the panels.
    fig.canvas.draw()
    legend_top = legend.get_window_extent().transformed(fig.transFigure.inverted()).y1
    label_h = 0.3 / fig.get_figheight()
    fig.supxlabel("Wall time (s, log scale)", color=INK_SECONDARY, fontsize=10, y=legend_top + 0.005, va="bottom")
    fig.tight_layout(rect=(0.01, legend_top + label_h + 0.01, 1, 0.96))
    path = out_dir / "scatter.png"
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return path


def dataset_colors(datasets: list[str]) -> dict[str, str]:
    """ColorBrewer qualitative palettes: Dark2 (8 colors), then Paired (12)."""
    name = "Dark2" if len(datasets) <= 8 else "Paired"
    palette = colormaps[name].colors
    return {d: to_hex(palette[i % len(palette)]) for i, d in enumerate(datasets)}


# Per-task bars -----------------------------------------------------------


def version_colors(versions: list[str]) -> dict[str, str]:
    """ColorBrewer Blues, light (oldest) to dark (newest), skipping the
    lightest steps that disappear on the surface."""
    blues = colormaps["Blues"]
    n = len(versions)
    if n == 1:
        return {versions[0]: to_hex(blues(0.85))}
    return {v: to_hex(blues(0.4 + 0.55 * i / (n - 1))) for i, v in enumerate(versions)}


def plot_task(task, versions, colors, dataset, out_dir) -> Path:
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
    fig.suptitle(f"{task['label']}  [{dataset}]", x=0.01, ha="left", color=INK, fontsize=11)
    fig.tight_layout()
    path = out_dir / f"{slug(task['id'])}.png"
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return path


# Shared ------------------------------------------------------------------


def style_axes(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_SECONDARY, labelsize=9)
    ax.set_axisbelow(True)


def set_log_ticks(axis, values: list[float]) -> None:
    """Plain-number ticks on a log axis (0.5, 1, 2 instead of 5×10⁻¹).
    Narrow ranges get a few even ticks, wider ones 1-2-5 steps, and very
    wide ones only powers of ten."""
    values = [v for v in values if v > 0]
    if not values:
        return
    spread = max(values) / min(values)
    if spread < 3:
        axis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=2))
    else:
        axis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5) if spread < 100 else (1,)))
    axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    axis.set_minor_formatter(NullFormatter())


def add_caption(fig, report) -> None:
    machine = report["machine"]
    caption = (
        f"{machine.get('cpu', '')} · {machine.get('cores')} cores · {machine.get('os')} {machine.get('arch')} · "
        f"{report['config'].get('replicates')} replicates · scale {report['config'].get('scale')} · {report['created'][:10]}"
    )
    fig.text(0.01, 0.005, caption, fontsize=7, color=INK_SECONDARY, ha="left", va="bottom")


def slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "-", text)
