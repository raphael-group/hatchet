"""Shared plotting helpers used by multiple plot-* commands."""

import os
import logging
from itertools import groupby

import matplotlib.pyplot as plt


def plot_bars(pdf, groups, title, ylabel, ylim, dpi=150, transparent=False):
    """Render one grouped bar plot page into an open PdfPages.

    ``groups`` is a list of ``(group_label, [value, ...])`` tuples. Each
    group draws adjacent bars; a single xtick sits at the group's center.
    Small gaps separate groups.
    """
    bar_step = 0.3
    bar_width = 0.25
    group_gap = 0.4
    positions, flat_values, centers, labels = [], [], [], []
    cursor = 0.0
    for label, values in groups:
        start = cursor
        for v in values:
            positions.append(cursor)
            flat_values.append(v)
            cursor += bar_step
        end = cursor - bar_step
        centers.append((start + end) / 2)
        labels.append(label)
        cursor += group_gap

    fig_w = max(6, 0.4 * len(flat_values) + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, 4.5))
    ax.bar(
        positions,
        flat_values,
        width=bar_width,
        color="#4c78a8",
        edgecolor="black",
        linewidth=0.4,
    )
    pad = (ylim[1] - ylim[0]) * 0.01
    for x, v in zip(positions, flat_values):
        ax.text(
            x,
            v + pad,
            f"{v:.2f}",
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=90,
        )
    ax.set_xticks(centers)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_xlim(-0.3, cursor - group_gap + 0.3 if cursor > 0 else 1)
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel)
    ax.set_title(title, pad=15)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    pdf.savefig(fig, dpi=dpi, bbox_inches="tight", transparent=transparent)
    plt.close(fig)


def plot_summary_pdf(summary_rows, out_file, dpi=150, transparent=False):
    """Emit per-sample tumor-purity and tumor-ploidy bar plots, one page per
    (metric, cancer_type), into a single PDF.

    ``summary_rows`` is a list of ``(cancer_type, row_label, sample_id,
    purity, ploidy)`` tuples in panel-row order. Consecutive entries
    sharing ``row_label`` render as one group under a shared xtick; order
    is preserved so same-seg-file samples stay adjacent.
    """
    from matplotlib.backends.backend_pdf import PdfPages

    if not summary_rows:
        logging.warning("plot_summary: no samples to plot")
        return
    base, _ = os.path.splitext(out_file)
    bar_file = f"{base}_summary.pdf"

    ct_groups = {}
    for ct, row_label, sid, purity, ploidy in summary_rows:
        ct_groups.setdefault(ct, []).append((row_label, sid, purity, ploidy))

    metrics = [
        ("tumor purity", 2, (0, 1.15)),
        ("tumor ploidy", 3, None),
    ]
    with PdfPages(bar_file) as pdf:
        for metric_label, idx, fixed_ylim in metrics:
            for ct in sorted(ct_groups):
                entries = ct_groups[ct]
                groups = [
                    (label, [e[idx] for e in run])
                    for label, run in groupby(entries, key=lambda e: e[0])
                ]
                all_values = [v for _, vs in groups for v in vs]
                ylim = fixed_ylim or (0, max(all_values + [1]) * 1.15)
                title = (
                    f"{ct} {metric_label} bar plot"
                    if ct
                    else f"{metric_label} bar plot"
                )
                plot_bars(
                    pdf,
                    groups,
                    title,
                    metric_label,
                    ylim,
                    dpi=dpi,
                    transparent=transparent,
                )
    logging.info(f"summary bar plot saved to {bar_file}")
