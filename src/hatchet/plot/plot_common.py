"""Shared plotting helpers used by multiple plot-* commands."""

import os
import logging

import matplotlib.pyplot as plt


def plot_bars(pdf, sids, values, title, ylabel, ylim, dpi=150, transparent=False):
    """Render one bar plot page (sample → value) into the open PdfPages."""
    fig_w = max(6, 0.4 * len(sids) + 2)
    fig, ax = plt.subplots(figsize=(fig_w, 4.5))
    ax.bar(range(len(sids)), values, color="#4c78a8", edgecolor="black", linewidth=0.4)
    pad = (ylim[1] - ylim[0]) * 0.01
    for x, v in enumerate(values):
        ax.text(x, v + pad, f"{v:.2f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(range(len(sids)))
    ax.set_xticklabels(sids, rotation=45, ha="right", fontsize=8)
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    pdf.savefig(fig, dpi=dpi, bbox_inches="tight", transparent=transparent)
    plt.close(fig)


def plot_summary_pdf(summary_rows, out_file, dpi=150, transparent=False):
    """Emit per-sample tumor-purity and tumor-ploidy bar plots, one page per
    (metric × cancer_type), into a single PDF.

    ``summary_rows`` is a list of ``(cancer_type, row_label, sample_id,
    purity, ploidy)`` tuples in panel-row order, keeping same-seg-file
    samples adjacent.
    """
    from matplotlib.backends.backend_pdf import PdfPages

    if not summary_rows:
        logging.warning("plot_summary: no samples to plot")
        return
    base, _ = os.path.splitext(out_file)
    bar_file = f"{base}_summary.pdf"

    groups = {}
    for ct, row_label, sid, purity, ploidy in summary_rows:
        groups.setdefault(ct, []).append((row_label, sid, purity, ploidy))

    metrics = [
        ("tumor purity", 2, (0, 1.05)),
        ("tumor ploidy", 3, None),
    ]
    with PdfPages(bar_file) as pdf:
        for metric_label, idx, fixed_ylim in metrics:
            for ct in sorted(groups):
                entries = groups[ct]
                sids = [
                    f"{rl}/{sid}" if rl and rl != sid else sid
                    for rl, sid, _, _ in entries
                ]
                values = [e[idx] for e in entries]
                ylim = fixed_ylim or (0, max(values + [1]) * 1.15)
                title = (
                    f"{ct} {metric_label} bar plot"
                    if ct
                    else f"{metric_label} bar plot"
                )
                plot_bars(
                    pdf,
                    sids,
                    values,
                    title,
                    metric_label,
                    ylim,
                    dpi=dpi,
                    transparent=transparent,
                )
    logging.info(f"summary bar plot saved to {bar_file}")
