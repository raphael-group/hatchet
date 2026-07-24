"""Standalone script to plot a CNP panel over multiple samples."""

import os
import logging
from itertools import groupby

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from hatchet.utils import (
    compute_clone_ploidies,
    compute_tumor_ploidy,
    normalize_args,
    setup_logging,
)
from hatchet.io_utils import override_solution, read_region_bed, read_seg_ucn_file
from cnplot import plot_cnv_profile
from hatchet.plot.plot_utils import build_genome_axis, use_editable_fonts


def run(args=None):
    args = normalize_args(args)
    setup_logging(args)
    logging.info("run hatchet plot_panel")

    title = args["title"]
    panel_file = args["panel_file"]
    genome_size = args["genome_size"]
    region_bed = args["region_bed"]
    out_file = args["out_file"]
    plot_1d2d = args["plot_1d2d"]
    plot_summary = args["plot_summary"]

    row_width = args["width"]
    row_height = args["height"]
    show_clone_name = args["show_clone_name"]
    show_prop = args["show_prop"]
    show_ploidy = args["show_ploidy"]
    min_prop = args["min_prop"]
    dpi = args["dpi"]
    transparent = args["transparent"]

    use_editable_fonts()

    regions = read_region_bed(region_bed)

    panel = pd.read_table(panel_file, sep="\t", index_col=False).fillna("")
    if plot_1d2d:
        assert "PATH_TO_BBC" in panel.columns, (
            "--plot_1d2d requires PATH_TO_BBC column in panel TSV"
        )

    # Restrict the axis to chromosomes present in the panel's seg files, so an
    # autosome-only run does not draw empty chrX/chrY.
    panel_chroms = set()
    for seg_path in panel["PATH_TO_SEG"]:
        panel_chroms.update(
            pd.read_table(seg_path, sep="\t", usecols=["#CHR"])["#CHR"].unique()
        )
    genome_axis = build_genome_axis(region_bed, genome_size, keep_chroms=panel_chroms)
    nrows = len(panel)
    fig, axes = plt.subplots(
        nrows=nrows + 1,
        ncols=1,
        figsize=(row_width, row_height * nrows),
        gridspec_kw={"height_ratios": [row_height] * nrows + [2 * row_height]},
    )
    main_axes = axes[:-1]
    ax_leg = axes[-1]

    summary_rows = []
    for i, row in panel.iterrows():
        label = row["SAMPLE"]
        seg_ucn = row["PATH_TO_SEG"]
        seg_info_all, clones = read_seg_ucn_file(seg_ucn)

        solfile = row.get("PATH_TO_SOLFILE", "")
        bbc_path = row.get("PATH_TO_BBC", "")
        if solfile and bbc_path:
            bbcs = read_seg_ucn_file(bbc_path)[0]
            samples_all = sorted(bbcs["SAMPLE"].unique().tolist())
            clusters = sorted(bbcs["CLUSTER"].unique().tolist())
            _, seg_info_all, n_clones, _, _ = override_solution(
                bbcs, samples_all, clusters, len(clones), solfile, regions
            )
            clones = ["normal"] + [f"clone{i}" for i in range(1, n_clones)]

        samples_in_file = seg_info_all["SAMPLE"].unique().tolist()
        seg_info = seg_info_all.loc[
            seg_info_all["SAMPLE"] == samples_in_file[0], :
        ].reset_index(drop=True)

        cancer_type = str(row.get("cancer_type", "")).strip()
        for sid in samples_in_file:
            sp = seg_info_all.loc[seg_info_all["SAMPLE"] == sid, :].reset_index(
                drop=True
            )
            cps = sp[[f"u_{c}" for c in clones]].iloc[0].tolist()
            purity = round(np.sum(cps[1:]), 2)
            ploidy = round(compute_tumor_ploidy(sp, clones, np.sum(cps[1:])), 2)
            summary_rows.append((cancer_type, label, sid, purity, ploidy))
            logging.info(f"{label} / {sid}: purity={purity}, ploidy={ploidy}")

        # Display clone names passed to cnplot via `clones=`: "Clone i" / "i", with
        # the per-sample proportion inlined as "i (xx.xx%)" when show_prop so the
        # label stays one line and aligns with its y-tick. cnplot reads the cn_/u_
        # columns by that name, so rename them (and the ploidy keys) to match.
        props0 = seg_info[[f"u_{c}" for c in clones]].iloc[0]
        disp = {}
        for c in clones[1:]:
            if props0[f"u_{c}"] < min_prop:
                continue  # hide clones below the display threshold
            name = f"Clone {c[len('clone'):]}" if show_clone_name else c[len("clone") :]
            if show_prop:
                name = f"{name} ({props0[f'u_{c}'] * 100:.2f}%)"
            disp[c] = name
        seg_disp = seg_info.rename(
            columns={f"cn_{c}": f"cn_{d}" for c, d in disp.items()}
            | {f"u_{c}": f"u_{d}" for c, d in disp.items()}
        )
        clone_ploidies = None
        if show_ploidy:
            per_clone = compute_clone_ploidies(seg_info, clones)
            clone_ploidies = {disp.get(k, k): v for k, v in per_clone.items()}
        plot_cnv_profile(
            main_axes[i],
            seg_disp,
            genome_axis,
            ax_leg=(ax_leg if i == nrows - 1 else None),
            plot_chrname=(i == 0),
            clones=list(disp.values()),
            show_prop=False,  # proportion is inlined into the clone name
            clone_ploidies=clone_ploidies,
        )
        main_axes[i].set_ylabel(label, rotation=0, ha="right", va="center")

    main_axes[0].set_title(title, pad=30)
    plt.savefig(out_file, dpi=dpi, bbox_inches="tight", transparent=transparent)
    plt.close()

    if plot_summary:
        plot_summary_pdf(summary_rows, out_file, dpi=dpi, transparent=transparent)

    if plot_1d2d:
        from hatchet.plot.plot_cn import run as run_plot_cn

        base_dir = os.path.dirname(os.path.abspath(out_file)) or "."
        for _, row in panel.iterrows():
            bbc_path = row["PATH_TO_BBC"]
            seg_path = row["PATH_TO_SEG"]
            solfile = row.get("PATH_TO_SOLFILE", "") or None
            gamma_file = os.path.join(os.path.dirname(bbc_path), "gammas.tsv")
            bn = os.path.basename(bbc_path)
            if "tetraploid" in bn:
                ploidy = "tetraploid"
            elif "diploid" in bn:
                ploidy = "diploid"
            else:
                raise ValueError(f"cannot infer ploidy from bbc filename: {bn}")
            label = str(row["SAMPLE"]).replace("/", "_").replace(" ", "_")
            plot_dir = os.path.join(base_dir, f"{label}_1d2d")
            logging.info(f"plot_1d2d: {label} → {plot_dir}")
            run_plot_cn(
                {
                    "bbc": bbc_path,
                    "seg": seg_path,
                    "solfile": solfile,
                    "gamma_file": gamma_file,
                    "genome_size": genome_size,
                    "region_bed": region_bed,
                    "plot_dir": plot_dir,
                    "ploidy": ploidy,
                    "dpi": dpi,
                    "img_type": "png",
                    "transparent": transparent,
                    "show_gap": False,
                    "tail_alpha": 0.8,
                    "center_alpha": 1.0,
                    "onetail_area": 0.025,
                    "maxlim_fcn": 30,
                }
            )

    logging.info("Done")
    return


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
