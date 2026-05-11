"""Standalone script to plot a CNP panel over multiple samples."""

import os
import logging

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from hatchet.utils import (
    compute_clone_ploidies,
    compute_tumor_ploidy,
    normalize_args,
    read_region_bed,
    read_seg_ucn_file,
    setup_logging,
)
from hatchet.plot.plot_cn_utils import (
    plot_ascn_legend,
    plot_ascn_profile,
    plot_cnv_legend,
    plot_cnv_profile,
)
from hatchet.plot.plot_common import plot_summary_pdf
from hatchet.plot.plot_utils import override_solution


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
    dpi = args["dpi"]
    transparent = args["transparent"]

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"

    regions = read_region_bed(region_bed)

    panel = pd.read_table(panel_file, sep="\t", index_col=False).fillna("")
    if plot_1d2d:
        assert "PATH_TO_BBC" in panel.columns, (
            "--plot_1d2d requires PATH_TO_BBC column in panel TSV"
        )
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
        seg_info_all, clones, _ = read_seg_ucn_file(seg_ucn)

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
        stats = []
        for sid in samples_in_file:
            sp = seg_info_all.loc[seg_info_all["SAMPLE"] == sid, :].reset_index(
                drop=True
            )
            cps = sp[[f"u_{c}" for c in clones]].iloc[0].tolist()
            purity = round(np.sum(cps[1:]), 2)
            ploidy = round(compute_tumor_ploidy(sp, clones, np.sum(cps[1:])), 2)
            stats.append((sid, purity, ploidy))
            summary_rows.append((cancer_type, label, sid, purity, ploidy))
            logging.info(f"{label} / {sid}: purity={purity}, ploidy={ploidy}")

        clone_ploidies = compute_clone_ploidies(seg_info, clones)
        _profile_fn = plot_ascn_profile if args["plot_ascn"] else plot_cnv_profile
        _profile_fn(
            main_axes[i],
            seg_info,
            regions,
            plot_chrname=(i == 0),
            width=row_width,
            height=row_height,
            show_clone_name=show_clone_name,
            show_prop=show_prop,
            clone_ploidies=clone_ploidies,
        )
        stats_lines = "\n".join(f"{sid}: p={p} pl={pl}" for sid, p, pl in stats)
        ylabel = f"{label}\n{stats_lines}"
        main_axes[i].set_ylabel(ylabel, rotation=0, ha="right", va="center")
    _legend_fn = plot_ascn_legend if args["plot_ascn"] else plot_cnv_legend
    _legend_fn(ax_leg)

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
                    "keep_gap": False,
                    "tail_alpha": 0.8,
                    "center_alpha": 1.0,
                    "onetail_area": 0.025,
                    "maxlim_fcn": 30,
                }
            )

    logging.info("Done")
    return
