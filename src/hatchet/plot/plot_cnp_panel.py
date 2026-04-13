"""Standalone script to plot a CNP panel over multiple samples."""

import os
import re
import logging
import argparse

import numpy as np
import pandas as pd

from hatchet.utils import *
from hatchet.hatchet_parser import add_arguments_plot_panel
from hatchet.plot.plot_cn_utils import *
from hatchet.plot.plot_utils import override_solution


def run(args=None):
    logging.info("run hatchet plot_panel")
    if isinstance(args, argparse.Namespace):
        args = vars(args)

    title = args["title"]
    panel_file = args["panel_file"]
    genome_size = args["genome_size"]
    region_bed = args["region_bed"]
    out_file = args["out_file"]
    plot_1d2d = args["plot_1d2d"]

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

        stats = []
        for sid in samples_in_file:
            sp = seg_info_all.loc[seg_info_all["SAMPLE"] == sid, :].reset_index(
                drop=True
            )
            cps = sp[[f"u_{c}" for c in clones]].iloc[0].tolist()
            purity = round(np.sum(cps[1:]), 2)
            ploidy = round(compute_tumor_ploidy(sp, clones, np.sum(cps[1:])), 2)
            stats.append((sid, purity, ploidy))
            logging.info(f"{label} / {sid}: purity={purity}, ploidy={ploidy}")

        clone_ploidies = compute_clone_ploidies(seg_info, clones)
        plot_cnv_profile(
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
    plot_cnv_legend(ax_leg)

    main_axes[0].set_title(title, pad=30)
    plt.savefig(out_file, dpi=dpi, bbox_inches="tight", transparent=transparent)
    plt.close()

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
                    "style": "cnv",
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


def _format_pool_label(tag):
    """Convert 'pool_p0.05_s1' to 'p=0.05,s=1'."""
    m = re.match(r"pool_p([^_]+)_s(\d+)", tag)
    if m:
        return f"p={m.group(1)},s={m.group(2)}"
    return tag


def plot_pool_cnp(
    pool_entries,
    region_bed,
    out_file,
    title=None,
    width=20,
    height=1,
    dpi=150,
    style="cnv",
):
    """Plot a multi-row CNP panel PDF, one row per Pareto-optimal pool solution.

    Args:
        pool_entries: List of ``(label, seg_df, imf_obj, is_pareto, is_selected)``
            tuples. Only entries where ``is_pareto`` is True are plotted.
        region_bed: Path to the whitelist region BED file.
        out_file: Output file path for the saved figure.
        title: Optional figure title placed above the top row.
        width: Figure width in inches.
        height: Height in inches per profile row (legend row is 2x this).
        dpi: Output resolution.
        style: ``"ascn"`` draws allele-specific A/B bars per clone;
            ``"cnv"`` draws total-CN colored bars.
    """
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"

    regions = read_region_bed(region_bed)

    valid = sorted(
        [
            (label, df, obj, selected)
            for label, df, obj, pareto, selected in pool_entries
        ],
        key=lambda x: x[2],  # sort by IMF-obj, low to high
    )
    if not valid:
        logging.warning(f"plot_pool_cnp: no Pareto solutions, skipping {out_file}")
        return

    nrows = len(valid)
    # Scale row height with number of clones and samples to avoid label overlap
    first_seg_df = valid[0][1]
    n_clones = len([c for c in first_seg_df.columns if c.startswith("cn_")])
    n_samples = first_seg_df["SAMPLE"].nunique()
    # Base height scales with clones; add extra per sample for ylabel text
    row_h = height * max(1, n_clones - 1) + 0.2 * max(0, n_samples - 1)
    fig, axes = plt.subplots(
        nrows=nrows + 1,
        ncols=1,
        figsize=(width, row_h * nrows),
        gridspec_kw={"height_ratios": [row_h] * nrows + [2 * height]},
    )
    fig.subplots_adjust(hspace=0.6 + 0.1 * max(0, n_samples - 1))
    main_axes = axes[:-1]
    ax_leg = axes[-1]

    for i, (label, seg_df, obj, is_selected) in enumerate(valid):
        seg_info_all, clones, _ = prepare_seg_ucn(seg_df)
        samples = seg_info_all["SAMPLE"].unique().tolist()
        sample_first = samples[0]
        seg_info = seg_info_all.loc[
            seg_info_all["SAMPLE"] == sample_first, :
        ].reset_index(drop=True)

        # Per-sample purity/ploidy
        sample_stats = []
        for sid in samples:
            sp = seg_info_all.loc[seg_info_all["SAMPLE"] == sid, :].reset_index(
                drop=True
            )
            cps = sp[[f"u_{c}" for c in clones]].iloc[0].tolist()
            purity = round(np.sum(cps[1:]), 2)
            ploidy = round(compute_tumor_ploidy(sp, clones, np.sum(cps[1:])), 2)
            sample_stats.append((sid, purity, ploidy))

        clone_ploidies = compute_clone_ploidies(seg_info, clones)
        profile_fn = plot_ascn_profile if style == "ascn" else plot_cnv_profile
        profile_fn(
            main_axes[i],
            seg_info,
            regions,
            plot_chrname=True,
            width=width,
            height=height,
            show_clone_name=False,
            show_prop=True,
            clone_ploidies=clone_ploidies,
        )
        short_label = _format_pool_label(label)
        if is_selected:
            short_label += " *"
        stats_lines = "\n".join(f"{sid}: p={p} pl={pl}" for sid, p, pl in sample_stats)
        ylabel = f"{short_label}\nimf {round(obj, 2)}\n{stats_lines}"
        color = "red" if is_selected else "black"
        main_axes[i].set_ylabel(
            ylabel, rotation=0, ha="right", va="center", color=color
        )

    legend_fn = plot_ascn_legend if style == "ascn" else plot_cnv_legend
    legend_fn(ax_leg)

    if title:
        main_axes[0].set_title(title)
    plt.savefig(out_file, dpi=dpi, bbox_inches="tight")
    plt.close()
    logging.info(f"pool CNP panel saved to {out_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="HATCHet plot_panel",
        description="plot panel",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    add_arguments_plot_panel(parser)
    args = parser.parse_args()
    setup_logging(args)
    run(args)
