"""Standalone script to plot a CNP panel over multiple samples."""

import re
import logging
import argparse

import numpy as np
import pandas as pd

from hatchet.utils import *
from hatchet.hatchet_parser import add_arguments_plot_panel
from hatchet.plot.plot_cn_utils import *


def run(args=None):
    logging.info("run hatchet plot_panel")
    if isinstance(args, argparse.Namespace):
        args = vars(args)

    title = args["title"]
    panel_file = args["panel_file"]
    region_bed = args["region_bed"]
    out_file = args["out_file"]

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
        sample = row["SAMPLE"]
        seg_ucn = row["PATH_TO_SEG"]
        seg_info, clones, clone_props = read_seg_ucn_file(seg_ucn)

        dummy_sample = seg_info["SAMPLE"].iloc[0]
        seg_info = seg_info.loc[seg_info["SAMPLE"] == dummy_sample, :].reset_index(
            drop=True
        )

        tumor_purity = round(np.sum(clone_props[1:]), 2)
        tumor_ploidy = round(
            compute_tumor_ploidy(seg_info, clones, np.sum(clone_props[1:])), 2
        )
        logging.info(f"sample={sample}, purity={tumor_purity}, ploidy={tumor_ploidy}")

        clone_ploidies = compute_clone_ploidies(seg_info, clones)
        plot_cnv_profile(
            main_axes[i],
            seg_info,
            regions,
            plot_chrname=True,
            width=row_width,
            height=row_height,
            show_clone_name=show_clone_name,
            show_prop=show_prop,
            clone_ploidies=clone_ploidies,
        )
        ylabel = f"{sample}\npurity {tumor_purity}\nploidy {tumor_ploidy}"
        main_axes[i].set_ylabel(ylabel, rotation=0, ha="right", va="center")
    plot_cnv_legend(ax_leg)

    main_axes[0].set_title(title)
    plt.savefig(out_file, dpi=dpi, bbox_inches="tight", transparent=transparent)
    plt.close()

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
            if pareto
        ],
        key=lambda x: x[2],  # sort by IMF-obj, low to high
    )
    if not valid:
        logging.warning(f"plot_pool_cnp: no Pareto solutions, skipping {out_file}")
        return

    nrows = len(valid)
    # Scale row height with number of clones to avoid label overlap
    first_seg_df = valid[0][1]
    n_clones = len([c for c in first_seg_df.columns if c.startswith("cn_")])
    row_h = height * max(1, n_clones - 1)
    fig, axes = plt.subplots(
        nrows=nrows + 1,
        ncols=1,
        figsize=(width, row_h * nrows),
        gridspec_kw={"height_ratios": [row_h] * nrows + [2 * height]},
    )
    fig.subplots_adjust(hspace=0.6)
    main_axes = axes[:-1]
    ax_leg = axes[-1]

    for i, (label, seg_df, obj, is_selected) in enumerate(valid):
        seg_info, clones, clone_props = prepare_seg_ucn(seg_df)
        dummy_sample = seg_info["SAMPLE"].iloc[0]
        seg_info = seg_info.loc[seg_info["SAMPLE"] == dummy_sample, :].reset_index(
            drop=True
        )

        tumor_purity = round(np.sum(clone_props[1:]), 2)
        tumor_ploidy = round(
            compute_tumor_ploidy(seg_info, clones, np.sum(clone_props[1:])), 2
        )

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
        ylabel = (
            f"{short_label}\nimf {round(obj, 2)}"
            f"\npurity {tumor_purity}\nploidy {tumor_ploidy}"
        )
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
