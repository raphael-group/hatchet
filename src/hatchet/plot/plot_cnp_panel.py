import os
import sys
import time
import logging
import argparse
import numpy as np
import pandas as pd

from hatchet.utils import *
from hatchet.hatchet_parser import add_arguments_plot_panel
from hatchet.plot.plot_cn_utils import *

"""
Standalone script, plot CNP panel over multiple samples    
"""


def run(args=None):
    logging.info("run hatchet plot_panel")
    if isinstance(args, argparse.Namespace):
        args = vars(args)

    title = args["title"]
    panel_file = args["panel_file"]
    genome_size = args["genome_size"]
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
        ploidy = row["PLOIDY"]
        seg_ucn = row["PATH_TO_SEG"]
        seg_info, clones, clone_props = read_seg_ucn_file(seg_ucn)

        # all replicates shared same CNP
        dummy_sample = seg_info["SAMPLE"].iloc[0]
        seg_info = seg_info.loc[seg_info["SAMPLE"] == dummy_sample, :].reset_index(
            drop=True
        )

        # assumes single-sample TODO
        tumor_purity = round(np.sum(clone_props[1:]), 2)
        tumor_ploidy = round(
            compute_tumor_ploidy(seg_info, clones, np.sum(clone_props[1:])), 2
        )
        logging.info(f"sample={sample}, purity={tumor_purity}, ploidy={tumor_ploidy}")

        plot_cnv_profile(
            main_axes[i],
            seg_info,
            regions,
            plot_chrname=i == 0,
            width=row_width,
            height=row_height,
            show_clone_name=show_clone_name,
            show_prop=show_prop,
        )
        ylabel = f"{sample}\npurity {tumor_purity}\nploidy {tumor_ploidy}"
        main_axes[i].set_ylabel(ylabel, rotation=0, ha="right", va="center")
    plot_cnv_legend(ax_leg)

    main_axes[0].set_title(title)
    plt.savefig(out_file, dpi=dpi, bbox_inches="tight", transparent=transparent)
    plt.close()

    logging.info("Done")
    return


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
