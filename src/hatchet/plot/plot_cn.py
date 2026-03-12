import os
import sys
import time
import logging
import argparse
import numpy as np
import pandas as pd

from hatchet.utils import *
from hatchet.hatchet_parser import add_arguments_plot_cn
from hatchet.plot.plot_utils import *
from hatchet.plot.plot_1d2d import *
from hatchet.plot.plot_cn_utils import *


def run(args=None):
    logging.info("run hatchet plot-cn one sample")
    if isinstance(args, argparse.Namespace):
        args = vars(args)

    ##################################################
    # files
    bbc_ucn = args["bbc"]
    seg_ucn = args["seg"]
    genome_size = args["genome_size"]
    region_bed = args["region_bed"]

    solfile = args["solfile"]
    gamma_file = args["gamma_file"]
    # is_diploid = not args["tetraploid"]
    plot_dir = args["plot_dir"]
    os.makedirs(plot_dir, exist_ok=True)
    add_file_logging(plot_dir, "plot-cn")

    ##################################################
    # parameters
    row_width = 20
    row_height = 6
    dpi = args["dpi"]

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"

    ignore_gap = not args["keep_gap"]
    dpi = args["dpi"]
    transparent = args["transparent"]
    file_type = args["img_type"]

    # TODO per-cluster transparency
    tail_alpha = args["tail_alpha"]
    nontail_alpha = args["center_alpha"]
    onetail_area = args["onetail_area"]

    # figure axis limits
    maxlim_fcn = args["maxlim_fcn"]

    tol = 1e-2
    get_title = lambda sid: f"Sample: {sid}"
    get_filename = lambda sid: str(sid) + (f".{solID}" if solID != "" else "")
    ##################################################
    # load files
    segs, clones, clone_props = read_seg_ucn_file(seg_ucn)
    n_clones = len(clones)
    n_tumors = n_clones - 1
    bbcs = read_seg_ucn_file(bbc_ucn)[0]
    chrom_sizes = read_genome_sizes(genome_size)
    regions = read_region_bed(region_bed)

    ##################################################
    chrs = segs["#CHR"].unique().tolist()
    assert len(chrs) > 0, "No chromosomes are found"
    assert all(str.startswith(ch, "chr") for ch in chrs), (
        "BBC #CHR must start with chr-prefix"
    )
    samples = segs["SAMPLE"].unique().tolist()
    clusters = sorted(bbcs["CLUSTER"].unique().tolist())

    solID = ""
    if solfile != None:
        bbcs, segs, n_clones, n_tumors, solID = override_solution(
            bbcs, segs, samples, clusters, n_clones, solfile
        )
    assert n_tumors > 0, "at least one tumor clone must present"

    clone_states = bbcs["CNP"].unique().tolist()
    is_diploid = args["ploidy"] == "diploid"
    gammas = load_gammas(gamma_file, is_diploid)
    logging.info(f"Gammas: {gammas}")
    assert len(gammas) == len(samples), "gamma-file doesn't match samples from BBC file"

    ##################################################
    # start plotting 1D and 2D
    state_style, _ = get_cn_colors()
    sns.set_style("whitegrid")
    for sample in samples:
        logging.info(f"plot {sample}")
        outfile_1d = os.path.join(plot_dir, f"{get_filename(sample)}.1D.{file_type}")
        outfile_2d = os.path.join(plot_dir, f"{get_filename(sample)}.2D.{file_type}")
        bin_info: pd.DataFrame = bbcs[bbcs["SAMPLE"] == sample].reset_index(drop=True)
        seg_info: pd.DataFrame = segs[segs["SAMPLE"] == sample].reset_index(drop=True)
        # CNP profile
        cnp_ids = bin_info["CNP"].to_numpy()
        clone_states = bin_info["CNP"].unique().tolist()
        clone_props = bin_info.iloc[0][[f"u_{clone}" for clone in clones]].to_numpy()
        palette = make_cnp_palette(clone_states, clone_props, state_style)

        # compute per-segment FCN
        bin_info["FCN"] = bin_info.apply(
            func=lambda row: row["RD"] * gammas[sample], axis=1
        )
        tumor_purity = round(np.sum(clone_props[1:]), 3)
        tumor_ploidy = round(
            compute_tumor_ploidy(seg_info, clones, np.sum(clone_props[1:])), 3
        )
        logging.info(f"purity={tumor_purity}, ploidy={tumor_ploidy}")

        # compute expected FCNs over clustres
        bin_info["exp-BAF"] = 0.0
        bin_info["exp-FCN"] = 0.0
        exp_bafs = np.zeros(len(clone_states), dtype=np.float32)
        exp_fcns = np.zeros(len(clone_states), dtype=np.float32)
        exp_labels = []
        for i, clone_state in enumerate(clone_states):
            states = [
                (int(x.split("|")[0]), int(x.split("|")[1]))
                for x in clone_state.split(";")
            ]
            exp_labels.append(states)
            _, _, exp_fcns[i], exp_bafs[i] = get_expected_baf_fcn(states, clone_props)
            bin_info.loc[bin_info["CNP"] == clone_state, "exp-BAF"] = exp_bafs[i]
            bin_info.loc[bin_info["CNP"] == clone_state, "exp-FCN"] = exp_fcns[i]

        lim_baf = (0, 1) if bin_info["BAF"].max() > 0.5 else (0, 0.55)
        max_fcn = np.round(bin_info["FCN"].max()).astype(int)
        if max_fcn > maxlim_fcn:
            num_exceeded = np.sum(bin_info["FCN"] >= maxlim_fcn)
            logging.warning(
                f"there are {num_exceeded} bins having FCN exceed maxlim_fcn={maxlim_fcn}"
            )
        lim_fcn = (0, min(max(3, max_fcn), maxlim_fcn))

        _, g0_colors = plot_2d(
            sample,
            bin_info,
            bin_info["BAF"].to_numpy(),
            bin_info["FCN"].to_numpy(),
            exp_bafs,
            exp_fcns,
            exp_labels,
            clone_props,
            alphas=None,
            hue=cnp_ids,
            palette=palette,
            label_clone=True,
            xlab="Minor haplotype B-allele frequency (mhBAF)",
            ylab="Fractional copy number (FCN)",
            xlim=lim_baf,
            ylim=lim_fcn,
            title=f"sample={sample}; purity={tumor_purity}; ploidy={tumor_ploidy}",
            dpi=dpi,
            transparent=transparent,
            out_file=outfile_2d,
        )
        fig, axes = plt.subplots(
            nrows=4,
            ncols=1,
            figsize=(row_width, row_height),
            gridspec_kw={"height_ratios": [3, 3, 2, 1]},
        )
        main_axes = axes[:-1]
        ax_leg = axes[-1]
        for i, [ctype, ylim] in enumerate([["FCN", lim_fcn], ["BAF", lim_baf]]):
            plot_1d(
                main_axes[i],
                sample,
                bin_info,
                bin_info[ctype],
                regions,
                chrom_sizes,
                exp_colname=f"exp-{ctype}",
                val_type=ctype,
                colors=g0_colors,
                hue=cnp_ids if i == 0 else None,
                palette=palette if i == 0 else None,
                ylim=ylim,
                ylab=ctype,
                plot_chrname=True,
                ignore_gap=ignore_gap,
                show_legend=False,
            )

        plot_cnv_profile(
            main_axes[2],
            seg_info,
            regions,
            width=row_width,
            height=1,
            plot_chrname=False,
            show_clone_name=True,
            show_prop=True,
        )
        plot_cnv_legend(ax_leg)

        fig.suptitle(f"sample={sample}; purity={tumor_purity}; ploidy={tumor_ploidy}")
        axes[1].set_ylabel("mhBAF")

        axes[0].grid(False)
        axes[1].grid(False)
        # axes[0].legend(markerscale=6)
        # sns.move_legend(axes[0], "upper left", bbox_to_anchor=(1, 1), title=None)
        plt.tight_layout()
        plt.savefig(outfile_1d, dpi=dpi, bbox_inches="tight", transparent=transparent)
        plt.close(fig)
        logging.info(f"finish {sample}")
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="HATCHet plot 1D2D",
        description="plot HATCHet results",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    add_arguments_plot_cn(parser)
    args = parser.parse_args()
    setup_logging(args)
    run(args)
