import os
import logging
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from hatchet.utils import (
    add_file_logging,
    compute_clone_ploidies,
    compute_tumor_ploidy,
    normalize_args,
    read_genome_sizes,
    read_region_bed,
    read_seg_ucn_file,
    setup_logging,
)
from hatchet.plot.plot_utils import (
    get_expected_baf_fcn,
    load_gammas,
    override_solution,
    set_palette,
)
from hatchet.plot.plot_1d2d import get_transparency, plot_1d, plot_2d
from hatchet.plot.plot_cn_utils import (
    plot_ascn_legend,
    plot_ascn_profile,
    plot_cnv_legend,
    plot_cnv_profile,
)


def run(args=None):
    args = normalize_args(args)
    setup_logging(args)
    logging.info("run hatchet plot-cn one sample")

    ##################################################
    # files
    bbc_ucn = args["bbc"]
    seg_ucn = args["seg"]
    genome_size = args["genome_size"]
    region_bed = args["region_bed"]

    solfile = args["solfile"]
    gamma_file = args["gamma_file"]
    plot_dir = args["plot_dir"]
    os.makedirs(plot_dir, exist_ok=True)
    add_file_logging(plot_dir, "plot-cn")

    ##################################################
    # parameters
    row_width = 20
    row_height = 6

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"

    ignore_gap = not args["keep_gap"]
    dpi = args["dpi"]
    transparent = args["transparent"]
    file_type = args["img_type"]

    tail_alpha = args["tail_alpha"]
    center_alpha = args["center_alpha"]
    onetail_area = args["onetail_area"]

    # figure axis limits
    maxlim_fcn = args["maxlim_fcn"]

    def get_filename(sample_id):
        suffix = f".{solID}" if solID != "" else ""
        return str(sample_id) + suffix

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
    if solfile is not None:
        bbcs, segs, n_clones, n_tumors, solID = override_solution(
            bbcs, samples, clusters, n_clones, solfile, regions
        )
    assert n_tumors > 0, "at least one tumor clone must present"

    clone_states = bbcs["CNP"].unique().tolist()
    is_diploid = args["ploidy"] == "diploid"
    gammas = load_gammas(gamma_file, is_diploid)
    logging.info(f"Gammas: {gammas}")
    assert len(gammas) == len(samples), "gamma-file doesn't match samples from BBC file"

    ##################################################
    # start plotting 1D and 2D
    from matplotlib.backends.backend_pdf import PdfPages

    sns.set_style("whitegrid")
    for sample in samples:
        logging.info(f"plot {sample}")
        outfile = os.path.join(plot_dir, f"{get_filename(sample)}.{file_type}")
        bin_info: pd.DataFrame = bbcs[bbcs["SAMPLE"] == sample].reset_index(drop=True)
        seg_info: pd.DataFrame = segs[segs["SAMPLE"] == sample].reset_index(drop=True)
        # CNP profile
        cnp_ids = bin_info["CNP"].to_numpy()
        clone_states = bin_info["CNP"].unique().tolist()
        clone_props = bin_info.iloc[0][[f"u_{clone}" for clone in clones]].to_numpy()
        palette = dict(zip(clone_states, set_palette(num_colors=len(clone_states))))

        # compute per-segment FCN
        bin_info["FCN"] = bin_info.apply(
            func=lambda row: row["RD"] * gammas[sample], axis=1
        )
        tumor_purity = round(np.sum(clone_props[1:]), 3)
        tumor_ploidy = round(
            compute_tumor_ploidy(seg_info, clones, np.sum(clone_props[1:])), 3
        )
        logging.info(f"purity={tumor_purity}, ploidy={tumor_ploidy}")

        # compute expected FCNs over clusters
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
        max_fcn = int(np.ceil(bin_info["FCN"].max()))
        if max_fcn > maxlim_fcn:
            num_exceeded = np.sum(bin_info["FCN"] >= maxlim_fcn)
            logging.warning(
                f"there are {num_exceeded} bins having FCN exceed maxlim_fcn={maxlim_fcn}"
            )
        lim_fcn = (0, min(max(3, max_fcn), maxlim_fcn))

        sample_title = f"sample={sample}; purity={tumor_purity}; ploidy={tumor_ploidy}"

        alphas = get_transparency(
            bin_info,
            by="CNP",
            one_tail=onetail_area,
            tail_alpha=tail_alpha,
            nontail_alpha=center_alpha,
        ).to_numpy()

        # Page 1: 2D scatter
        fig_2d, g0_colors = plot_2d(
            sample,
            bin_info,
            bin_info["BAF"].to_numpy(),
            bin_info["FCN"].to_numpy(),
            exp_bafs,
            exp_fcns,
            exp_labels,
            clone_props,
            alphas=alphas,
            hue=cnp_ids,
            palette=palette,
            label_clone=True,
            xlab="Minor haplotype B-allele frequency (mhBAF)",
            ylab="Fractional copy number (FCN)",
            xlim=lim_baf,
            ylim=lim_fcn,
            title=sample_title,
            dpi=dpi,
            transparent=transparent,
        )

        # Page 2: 1D scatter + CNP profile
        cnp_h = max(3, n_tumors * 2)  # scale with clones, same as pool panel
        fig_1d, axes = plt.subplots(
            nrows=4,
            ncols=1,
            figsize=(row_width, row_height + cnp_h * 0.5),
            gridspec_kw={"height_ratios": [3, 3, cnp_h, 1]},
        )
        for i, [ctype, ylim] in enumerate([["FCN", lim_fcn], ["BAF", lim_baf]]):
            plot_1d(
                axes[i],
                sample,
                bin_info,
                bin_info[ctype],
                regions,
                chrom_sizes,
                exp_colname=f"exp-{ctype}",
                exp_groups=cnp_ids,
                val_type=ctype,
                colors=g0_colors,
                hue=cnp_ids if i == 0 else None,
                palette=palette if i == 0 else None,
                alphas=alphas,
                ylim=ylim,
                ylab=ctype,
                plot_chrname=True,
                ignore_gap=ignore_gap,
                show_legend=False,
            )

        clone_ploidies = compute_clone_ploidies(seg_info, clones)
        _profile_fn = plot_ascn_profile if args["plot_ascn"] else plot_cnv_profile
        _profile_fn(
            axes[2],
            seg_info,
            regions,
            width=row_width,
            height=1,
            plot_chrname=True,
            show_clone_name=True,
            show_prop=True,
            clone_ploidies=clone_ploidies,
        )
        _legend_fn = plot_ascn_legend if args["plot_ascn"] else plot_cnv_legend
        _legend_fn(axes[-1])

        fig_1d.suptitle(sample_title)
        axes[1].set_ylabel("mhBAF")
        axes[0].grid(False)
        axes[1].grid(False)
        plt.tight_layout()

        # Save: single PDF (page1=2D, page2=1D+CNP) or separate files
        if file_type == "pdf":
            with PdfPages(outfile) as pdf:
                pdf.savefig(
                    fig_2d, dpi=dpi, bbox_inches="tight", transparent=transparent
                )
                pdf.savefig(
                    fig_1d, dpi=dpi, bbox_inches="tight", transparent=transparent
                )
        else:
            fig_2d.savefig(
                outfile.replace(f".{file_type}", f".2D.{file_type}"),
                dpi=dpi,
                bbox_inches="tight",
                transparent=transparent,
            )
            fig_1d.savefig(
                outfile.replace(f".{file_type}", f".1D.{file_type}"),
                dpi=dpi,
                bbox_inches="tight",
                transparent=transparent,
            )
        plt.close(fig_2d)
        plt.close(fig_1d)
        logging.info(f"finish {sample}")
    return
