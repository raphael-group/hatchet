import os
import logging
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from hatchet.utils import (
    add_file_logging,
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
    # parameters (styling — defaults in hatchet.yaml)
    row_width = args["plot_row_width"]
    row_height = args["plot_row_height"]
    inter_sample_hspace = args["plot_inter_sample_hspace"]
    intra_sample_hspace = args["plot_intra_sample_hspace"]
    baf_cnp_hspace = args["plot_baf_cnp_hspace"]
    cnp_legend_hspace = args["plot_cnp_legend_hspace"]
    title_fs = args["plot_title_fontsize"]
    ylabel_fs = args["plot_ylabel_fontsize"]
    chrname_fs = args["plot_chrname_fontsize"]
    cnp_label_fs = args["plot_cnp_label_fontsize"]
    cnp_ab_fs = args["plot_cnp_ab_fontsize"]

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"

    ignore_gap = not args["keep_gap"]
    dpi = args["dpi"]
    transparent = args["transparent"]

    tail_alpha = args["tail_alpha"]
    center_alpha = args["center_alpha"]
    onetail_area = args["onetail_area"]

    # figure axis limits
    maxlim_fcn = args["maxlim_fcn"]

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
    # Per-sample setup: compute everything needed for FCN/BAF rows and 2D scatters.
    from matplotlib.backends.backend_pdf import PdfPages

    sns.set_style("whitegrid")

    per_sample = {}
    for sample in samples:
        bin_info: pd.DataFrame = bbcs[bbcs["SAMPLE"] == sample].reset_index(drop=True)
        seg_info: pd.DataFrame = segs[segs["SAMPLE"] == sample].reset_index(drop=True)
        cnp_ids = bin_info["CNP"].to_numpy()
        clone_states = bin_info["CNP"].unique().tolist()
        clone_props = bin_info.iloc[0][[f"u_{clone}" for clone in clones]].to_numpy()
        palette = dict(zip(clone_states, set_palette(num_colors=len(clone_states))))

        bin_info["FCN"] = bin_info["RD"] * gammas[sample]
        tumor_purity = round(np.sum(clone_props[1:]), 3)
        tumor_ploidy = round(
            compute_tumor_ploidy(seg_info, clones, np.sum(clone_props[1:])), 3
        )
        logging.info(f"{sample}: purity={tumor_purity}, ploidy={tumor_ploidy}")

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
                f"{sample}: {num_exceeded} bins have FCN > maxlim_fcn={maxlim_fcn}"
            )
        lim_fcn = (0, min(max(3, max_fcn), maxlim_fcn))

        alphas = get_transparency(
            bin_info,
            by="CNP",
            one_tail=onetail_area,
            tail_alpha=tail_alpha,
            nontail_alpha=center_alpha,
        ).to_numpy()

        per_sample[sample] = {
            "bin_info": bin_info,
            "seg_info": seg_info,
            "cnp_ids": cnp_ids,
            "clone_props": clone_props,
            "palette": palette,
            "exp_bafs": exp_bafs,
            "exp_fcns": exp_fcns,
            "exp_labels": exp_labels,
            "alphas": alphas,
            "lim_fcn": lim_fcn,
            "lim_baf": lim_baf,
            "title": f"{sample}; purity={tumor_purity}; ploidy={tumor_ploidy}",
        }

    patient_id = args["patient_id"] or "panel"
    out_1d = os.path.join(plot_dir, f"{patient_id}{solID and '.' + solID}.1D.pdf")
    out_2d = os.path.join(plot_dir, f"{patient_id}{solID and '.' + solID}.2D.pdf")

    ##################################################
    # 2D scatter: one multi-page PDF, one page per sample. Capture g0_colors
    # for the 1D expected-value overlay so we don't re-run plot_2d below.
    logging.info(f"writing combined 2D PDF: {out_2d}")
    sample_g0 = {}
    with PdfPages(out_2d) as pdf:
        for sample in samples:
            d = per_sample[sample]
            bi = d["bin_info"]
            fig_2d, g0c = plot_2d(
                sample,
                bi,
                bi["BAF"].to_numpy(),
                bi["FCN"].to_numpy(),
                d["exp_bafs"],
                d["exp_fcns"],
                d["exp_labels"],
                d["clone_props"],
                alphas=d["alphas"],
                hue=d["cnp_ids"],
                palette=d["palette"],
                label_clone=True,
                xlab="Minor haplotype B-allele frequency (mhBAF)",
                ylab="Fractional copy number (FCN)",
                xlim=d["lim_baf"],
                ylim=d["lim_fcn"],
                title=d["title"],
                dpi=dpi,
                transparent=transparent,
            )
            sample_g0[sample] = g0c
            pdf.savefig(fig_2d, dpi=dpi, bbox_inches="tight", transparent=transparent)
            plt.close(fig_2d)

    ##################################################
    # Combined 1D PDF: per-sample (FCN, BAF) pairs + shared ASCN CNP + legend.
    # Layout mirrors copytyping plot_rdr_baf_1d_pseudobulk: nested GridSpec
    # gives a tight FCN/BAF pairing within each sample and a larger gap
    # between samples.
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

    k = len(samples)
    cnp_h = max(1.5, n_tumors)
    fig_h = row_height * k + cnp_h * 0.5 + 0.5
    fig_1d = plt.figure(figsize=(row_width, fig_h))
    # Two-section layout: samples_block above, cnp_block below.
    # Inter-sample gap (large) and BAF↔CNP gap (small) are controlled separately.
    outer = GridSpec(
        2,
        1,
        figure=fig_1d,
        height_ratios=[2 * k, cnp_h * 0.5 + 0.3],
        hspace=baf_cnp_hspace,
        top=0.97,
    )
    samples_gs = GridSpecFromSubplotSpec(
        k,
        1,
        subplot_spec=outer[0],
        height_ratios=[2] * k,
        hspace=inter_sample_hspace,
    )
    sample_axes = []
    for si in range(k):
        inner = GridSpecFromSubplotSpec(
            2,
            1,
            subplot_spec=samples_gs[si],
            height_ratios=[1, 1],
            hspace=intra_sample_hspace,
        )
        sample_axes.append((fig_1d.add_subplot(inner[0]), fig_1d.add_subplot(inner[1])))
    inner_bot = GridSpecFromSubplotSpec(
        2, 1, subplot_spec=outer[1], height_ratios=[cnp_h, 1], hspace=cnp_legend_hspace
    )
    ax_cnp = fig_1d.add_subplot(inner_bot[0])
    ax_leg = fig_1d.add_subplot(inner_bot[1])

    for si, sample in enumerate(samples):
        d = per_sample[sample]
        bi = d["bin_info"]
        ax_fcn, ax_baf = sample_axes[si]

        # FCN row: title with sample info, no chrname, no x ticks.
        plot_1d(
            ax_fcn,
            sample,
            bi,
            bi["FCN"],
            regions,
            chrom_sizes,
            exp_colname="exp-FCN",
            exp_groups=d["cnp_ids"],
            val_type="FCN",
            colors=sample_g0[sample],
            hue=d["cnp_ids"],
            palette=d["palette"],
            alphas=d["alphas"],
            ylim=d["lim_fcn"],
            ylab="FCN",
            plot_chrname=False,
            ignore_gap=ignore_gap,
            show_legend=False,
            chr_shift=0,
        )
        ax_fcn.set_xticklabels([])
        ax_fcn.tick_params(axis="x", bottom=False)
        ax_fcn.tick_params(axis="y", left=True, length=4)
        ax_fcn.set_title(d["title"], fontsize=title_fs, fontweight="bold", loc="left")
        ax_fcn.grid(False)
        ax_fcn.yaxis.label.set_fontweight("bold")
        ax_fcn.yaxis.label.set_fontsize(ylabel_fs)
        for spine in ax_fcn.spines.values():
            spine.set_color("black")

        # BAF row: chrname labels shown, no title.
        plot_1d(
            ax_baf,
            sample,
            bi,
            bi["BAF"],
            regions,
            chrom_sizes,
            exp_colname="exp-BAF",
            exp_groups=d["cnp_ids"],
            val_type="BAF",
            colors=sample_g0[sample],
            hue=None,
            palette=None,
            alphas=d["alphas"],
            ylim=d["lim_baf"],
            ylab="mhBAF",
            plot_chrname=True,
            ignore_gap=ignore_gap,
            show_legend=False,
            chr_shift=0,
        )
        ax_baf.set_ylim(-0.05, 1.05)
        ax_baf.grid(False)
        ax_baf.tick_params(axis="x", bottom=True, length=4)
        ax_baf.tick_params(axis="y", left=True, length=4)
        plt.setp(ax_baf.get_xticklabels(), fontweight="bold", fontsize=chrname_fs)
        ax_baf.yaxis.label.set_fontweight("bold")
        ax_baf.yaxis.label.set_fontsize(ylabel_fs)
        for spine in ax_baf.spines.values():
            spine.set_color("black")

    # Shared CNP profile + legend at the bottom (clonal CN, identical across samples)
    _profile_fn = plot_ascn_profile if args["plot_ascn"] else plot_cnv_profile
    _profile_fn(
        ax_cnp,
        per_sample[samples[0]]["seg_info"],
        regions,
        width=row_width,
        height=1,
        plot_chrname=False,
        show_clone_name=True,
        show_prop=False,
    )
    plt.setp(ax_cnp.get_xticklabels(), fontweight="bold", fontsize=cnp_label_fs)
    plt.setp(ax_cnp.get_yticklabels(), fontweight="bold", fontsize=cnp_label_fs)
    plt.setp(ax_cnp.get_yticklabels(minor=True), fontweight="bold", fontsize=cnp_ab_fs)
    if args["plot_ascn"]:
        plot_ascn_legend(
            ax_leg,
            box_w=args["plot_legend_box_w"],
            box_h=args["plot_legend_box_h"],
            tick_len=args["plot_legend_tick_len"],
            label_fontsize=args["plot_legend_label_fontsize"],
        )
    else:
        plot_cnv_legend(ax_leg)

    logging.info(f"writing combined 1D PDF: {out_1d}")
    fig_1d.savefig(out_1d, dpi=dpi, bbox_inches="tight", transparent=transparent)
    plt.close(fig_1d)
    return
