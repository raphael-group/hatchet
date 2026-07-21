import os
import argparse
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
    plot_combined_1d,
    _row_def,
)
from hatchet.plot.plot_1d2d import get_transparency, plot_2d
from hatchet.plot.plot_cn_utils import build_cnp_palette


def run(args=None):
    is_cli = isinstance(args, argparse.Namespace)
    args = normalize_args(args)
    if is_cli:
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
    if is_cli:
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

    # tumor clones below this per-sample prop are hidden from legends/labels
    display_min_clone_prop = args["min_prop"]

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
        palette = build_cnp_palette(clone_states, clone_props, display_min_clone_prop)

        bin_info["FCN"] = bin_info["RD"] * gammas[sample]
        # Allele-specific observed FCN: minor mhBAF splits total FCN into B (minor) and A (major).
        bin_info["FCN-B"] = bin_info["FCN"] * bin_info["BAF"]
        bin_info["FCN-A"] = bin_info["FCN"] * (1.0 - bin_info["BAF"])
        tumor_purity = round(np.sum(clone_props[1:]), 3)
        tumor_ploidy = round(
            compute_tumor_ploidy(seg_info, clones, np.sum(clone_props[1:])), 3
        )
        logging.info(f"{sample}: purity={tumor_purity}, ploidy={tumor_ploidy}")

        bin_info["exp-BAF"] = 0.0
        bin_info["exp-FCN"] = 0.0
        bin_info["exp-FCN-A"] = 0.0
        bin_info["exp-FCN-B"] = 0.0
        exp_bafs = np.zeros(len(clone_states), dtype=np.float32)
        exp_fcns = np.zeros(len(clone_states), dtype=np.float32)
        exp_labels = []
        for i, clone_state in enumerate(clone_states):
            states = [
                (int(x.split("|")[0]), int(x.split("|")[1]))
                for x in clone_state.split(";")
            ]
            exp_labels.append(states)
            fcn_a, fcn_b, exp_fcns[i], exp_bafs[i] = get_expected_baf_fcn(
                states, clone_props
            )
            mask = bin_info["CNP"] == clone_state
            bin_info.loc[mask, "exp-BAF"] = exp_bafs[i]
            bin_info.loc[mask, "exp-FCN"] = exp_fcns[i]
            bin_info.loc[mask, "exp-FCN-A"] = fcn_a
            bin_info.loc[mask, "exp-FCN-B"] = fcn_b

        lim_baf = (0, 1) if bin_info["BAF"].max() > 0.5 else (0, 0.55)
        max_fcn = int(np.ceil(bin_info["FCN"].max()))
        if max_fcn > maxlim_fcn:
            num_exceeded = np.sum(bin_info["FCN"] >= maxlim_fcn)
            logging.warning(
                f"{sample}: {num_exceeded} bins have FCN > maxlim_fcn={maxlim_fcn}"
            )
        top_fcn = min(max(2, max_fcn), maxlim_fcn)
        lim_fcn = (-0.05 * top_fcn, top_fcn)

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
    ext = args["img_type"]
    out_1d = os.path.join(plot_dir, f"{patient_id}{solID and '.' + solID}.1D.{ext}")
    out_2d = os.path.join(plot_dir, f"{patient_id}{solID and '.' + solID}.2D.{ext}")

    ##################################################
    # 2D scatter: one page per sample. For pdf, a single multi-page file; for
    # other formats, one file per sample. Capture g0_colors for the 1D
    # expected-value overlay so we don't re-run plot_2d below.
    logging.info(f"writing combined 2D: {out_2d}")
    sample_g0 = {}
    pdf = PdfPages(out_2d) if ext == "pdf" else None
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
            display_min_clone_prop=display_min_clone_prop,
            xlab="Minor haplotype B-allele frequency (mhBAF)",
            ylab="Fractional copy number (FCN)",
            xlim=d["lim_baf"],
            ylim=d["lim_fcn"],
            title=d["title"],
            dpi=dpi,
            transparent=transparent,
        )
        sample_g0[sample] = g0c
        if pdf is not None:
            pdf.savefig(fig_2d, dpi=dpi, bbox_inches="tight", transparent=transparent)
        else:
            path = (
                out_2d
                if len(samples) == 1
                else os.path.join(
                    plot_dir,
                    f"{patient_id}{solID and '.' + solID}.2D.{sample}.{ext}",
                )
            )
            fig_2d.savefig(path, dpi=dpi, bbox_inches="tight", transparent=transparent)
        plt.close(fig_2d)
    if pdf is not None:
        pdf.close()

    ##################################################
    # Combined 1D figures: per-sample two-row pairs + shared CNP + legend.
    # Row styling (title/chrname/prop-legend) is fixed by row position; the
    # value/expected columns are supplied per figure via row_defs.
    style = {
        "row_width": row_width,
        "row_height": row_height,
        "inter_sample_hspace": inter_sample_hspace,
        "intra_sample_hspace": intra_sample_hspace,
        "baf_cnp_hspace": baf_cnp_hspace,
        "cnp_legend_hspace": cnp_legend_hspace,
        "title_fs": title_fs,
        "ylabel_fs": ylabel_fs,
        "chrname_fs": chrname_fs,
        "cnp_label_fs": cnp_label_fs,
        "cnp_ab_fs": cnp_ab_fs,
    }
    common = dict(
        regions=regions,
        chrom_sizes=chrom_sizes,
        n_tumors=n_tumors,
        style=style,
        args=args,
        dpi=dpi,
        transparent=transparent,
        ignore_gap=ignore_gap,
        display_min_clone_prop=display_min_clone_prop,
    )

    # (FCN, BAF): total FCN over minor-haplotype BAF.
    fcn_baf_rows = [
        _row_def("FCN", "exp-FCN", "FCN", "FCN", "lim_fcn", top=True, href=2.0),
        _row_def(
            "BAF",
            "exp-BAF",
            "BAF",
            "mhBAF",
            "lim_baf",
            top=False,
            force_ylim=(-0.05, 1.05),
        ),
    ]
    plot_combined_1d(out_1d, samples, per_sample, sample_g0, fcn_baf_rows, **common)

    # (FCN-A, FCN-B): allele-specific observed vs CN-expected fractional copy number.
    out_1d_ab = os.path.join(
        plot_dir, f"{patient_id}{solID and '.' + solID}.1D.FCN_AB.{ext}"
    )
    fcn_ab_rows = [
        _row_def("FCN-A", "exp-FCN-A", "FCN", "FCN-A", "lim_fcn", top=True, href=1.0),
        _row_def(
            "FCN-B",
            "exp-FCN-B",
            "FCN",
            "FCN-B",
            "lim_fcn",
            top=False,
            href=1.0,
            reverse_y=True,
        ),
    ]
    plot_combined_1d(out_1d_ab, samples, per_sample, sample_g0, fcn_ab_rows, **common)
    return
