import os
import argparse
import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from cnplot import (
    get_mixcn_cmap,
    get_transparency,
    make_row_spec,
    plot_scatter_1d_multisample,
    plot_scatter_2d,
)

from hatchet.utils import (
    add_file_logging,
    compute_expected_baf_fcn,
    compute_tumor_ploidy,
    normalize_args,
    setup_logging,
)
from hatchet.io_utils import (
    override_solution,
    read_gamma_file,
    read_region_bed,
    read_seg_ucn_file,
)
from hatchet.plot.plot_utils import build_genome_axis, use_editable_fonts


def run(args=None):
    is_cli = isinstance(args, argparse.Namespace)
    args = normalize_args(args)
    if is_cli:
        setup_logging(args)
    logging.info("run hatchet plot-cn")

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

    use_editable_fonts()

    ignore_gap = not args["show_gap"]
    dpi = args["dpi"]
    transparent = args["transparent"]

    tail_alpha = args["tail_alpha"]
    center_alpha = args["center_alpha"]
    onetail_area = args["onetail_area"]

    maxlim_fcn = args["maxlim_fcn"]

    # tumor clones below this per-sample prop are hidden from legends/labels
    display_min_clone_prop = args["min_prop"]

    ##################################################
    # load files
    segs, clones = read_seg_ucn_file(seg_ucn)
    n_clones = len(clones)
    n_tumors = n_clones - 1
    bbcs = read_seg_ucn_file(bbc_ucn)[0]

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
        regions = read_region_bed(region_bed)
        bbcs, segs, n_clones, n_tumors, solID = override_solution(
            bbcs, samples, clusters, n_clones, solfile, regions
        )
    assert n_tumors > 0, "at least one tumor clone must present"

    # Joint per-clone copy-number state, built from the cn_<clone> columns
    # (normal first) and used as the scatter hue.
    clones = ["normal"] + [f"clone{i}" for i in range(1, n_clones)]
    cn_cols = [f"cn_{c}" for c in clones]
    bbcs["cnp"] = bbcs[cn_cols].astype(str).apply(";".join, axis=1)

    is_diploid = args["ploidy"] == "diploid"
    gammas = read_gamma_file(gamma_file, is_diploid)
    logging.info(f"Gammas: {gammas}")
    assert len(gammas) == len(samples), "gamma-file doesn't match samples from BBC file"

    # One genome axis shared by every 1D row, panel, and 2D landmark map,
    # restricted to the chromosomes present in the data.
    genome_axis = build_genome_axis(
        region_bed, genome_size, keep_chroms=chrs, collapse_gaps=ignore_gap
    )

    ##################################################
    # Per-sample domain compute: FCN/FCN-A/FCN-B, per-point alpha, per-state
    # expected values (feed the 1D overlay and the 2D landmarks).
    obs_parts = []
    exp_parts = []
    land_rows = []
    titles = {}
    for sample in samples:
        bin_info: pd.DataFrame = bbcs[bbcs["SAMPLE"] == sample].reset_index(drop=True)
        seg_info: pd.DataFrame = segs[segs["SAMPLE"] == sample].reset_index(drop=True)
        clone_states = bin_info["cnp"].unique().tolist()
        clone_props = bin_info.iloc[0][[f"u_{clone}" for clone in clones]].to_numpy()

        bin_info["FCN"] = bin_info["RD"] * gammas[sample]
        # Allele-specific observed FCN: BAF splits total FCN into the B and A alleles.
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
        for clone_state in clone_states:
            states = [
                (int(x.split("|")[0]), int(x.split("|")[1]))
                for x in clone_state.split(";")
            ]
            fcn_a, fcn_b, fcn, baf = compute_expected_baf_fcn(states, clone_props)
            mask = bin_info["cnp"] == clone_state
            bin_info.loc[mask, "exp-BAF"] = baf
            bin_info.loc[mask, "exp-FCN"] = fcn
            bin_info.loc[mask, "exp-FCN-A"] = fcn_a
            bin_info.loc[mask, "exp-FCN-B"] = fcn_b
            # 2D landmark: one row per distinct joint state (seg.ucn layout).
            row = {"SAMPLE": sample, "exp_BAF": baf, "exp_FCN": fcn}
            for clone, ab in zip(clones, clone_state.split(";")):
                row[f"cn_{clone}"] = ab
                row[f"u_{clone}"] = clone_props[clones.index(clone)]
            land_rows.append(row)

        alphas = get_transparency(
            bin_info,
            by="cnp",
            cols=("RD", "BAF"),
            one_tail=onetail_area,
            tail_alpha=tail_alpha,
            nontail_alpha=center_alpha,
        )
        bin_info["alpha"] = alphas.to_numpy()

        obs_parts.append(bin_info)

        exp = bin_info[["#CHR", "START", "END"]].copy()
        exp[f"exp_FCN_{sample}"] = bin_info["exp-FCN"].to_numpy()
        exp[f"exp_BAF_{sample}"] = bin_info["exp-BAF"].to_numpy()
        exp[f"exp_FCN-A_{sample}"] = bin_info["exp-FCN-A"].to_numpy()
        exp[f"exp_FCN-B_{sample}"] = bin_info["exp-FCN-B"].to_numpy()
        exp_parts.append(exp)

        titles[sample] = f"{sample}; purity={tumor_purity}; ploidy={tumor_ploidy}"

    obs_df = pd.concat(obs_parts, ignore_index=True)
    exp_1d = exp_parts[0]
    for part in exp_parts[1:]:
        exp_1d = exp_1d.merge(part, on=["#CHR", "START", "END"], how="outer")
    land_2d = pd.DataFrame(land_rows)

    # One global palette so 1D and 2D color every joint state identically.
    palette = get_mixcn_cmap(pd.unique(obs_df["cnp"]).tolist())

    # Global limits (panel-wide, shared across samples).
    lim_baf = (0, 1) if obs_df["BAF"].max() > 0.5 else (0, 0.55)
    max_fcn = int(np.ceil(obs_df["FCN"].max()))
    if max_fcn > maxlim_fcn:
        num_exceeded = int(np.sum(obs_df["FCN"] >= maxlim_fcn))
        logging.warning(f"{num_exceeded} bins have FCN > maxlim_fcn={maxlim_fcn}")
    top_fcn = min(max(2, max_fcn), maxlim_fcn)
    lim_fcn = (-0.05 * top_fcn, top_fcn)

    ##################################################
    patient_id = args["patient_id"] or "panel"
    ext = args["img_type"]
    tag = solID and "." + solID
    out_1d = os.path.join(plot_dir, f"{patient_id}{tag}.1D.{ext}")
    out_1d_ab = os.path.join(plot_dir, f"{patient_id}{tag}.1D.FCN_AB.{ext}")
    out_2d = os.path.join(plot_dir, f"{patient_id}{tag}.2D.{ext}")

    common_1d = dict(
        obs_df=obs_df,
        genome_axis=genome_axis,
        groups=samples,
        group_col="SAMPLE",
        expected_df=exp_1d,
        hue="cnp",
        palette=palette,
        alphas="alpha",
        seg_df=segs,
        titles=titles,
        display_min_clone_prop=display_min_clone_prop,
        row_width=row_width,
        row_height=row_height,
        intra_group_hspace=intra_sample_hspace,
        inter_group_hspace=inter_sample_hspace,
        profile_hspace=baf_cnp_hspace,
        markersize=2.0,
    )

    # (FCN, BAF): total FCN over minor-haplotype BAF.
    logging.info(f"writing combined 1D: {out_1d}")
    fig = plot_scatter_1d_multisample(
        row_specs=[
            make_row_spec("FCN", ylabel="FCN", ylim=lim_fcn, href=2.0),
            make_row_spec("BAF", ylabel="mhBAF", ylim=(-0.05, 1.05), href=0.5),
        ],
        **common_1d,
    )
    fig.savefig(out_1d, dpi=dpi, bbox_inches="tight", transparent=transparent)
    plt.close(fig)

    # (FCN-A, FCN-B): allele-specific observed vs CN-expected fractional copy number.
    logging.info(f"writing combined 1D: {out_1d_ab}")
    fig = plot_scatter_1d_multisample(
        row_specs=[
            make_row_spec("FCN-A", ylabel="FCN-A", ylim=lim_fcn, href=1.0),
            make_row_spec(
                "FCN-B", ylabel="FCN-B", ylim=lim_fcn, href=1.0, reverse_y=True
            ),
        ],
        **common_1d,
    )
    fig.savefig(out_1d_ab, dpi=dpi, bbox_inches="tight", transparent=transparent)
    plt.close(fig)

    ##################################################
    # 2D scatter: one page per sample. For pdf, a single multi-page file;
    # otherwise one file per sample.
    logging.info(f"writing combined 2D: {out_2d}")
    pdf = PdfPages(out_2d) if ext == "pdf" else None
    for sample in samples:
        grid = plot_scatter_2d(
            obs_df,
            xcol="BAF",
            ycol="FCN",
            expected_df=land_2d,
            group=sample,
            group_col="SAMPLE",
            normal="normal",
            hue="cnp",
            palette=palette,
            alphas="alpha",
            xlim=lim_baf,
            ylim=lim_fcn,
            xlabel="Minor haplotype B-allele frequency (mhBAF)",
            ylabel="Fractional copy number (FCN)",
            title=titles[sample],
            refline_x=0.5,
            display_min_clone_prop=display_min_clone_prop,
            markersize=2.0,
        )
        fig_2d = grid.figure
        if pdf is not None:
            pdf.savefig(fig_2d, dpi=dpi, bbox_inches="tight", transparent=transparent)
        else:
            path = (
                out_2d
                if len(samples) == 1
                else os.path.join(plot_dir, f"{patient_id}{tag}.2D.{sample}.{ext}")
            )
            fig_2d.savefig(path, dpi=dpi, bbox_inches="tight", transparent=transparent)
        plt.close(fig_2d)
    if pdf is not None:
        pdf.close()
    return
