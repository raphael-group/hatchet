import logging
import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def _row_def(
    val_col,
    exp_col,
    val_type,
    ylab,
    ylim_key,
    top,
    force_ylim=None,
    href=None,
    reverse_y=False,
):
    """Build a per-row spec for plot_combined_1d.

    The top row shows the sample title and hides chrome labels; the bottom row
    shows chromosome names and the per-sample clone-proportion legend.

    Args:
        val_col: bin_info column with observed values to scatter.
        exp_col: bin_info column with CN-expected values for the step overlay.
        val_type: "FCN" or "BAF" (BAF adds the 0.5 reference line).
        ylab: y-axis label.
        ylim_key: per_sample key giving the (lo, hi) y-limit.
        top: True for the upper row (title, no chrname), False for the lower.
        force_ylim: optional (lo, hi) applied after plot_1d.
        href: optional y-value for a dotted grey reference line.
        reverse_y: invert the y-axis so it counts down from 0 at the top.

    Returns:
        dict row spec consumed by plot_combined_1d.
    """
    return {
        "val_col": val_col,
        "exp_col": exp_col,
        "val_type": val_type,
        "ylab": ylab,
        "ylim_key": ylim_key,
        "plot_chrname": not top,
        "set_title": top,
        "use_hue": top,
        "add_prop_legend": not top,
        "force_ylim": force_ylim,
        "href": href,
        "reverse_y": reverse_y,
    }


def plot_combined_1d(
    out_path,
    samples,
    per_sample,
    sample_g0,
    row_defs,
    regions,
    chrom_sizes,
    n_tumors,
    style,
    args,
    dpi,
    transparent,
    ignore_gap,
    display_min_clone_prop,
):
    """Draw the combined 1D figure: per-sample two-row pairs + CNP + legend.

    Layout mirrors copytyping plot_rdr_baf_1d_pseudobulk: a nested GridSpec
    gives a tight within-sample pairing and a larger gap between samples, with
    a shared clonal-CN profile and legend at the bottom.

    Args:
        out_path: output image path.
        samples: ordered sample names.
        per_sample: per-sample dict (bin_info, cnp_ids, palette, alphas,
            clone_props, seg_info, lim_*, title).
        sample_g0: per-sample scatter face colors reused for point coloring.
        row_defs: two row specs from _row_def (top, bottom).
        regions, chrom_sizes: genome layout inputs for plot_1d.
        n_tumors: number of tumor clones (sets CNP profile height).
        style: dict of figure sizes, hspaces, and font sizes.
        args: normalized args (plot_ascn and legend styling).
        dpi, transparent: savefig options.
        ignore_gap: drop bins outside whitelist regions.
        display_min_clone_prop: hide tumor clones below this prop from the legend.
    """
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    from hatchet.plot.plot_1d2d import plot_1d
    from hatchet.plot.plot_cn_utils import (
        plot_ascn_legend,
        plot_ascn_profile,
        plot_cnv_legend,
        plot_cnv_profile,
    )

    k = len(samples)
    cnp_h = max(1.5, n_tumors)
    fig_h = style["row_height"] * k + cnp_h * 0.5 + 0.5
    fig = plt.figure(figsize=(style["row_width"], fig_h))
    outer = GridSpec(
        2,
        1,
        figure=fig,
        height_ratios=[2 * k, cnp_h * 0.5 + 0.3],
        hspace=style["baf_cnp_hspace"],
        top=0.97,
    )
    samples_gs = GridSpecFromSubplotSpec(
        k,
        1,
        subplot_spec=outer[0],
        height_ratios=[2] * k,
        hspace=style["inter_sample_hspace"],
    )
    sample_axes = []
    for si in range(k):
        inner = GridSpecFromSubplotSpec(
            2,
            1,
            subplot_spec=samples_gs[si],
            height_ratios=[1, 1],
            hspace=style["intra_sample_hspace"],
        )
        sample_axes.append((fig.add_subplot(inner[0]), fig.add_subplot(inner[1])))
    inner_bot = GridSpecFromSubplotSpec(
        2,
        1,
        subplot_spec=outer[1],
        height_ratios=[cnp_h, 1],
        hspace=style["cnp_legend_hspace"],
    )
    ax_cnp = fig.add_subplot(inner_bot[0])
    ax_leg = fig.add_subplot(inner_bot[1])

    for si, sample in enumerate(samples):
        d = per_sample[sample]
        bi = d["bin_info"]
        for ridx, rd in enumerate(row_defs):
            ax = sample_axes[si][ridx]
            plot_1d(
                ax,
                sample,
                bi,
                bi[rd["val_col"]],
                regions,
                chrom_sizes,
                exp_colname=rd["exp_col"],
                exp_groups=d["cnp_ids"],
                val_type=rd["val_type"],
                colors=sample_g0[sample],
                hue=d["cnp_ids"] if rd["use_hue"] else None,
                palette=d["palette"] if rd["use_hue"] else None,
                alphas=d["alphas"],
                ylim=d[rd["ylim_key"]],
                ylab=rd["ylab"],
                plot_chrname=rd["plot_chrname"],
                ignore_gap=ignore_gap,
                show_legend=False,
                chr_shift=0,
            )
            if rd["force_ylim"] is not None:
                ax.set_ylim(*rd["force_ylim"])
            if rd["href"] is not None:
                ax.axhline(
                    rd["href"], color="grey", linestyle=":", linewidth=1, zorder=0
                )
            if rd["reverse_y"]:
                ax.invert_yaxis()
            ax.grid(False)
            if rd["plot_chrname"]:
                ax.tick_params(axis="x", bottom=True, length=4)
                plt.setp(
                    ax.get_xticklabels(),
                    fontweight="bold",
                    fontsize=style["chrname_fs"],
                )
            else:
                ax.set_xticklabels([])
                ax.tick_params(axis="x", bottom=False)
            ax.tick_params(axis="y", left=True, length=4)
            ax.yaxis.label.set_fontweight("bold")
            ax.yaxis.label.set_fontsize(style["ylabel_fs"])
            for spine in ax.spines.values():
                spine.set_color("black")
            if rd["set_title"]:
                ax.set_title(
                    d["title"],
                    fontsize=style["title_fs"],
                    fontweight="bold",
                    loc="left",
                )
            if rd["add_prop_legend"]:
                prop_handles = [
                    Line2D(
                        [0],
                        [0],
                        alpha=0,
                        label=(f"Normal: {p:.3f}" if i == 0 else f"Clone {i}: {p:.3f}"),
                    )
                    for i, p in enumerate(d["clone_props"])
                    if i == 0 or p >= display_min_clone_prop
                ]
                ax.legend(
                    handles=prop_handles,
                    loc="center left",
                    bbox_to_anchor=(1.01, 0.5),
                    fontsize="small",
                    fancybox=True,
                    framealpha=0.7,
                    handlelength=0,
                    handletextpad=0,
                )

    # Shared CNP profile + legend at the bottom (clonal CN, identical across samples)
    _profile_fn = plot_ascn_profile if args["plot_ascn"] else plot_cnv_profile
    _profile_fn(
        ax_cnp,
        per_sample[samples[0]]["seg_info"],
        regions,
        width=style["row_width"],
        height=1,
        plot_chrname=True,
        show_clone_name=True,
        show_prop=False,
    )
    plt.setp(
        ax_cnp.get_xticklabels(), fontweight="bold", fontsize=style["cnp_label_fs"]
    )
    plt.setp(
        ax_cnp.get_yticklabels(), fontweight="bold", fontsize=style["cnp_label_fs"]
    )
    plt.setp(
        ax_cnp.get_yticklabels(minor=True),
        fontweight="bold",
        fontsize=style["cnp_ab_fs"],
    )
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

    logging.info(f"writing combined 1D: {out_path}")
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", transparent=transparent)
    plt.close(fig)


def load_gammas(gamma_file: str, is_diploid=True):
    gammas = {}
    with open(gamma_file, "r") as fd:
        for line in fd.readlines():
            sample, gamma_dip, gamma_tet = line.strip().split("\t")
            if is_diploid:
                gammas[sample] = float(gamma_dip)
            else:
                gammas[sample] = float(gamma_tet)
        fd.close()
    return gammas


def override_solution(
    bbcs: pd.DataFrame,
    samples: list,
    clusters: list,
    n_clones: int,
    solfile: str,
    regions: pd.DataFrame,
):
    from hatchet.utils import build_seg_from_bbc

    solID = solfile[str.rindex(solfile, "/") + 1 : -len(".tsv")]
    logging.info(f"overwrite BBC fields with solution {solID}!")
    sol = pd.read_table(solfile)
    assert sorted(sol.CLUSTER.unique().tolist()) == clusters
    assert sorted(sol.SAMPLE.unique().tolist()) == samples

    clones = ["normal"] + [f"clone{i}" for i in range(1, n_clones)]
    for clone in clones:
        bbcs.drop(columns=[f"u_{clone}", f"cn_{clone}"], inplace=True)
    bbcs.drop(columns=["CNP", "PROPS"], inplace=True, errors="ignore")

    bbcs = pd.merge(
        left=bbcs,
        right=sol,
        on=["SAMPLE", "CLUSTER"],
        how="left",
        validate="m:1",
        sort=False,
    )

    # Rebuild derived columns
    n_tumors = len([c for c in bbcs.columns.tolist() if str.startswith(c, "cn_clone")])
    n_clones = n_tumors + 1
    clones = ["normal"] + [f"clone{i}" for i in range(1, n_clones)]
    bbcs["CNP"] = bbcs.apply(
        func=lambda r: ";".join(r[f"cn_{c}"] for c in clones), axis=1
    )
    bbcs["PROPS"] = bbcs.apply(
        func=lambda r: ";".join(str(r[f"u_{c}"]) for c in clones), axis=1
    )

    segs = build_seg_from_bbc(bbcs, regions)
    segs["CNP"] = segs.apply(
        func=lambda r: ";".join(r[f"cn_{c}"] for c in clones), axis=1
    )
    segs["PROPS"] = segs.apply(
        func=lambda r: ";".join(str(r[f"u_{c}"]) for c in clones), axis=1
    )

    return bbcs, segs, n_clones, n_tumors, solID


def get_expected_baf_fcn(cns, props):
    assert len(cns) == len(props)
    A = np.array([x[0] for x in cns])
    B = np.array([x[1] for x in cns])
    y_fcn_a = np.sum(A * props)
    y_fcn_b = np.sum(B * props)
    y_fcn = y_fcn_a + y_fcn_b
    y_baf = np.sum(B * props) / y_fcn

    return y_fcn_a, y_fcn_b, y_fcn, y_baf


def set_palette(num_colors=8, style="whitegrid"):
    sns.set_style(style)
    if num_colors > 8:
        palette = sns.color_palette("husl", n_colors=num_colors)
    else:
        palette = sns.color_palette("Set2", n_colors=num_colors)
    sns.set_palette(palette)
    return palette
